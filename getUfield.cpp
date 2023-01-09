//BIMcurvedLinU contains all necessary header files and variables.
#include"./Headers/BIMcurvedLinUinteractions.h"
#include<mpi/mpi.h>

double phi = 1.6180339887499;   //golden ratio

int nspheroids = 2;

//mpi variables.
int worldSize, myRank, myStartGC, myEndGC;

//Error bounds:
double PicardErrorTolerance = 0.0001;

//sets velocity after translating every body by corresponding delx's and rotating by dTheta's.
void setV(vector<BIMobjects> & spheroids)
{
    //uSNxt holds last iterated value, set up Prb and P1 for next iterations now:
    for (int iObj = 0; iObj < nspheroids; iObj++)
    {
        spheroids[iObj].refreshuS();
    }

    for (int iter = 0; iter < 250; iter++)
    {
        double totalError = 0.0;
        double error = 0.0;
        for (int iObj = 0; iObj < nspheroids; iObj++)
        {
            spheroids[iObj].resetUsNxt();//all cpus start with zeroes and fill only there workloads.
            for (int iGC = myStartGC; iGC < myEndGC; iGC++)
            {
                error += spheroids[iObj].picardIterate(iGC, spheroids, iObj);
            }
            // get correct uSNxt for all processes. 
            for (int iGC = 0; iGC < spheroids[iObj].nCoordFlat; iGC++)
            {
                double senduSNxt[3] = {spheroids[iObj].uSNxt[iGC].x[0], spheroids[iObj].uSNxt[iGC].x[1], spheroids[iObj].uSNxt[iGC].x[2]};           
                double getuSNxt[3] = {0.0, 0.0, 0.0};
                    
                MPI_Allreduce(&senduSNxt, &getuSNxt, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                spheroids[iObj].uSNxt[iGC].set(getuSNxt[0], getuSNxt[1], getuSNxt[2]);
            }
            
            MPI_Allreduce(&error, &totalError, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            spheroids[iObj].refreshuS();  //use new values as soon as you get them.
        }
        
        if(myRank==0 && iter%10==0)    cout<<"iteration no:"<<iter+1<<"; AvgError: "<<totalError/spheroids[0].nCoordFlat<<endl;
            
        if( totalError/spheroids[0].nCoordFlat <= PicardErrorTolerance )  
        {
            if(myRank==0)    cout<<"iteration no:"<<iter+1<<"; AvgError: "<<totalError/spheroids[0].nCoordFlat<<endl;
            break;
        }
    }
        
    for (int iObj = 0; iObj < nspheroids; iObj++)
    {
        spheroids[iObj].projectRB(); 

        spheroids[iObj].uRB = spheroids[iObj].integrateVectorfunc(&BIMobjects::getURB);
        spheroids[iObj].omega = spheroids[iObj].integrateVectorfunc(&BIMobjects::getOmegaRB);
    }
}

int main(int argc, char **argv)
{
    cout.precision(nPrecision);
    
    //Initialize MPI: 
    int workloadVCalc[worldSize];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    ThreeDVector initPts[] = {ThreeDVector(0, 1, phi), ThreeDVector(0, -1, phi), ThreeDVector(0, 1, -phi), ThreeDVector(0, -1, -phi),
                                ThreeDVector(1, phi, 0), ThreeDVector(-1, phi, 0), ThreeDVector(1, -phi, 0), ThreeDVector(-1, -phi, 0),
                                ThreeDVector(phi, 0, 1), ThreeDVector(-phi, 0, 1), ThreeDVector(phi, 0, -1), ThreeDVector(-phi, 0, -1)}; //vertices of icosahedron.
    
    BIMobjects spheroidTemplate(initPts, 12);    // no. of vertices of icosahedron = 12.
    
    spheroidTemplate.refineMesh(3);

    spheroidTemplate.scale(0.6, 1.0, 1.0);
    spheroidTemplate.translate(ThreeDVector(-1.5, 0.0, 0.0));   //choose COM of two spheroids at origin, for streamplot.
    spheroidTemplate.rotate(ThreeDVector(0.0, 0.0, 1.0), M_PI/4.0);

    vector<BIMobjects> spheroids;
    spheroids.push_back(spheroidTemplate);

    spheroidTemplate.translate(ThreeDVector(3.0, 0.0, 0.0));
    spheroidTemplate.rotate(ThreeDVector(0.0, 0.0, 1.0), -M_PI/2.0);
    spheroids.push_back(spheroidTemplate);

    if(myRank==0)   cout<<"number of elements: "<<spheroids[0].getElementSize()<<endl;
    
   // if(myRank==0)   { spheroids[0].storeElemDat(); }
    
    //determine workloads for each core:
    for (int i = 0; i < worldSize; i++)
    {
        workloadVCalc[i] = spheroids[0].nCoordFlat/worldSize; //nCoordFlat doesnt change for different objects made out of sphere.
        if(i < spheroids[0].nCoordFlat%worldSize) workloadVCalc[i]++; // take care of remainders.
    }

    myStartGC = 0;
    for (int i = 0; i < myRank; i++)
    {
        myStartGC += workloadVCalc[i];
    }
    myEndGC = myStartGC + workloadVCalc[myRank];

    cout<<"myRank:"<<myRank<<" myStartGC:"<<myStartGC<<" myEndGC:"<<myEndGC<<endl;

    //******** write double layer density for each spheroid to a file***********//
    ofstream outField;
    if(myRank==0)
    {
        outField.open("./OutputData/uFieldPoz.txt", ios::out);
        outField.precision(nPrecision);
    }
    
  //  setV(spheroids, DelPos, DelRot); // translate and rotate each body and calculate the velocity of this new configuration.

    double startTime = MPI_Wtime();
    
    setV(spheroids);
    
    //create discrete velocity field:
    int nCells = 100;   double domainSize = 2.0*6.0, cellSize = domainSize/(double)nCells;    
    cout<<"cellSize: "<<cellSize<<endl;

    if(myRank==0)    
    {
        for (int j = 0; j < nCells; j++)
        {
            for (int i = 0; i < nCells; i++)
            {
                ThreeDVector uNet(0.0, 0.0, 0.0);   //no background flow.
                ThreeDVector pt(i*cellSize - domainSize/2.0, j*cellSize - domainSize/2.0, 0.0);
                for (int iObj = 0; iObj < nspheroids; iObj++)
                {
                    uNet = uNet + spheroids[iObj].getNetFlow(pt);
                }
                //store 2D data of uNet in file:
                outField<<pt.x[0]<<'\t'<<pt.x[1]<<'\t'<<uNet.x[0]<<'\t'<<uNet.x[1]<<'\n';
            }
            
        }
        outField.close();
    }

    double endTime = MPI_Wtime();
    if(myRank==0)   cout<<"time taken: "<<(endTime-startTime)<<endl;
    
    MPI_Finalize();

    return 0;
}