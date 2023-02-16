//BIMcurvedLinU contains all necessary header files and variables.
#include"./Headers/BIMcurvedLinUinteractions.h"
#include<mpi/mpi.h>
#include<string>

double phi = 1.6180339887499;   //golden ratio

int nspheroids = 2;

//mpi variables.
int worldSize, myRank, myStartGC, myEndGC;

//Error bounds:
double PicardErrorTolerance = 0.00005;

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
    //Initialize MPI: 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int workloadVCalc[worldSize];
    
    cout.precision(nPrecision);
 
    
    ThreeDVector initPts[] = {ThreeDVector(0, 1, phi), ThreeDVector(0, -1, phi), ThreeDVector(0, 1, -phi), ThreeDVector(0, -1, -phi),
                                ThreeDVector(1, phi, 0), ThreeDVector(-1, phi, 0), ThreeDVector(1, -phi, 0), ThreeDVector(-1, -phi, 0),
                                ThreeDVector(phi, 0, 1), ThreeDVector(-phi, 0, 1), ThreeDVector(phi, 0, -1), ThreeDVector(-phi, 0, -1)}; //vertices of icosahedron.
    
    BIMobjects spheroidTemplate(initPts, 12);    // no. of vertices of icosahedron = 12.
    
    spheroidTemplate.refineMesh(4);

    double a = 1.0, e = 0.9922, b = a*sqrt(1-e*e);
    spheroidTemplate.scale(b, a, a);
    spheroidTemplate.translate(ThreeDVector(0.0, 1.1, 0.0));   //choose COM of two spheroids at origin, for streamplot.
    //spheroidTemplate.rotate(ThreeDVector(0.0, 0.0, 1.0), M_PI/4.0);

    vector<BIMobjects> spheroids;
    spheroids.push_back(spheroidTemplate);

    spheroidTemplate.translate(ThreeDVector(0.0, -1.2 - b, 0.0));
    spheroidTemplate.rotate(ThreeDVector(0.0, 0.0, 1.0), M_PI/2.0);
    spheroids.push_back(spheroidTemplate);

    //initialize orientation vectors of bodies:
    vector<ThreeDVector> dOrient;
    dOrient.push_back(ThreeDVector(cos(0.0*M_PI/4.0), sin(0.0*M_PI/4.0), 0.0));
    dOrient.push_back(ThreeDVector(cos(M_PI/2.0), sin(M_PI/2.0), 0.0));
    
    if(myRank==0)   cout<<"number of elements: "<<spheroids[0].getElementSize()<<endl;
    
    if(myRank==0)   { spheroids[1].storeElemDat("./meshData/sphere.txt"); cout<<"Mesh writing complete"<<endl;}
    
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
    ofstream outPos, outVel, outField;

    setV(spheroids); // set the velocity using current positions and orientations.

    if(myRank==0)
    {
        outPos.open("./OutputData/posSample.txt", ios::out);
        outPos.precision(nPrecision);
        outVel.open("./OutputData/velSample.txt", ios::out);
        outVel.precision(nPrecision);
        outField.open("./OutputData/uField.txt", ios::out);
        outField.precision(nPrecision);
        for (int iObj = 0; iObj < nspheroids; iObj++)
        {
            outPos<<spheroids[iObj].X0().x[0]<<'\t'<<spheroids[iObj].X0().x[1]<<'\t'<<spheroids[iObj].X0().x[2]<<'\t';       
            outPos<<dOrient[iObj].x[0]<<'\t'<<dOrient[iObj].x[1]<<'\t'<<dOrient[iObj].x[2]<<'\t';       
            //write veloicty:
            outVel<<spheroids[iObj].uRB.x[0]<<'\t'<<spheroids[iObj].uRB.x[1]<<'\t'<<spheroids[iObj].uRB.x[2]<<'\t';
            outVel<<spheroids[iObj].omega.x[0]<<'\t'<<spheroids[iObj].omega.x[1]<<'\t'<<spheroids[iObj].omega.x[2]<<'\t';
        }    
        outPos.close();
        outVel.close();
    }
    double startTime = MPI_Wtime();
    
    if(myRank==0)   cout<<"Writing velocity field data"<<endl;
    //create discrete velocity field:
    int nCellsX = 20, nCellsY = 20;   double domainSizeX = 2.0*1.2, domainSizeY = 2.0*0.3, cellSizeX = domainSizeX/(double)nCellsX, cellSizeY = domainSizeY/(double)nCellsY;    
    if(myRank==0)    cout<<"cellSize: "<<cellSizeX<<", "<<cellSizeY<<endl;
    
    if(myRank==0)
    {
        for (int j = 0; j < nCellsY; j++)
        {
            for (int i = 0; i < nCellsX; i++)
            {
            ThreeDVector uTmp(0.0, 0.0, 0.0);   //no background flow.
            ThreeDVector pt(i*cellSizeX - domainSizeX/2.0, j*cellSizeY - domainSizeY/2.0, 0.0);
            //if( (pt-spheroids[0].X0()).norm()<0.9 || (pt-spheroids[1].X0()).norm()<0.9) continue;
            for (int iObj = 0; iObj < nspheroids; iObj++)
            {
                uTmp = uTmp + spheroids[iObj].getNetFlow(pt);
            }
            outField<<pt.x[0]<<'\t'<<pt.x[1]<<'\t'<<uTmp.x[0]<<'\t'<<uTmp.x[1]<<endl;
            }
        }    
        outField.close();
    }
    
    double endTime = MPI_Wtime();
    if(myRank==0)   cout<<"time taken: "<<(endTime-startTime)<<endl;
    
    if(myRank==0)   cout<<"verticalVel "<<spheroids[0].uRB.x[1]<<endl;
    
    MPI_Finalize();

    return 0;
}