//BIMcurvedLinU contains all necessary header files and variables.
#include"./Headers/BIMcurvedLinUinteractions.h"
#include<mpi/mpi.h>

double phi = 1.6180339887499;   //golden ratio

int nspheroids = 2;

//Error bounds:
double PicardErrorTolerance = 0.001;

//mpi variables.
int worldSize, myRank, myStartGC, myEndGC;

//sets velocity after translating every body by corresponding delx's and rotating by dTheta's.
void setV(vector<BIMobjects> & spheroids, ThreeDVector *posDx, ThreeDVector *dTheta)
{
    //uSNxt holds last iterated value, set up Prb and P1 for next iterations now:
    for (int iObj = 0; iObj < nspheroids; iObj++)
    {
        if(posDx[iObj].norm() != 0.0)    spheroids[iObj].translate(posDx[iObj]);
        if(dTheta[iObj].norm() != 0.0)   spheroids[iObj].rotate(dTheta[iObj].normalize(), dTheta[iObj].norm());

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
        
      //      if(myRank==0 && iter%10==0)    cout<<"iteration no:"<<iter+1<<"; AvgError: "<<totalError/spheroids[0].nCoordFlat<<endl;
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
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    ThreeDVector initPts[] = {ThreeDVector(0, 1, phi), ThreeDVector(0, -1, phi), ThreeDVector(0, 1, -phi), ThreeDVector(0, -1, -phi),
                                ThreeDVector(1, phi, 0), ThreeDVector(-1, phi, 0), ThreeDVector(1, -phi, 0), ThreeDVector(-1, -phi, 0),
                                ThreeDVector(phi, 0, 1), ThreeDVector(-phi, 0, 1), ThreeDVector(phi, 0, -1), ThreeDVector(-phi, 0, -1)}; //vertices of icosahedron.
    
    int workloadVCalc[worldSize];

    BIMobjects spheroidTemplate(initPts, 12);    // no. of vertices of icosahedron = 12.
    
    spheroidTemplate.refineMesh(2);

    spheroidTemplate.scale(0.6, 1.0, 1.0);
    spheroidTemplate.rotate(ThreeDVector(0.0, 0.0, 1.0), M_PI/2.0);

    vector<BIMobjects> spheroids;
    spheroids.push_back(spheroidTemplate);

    spheroidTemplate.translate(ThreeDVector(6.0, 0.0, 0.0));
  //  spheroidTemplate.rotate(ThreeDVector(0.0, 0.0, 1.0), -M_PI/2.0);
    spheroids.push_back(spheroidTemplate);

    //initialize orientation vectors of bodies:
    vector<ThreeDVector> dOrient;
    dOrient.push_back(ThreeDVector(cos(M_PI/2.0), sin(M_PI/2.0), 0.0));
    dOrient.push_back(ThreeDVector(cos(M_PI/2.0), sin(M_PI/2.0), 0.0));
    
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

    //******** write position and velocities of spheroids to a file***********//
    ofstream outPos, outVel;
    if(myRank==0)
    {
        outPos.open("./OutputData/M2bPosPoz.txt", ios::out);
        outPos.precision(nPrecision);
        outVel.open("./OutputData/M2bVelPoz.txt", ios::out);
        outVel.precision(nPrecision);
        for (int iObj = 0; iObj < nspheroids; iObj++)
        {
            outPos<<spheroids[iObj].X0().x[0]<<'\t'<<spheroids[iObj].X0().x[1]<<'\t'<<spheroids[iObj].X0().x[2]<<'\t';       
            outPos<<dOrient[iObj].x[0]<<'\t'<<dOrient[iObj].x[1]<<'\t'<<dOrient[iObj].x[2]<<'\t';       
        }    
        outPos<<endl;
    }
    
    //setup variables for RK45 adaptive time stepping:
    ThreeDVector DelPos[nspheroids], DelRot[nspheroids]; 
    for (int iS = 0; iS < nspheroids; iS++) 
    {
        DelPos[iS] = ThreeDVector(0.0, 0.0, 0.0);
        DelRot[iS] = ThreeDVector(0.0, 0.0, 0.0);
    }
    setV(spheroids, DelPos, DelRot); // translate and rotate each body and calculate the velocity of this new configuration.
    //store velocity:
    if(myRank==0)
    {
        for (int iObj = 0; iObj < nspheroids; iObj++)
        {    
            outVel<<spheroids[iObj].uRB.x[0]<<'\t'<<spheroids[iObj].uRB.x[1]<<'\t'<<spheroids[iObj].uRB.x[2]<<'\t';
            outVel<<spheroids[iObj].omega.x[0]<<'\t'<<spheroids[iObj].omega.x[1]<<'\t'<<spheroids[iObj].omega.x[2]<<'\t';
        }
        outVel<<endl;
    }

    double startTime = MPI_Wtime(), dxStep = 0.1;
        
    int nSteps = floor((4.0-0.1)/dxStep);
    
    DelPos[1] = DelPos[1] - ThreeDVector(dxStep, 0.0, 0.0);
        
    for (int k = 0; k < nSteps; k++)
    {
        if(myRank==0)   cout<<"dx: "<<spheroids[1].X0().x[0] - spheroids[0].X0().x[0]<<endl;
    
        setV(spheroids, DelPos, DelRot); // translate and rotate each body and calculate the velocity of this new configuration.
        
        for (int i = 0; i < nspheroids; i++)
        {
            //print result to a file:
            if(myRank==0)   
            {
                outPos<<spheroids[i].X0().x[0]<<'\t'<<spheroids[i].X0().x[1]<<'\t'<<spheroids[i].X0().x[2]<<'\t';
                outPos<<dOrient[i].x[0]<<'\t'<<dOrient[i].x[1]<<'\t'<<dOrient[i].x[2]<<'\t';       
                outVel<<spheroids[i].uRB.x[0]<<'\t'<<spheroids[i].uRB.x[1]<<'\t'<<spheroids[i].uRB.x[2]<<'\t';
                outVel<<spheroids[i].omega.x[0]<<'\t'<<spheroids[i].omega.x[1]<<'\t'<<spheroids[i].omega.x[2]<<'\t';
            }
        }

        if(myRank==0) { outPos<<endl; outVel<<endl;  } 

    }

    double endTime = MPI_Wtime();
    if(myRank==0)   cout<<"time taken: "<<(endTime-startTime)<<endl;
    
    if(myRank==0)    
    {
        outPos.close();
        outVel.close();
    }

    MPI_Finalize();

    return 0;
}