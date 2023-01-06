//BIMcurvedLinU contains all necessary header files and variables.
#include"./Headers/BIMcurvedLinUinteractions.h"
#include<mpi/mpi.h>

double phi = 1.6180339887499;   //golden ratio
double dt = 0.1, tFinal = 80.0, dtMin = 0.01;

int nspheroids = 2;

//Error bounds:
double PicardErrorTolerance = 0.0001, dtTolerance = 0.00002;

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

double RKF4Evolve(vector<BIMobjects> & spheroids, ThreeDVector *DelPos, ThreeDVector *DelRot, double delT)
{
    double TE = 0.0;
    
    ThreeDVector kV[6][nspheroids], kOmega[6][nspheroids];
    ThreeDVector delx[nspheroids], dTheta[nspheroids];  

    for (int iS = 0; iS < nspheroids; iS++) 
    {
        delx[iS] = ThreeDVector(0.0, 0.0, 0.0);
        dTheta[iS] = ThreeDVector(0.0, 0.0, 0.0);
    
        kV[0][iS] = spheroids[iS].uRB*delT;
        kOmega[0][iS] = spheroids[iS].omega*delT;
    }

    for (int m = 1; m < 6; m++)
    {
        //get the k's for RKF4:
        for (int i = 0; i < nspheroids; i++)
        {
            if(m==1) 
            {
                delx[i] = kV[m-1][i]*(2.0/9.0);
                dTheta[i] = kOmega[m-1][i]*(2.0/9.0);
            }

            else if(m==2)   
            {
                delx[i] = kV[m-2][i]*(1.0/12.0) + kV[m-1][i]*(1.0/4.0); 
                dTheta[i] = kOmega[m-2][i]*(1.0/12.0) + kOmega[m-1][i]*(1.0/4.0);
            }

            else if(m==3)   
            {
                delx[i] = kV[m-3][i]*(69.0/128.0) - kV[m-2][i]*(243.0/128.0) + kV[m-1][i]*(135.0/64.0);
                dTheta[i] = kOmega[m-3][i]*(69.0/128.0) - kOmega[m-2][i]*(243.0/128.0) + kOmega[m-1][i]*(135.0/64.0);
            }

            else if(m==4)   
            {
                delx[i] = kV[m-4][i]*(-17.0/12.0) + kV[m-3][i]*(27.0/4.0) + kV[m-2][i]*(-27.0/5.0) + kV[m-1][i]*(16.0/15.0);
                dTheta[i] = kOmega[m-4][i]*(-17.0/12.0) + kOmega[m-3][i]*(27.0/4.0) + kOmega[m-2][i]*(-27.0/5.0) + kOmega[m-1][i]*(16.0/15.0);
            }

            else if(m==5)   
            {
                delx[i] = kV[m-5][i]*(65.0/432.0) + kV[m-4][i]*(-5.0/16.0) + kV[m-3][i]*(13.0/16.0) + kV[m-2][i]*(4.0/27.0) + kV[m-1][i]*(5.0/144.0);
                dTheta[i] = kOmega[m-5][i]*(65.0/432.0) + kOmega[m-4][i]*(-5.0/16.0) + kOmega[m-3][i]*(13.0/16.0) + kOmega[m-2][i]*(4.0/27.0) + kOmega[m-1][i]*(5.0/144.0);
            }
        }
        
        setV(spheroids, delx, dTheta); //sets velocity after translating every body by corresponding delx's and rotating by dTheta's.
        for (int i = 0; i < nspheroids; i++)
        {
            //bring every sphere back to their previous configuration:
            spheroids[i].translate(delx[i]*(-1.0));
            spheroids[i].rotate(dTheta[i].normalize()*(-1.0), dTheta[i].norm());

            kV[m][i] = spheroids[i].uRB*delT;
            kOmega[m][i] = spheroids[i].omega*delT;
        }
    }
    
    for (int i = 0; i < nspheroids; i++)
    {
        // think of DelPos as posNext - posCurrent;
        DelPos[i] = ( kV[0][i]*(47.0/450.0) + kV[2][i]*(12.0/25.0) + kV[3][i]*(32.0/225.0) + kV[4][i]*(1.0/30.0) + kV[5][i]*(6.0/25.0) );   //no contribution from kV[1]
        DelRot[i] = ( kOmega[0][i]*(47.0/450.0) + kOmega[2][i]*(12.0/25.0) + kOmega[3][i]*(32.0/225.0) + kOmega[4][i]*(1.0/30.0) + kOmega[5][i]*(6.0/25.0) );   //no contribution from kOmega[1]
    
        TE += pow( (kV[0][i]*(-1.0/150.0) + kV[2][i]*(3.0/100.0) + kV[3][i]*(-16.0/75.0) +kV[4][i]*(-1.0/20.0) + kV[5][i]*(6.0/25.0)).norm(), 2.0);
        TE += pow( (kOmega[0][i]*(-1.0/150.0) + kOmega[2][i]*(3.0/100.0) + kOmega[3][i]*(-16.0/75.0) +kOmega[4][i]*(-1.0/20.0) + kOmega[5][i]*(6.0/25.0)).norm(), 2.0);
    }

    return sqrt(TE);
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
    spheroidTemplate.rotate(ThreeDVector(0.0, 0.0, 1.0), M_PI/4.0);

    vector<BIMobjects> spheroids;
    spheroids.push_back(spheroidTemplate);

    spheroidTemplate.translate(ThreeDVector(3.0, 0.0, 0.0));
    spheroidTemplate.rotate(ThreeDVector(0.0, 0.0, 1.0), -M_PI/2.0);
    spheroids.push_back(spheroidTemplate);

    //initialize orientation vectors of bodies:
    vector<ThreeDVector> dOrient;
    dOrient.push_back(ThreeDVector(cos(M_PI/4.0), sin(M_PI/4.0), 0.0));
    dOrient.push_back(ThreeDVector(cos(M_PI/4.0), -sin(M_PI/4.0), 0.0));
    
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
    ofstream outPos, outVel, outTym;
    if(myRank==0)
    {
        outPos.open("./OutputData/posSpheroidPoz.txt", ios::out);
        outPos.precision(nPrecision);
        outVel.open("./OutputData/velSpheroidPoz.txt", ios::out);
        outVel.precision(nPrecision);
        outTym.open("./OutputData/timePoz.txt", ios::out);
        outTym.precision(nPrecision);
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

    vector<double> tList;
    tList.push_back(0.0);
    if(myRank==0)    outTym<<tList[0]<<endl;

    double startTime = MPI_Wtime(), tCurr = 0.0, dtNext = dt, tolerance = dtTolerance;
    while(tCurr < tFinal)
    {    
        //RKF4Evolve assigns DelPos and DelRot with which you should translate and rotate the bodies:

        double TE = RKF4Evolve(spheroids, DelPos, DelRot, dt);
        dtNext = 0.9*dt*pow(tolerance/TE, 0.2);
        
        while(TE > tolerance)
        {
            if(dtNext<=dtMin)   
            {
                dt = dtMin;
                TE = RKF4Evolve(spheroids, DelPos, DelRot, dt);
                dtNext = 0.9*dt*pow(tolerance/TE, 0.2);
                dtNext = dtNext<=dtMin?dtMin:dtNext;
                break;
            }
             
            dt = dtNext;
            TE = RKF4Evolve(spheroids, DelPos, DelRot, dt);
            dtNext = 0.9*dt*pow(tolerance/TE, 0.2);   
        }
        
        setV(spheroids, DelPos, DelRot); // translate and rotate each body and calculate the velocity of this new configuration.
        
        //update tCurr once correct dt is calculated:
        if(dtNext>2.0)  dtNext = 2.0;
        tCurr += dt;
        if(tCurr + dtNext >= tFinal)  dtNext = tFinal-tCurr;
        
        tList.push_back(tCurr);

        //update next step and print in file:
        for (int i = 0; i < nspheroids; i++)
        {
            ThreeDVector nHat = DelRot[i].normalize();
            double dPhi = DelRot[i].norm();
            dOrient[i] = dOrient[i].rotate(nHat, dPhi);
            //print result to a file:
            if(myRank==0)   
            {
                outPos<<spheroids[i].X0().x[0]<<'\t'<<spheroids[i].X0().x[1]<<'\t'<<spheroids[i].X0().x[2]<<'\t';
                outPos<<dOrient[i].x[0]<<'\t'<<dOrient[i].x[1]<<'\t'<<dOrient[i].x[2]<<'\t';       
                outVel<<spheroids[i].uRB.x[0]<<'\t'<<spheroids[i].uRB.x[1]<<'\t'<<spheroids[i].uRB.x[2]<<'\t';
                outVel<<spheroids[i].omega.x[0]<<'\t'<<spheroids[i].omega.x[1]<<'\t'<<spheroids[i].omega.x[2]<<'\t';
            }
        }
        if(myRank==0) { outPos<<endl; outVel<<endl; outTym<<tList[tList.size()-1]<<endl;  } 
        dt = dtNext;
        if(myRank==0)   cout<<"dtNext: "<<dtNext<<", Time evolved: "<<tCurr<<endl;
    
    }


    double endTime = MPI_Wtime();
    if(myRank==0)   cout<<"time taken: "<<(endTime-startTime)<<endl;
    
    if(myRank==0)    
    {
        outPos.close();
        outVel.close();
        outTym.close();
    }

    MPI_Finalize();

    return 0;
}