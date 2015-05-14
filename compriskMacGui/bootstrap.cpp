//
//  bootstrap.cpp
//  trial_comprisk
//
//  Created by Piet on 11-03-15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "icm.h"


void bootstrap_SMLE(int iter, double A, double B, int N, int n, int K, int ngrid, double grid[], double data[],double x2[], double data2[], int delta[], int delta2[], int **freq2, double **p2,double **jumploc2, int njumps2[], double **F2, int **ind2, double **dens2, double **Fsmooth2, double ***f3, double **Fsmooth)
{
    int i,j,k,n2,iterations,n_Iterations=1000;
    double h1,h2,phi,c;
    
    c=B;
    h1=c*pow(N,-1.0/5);
    h2=fmin(c*pow(N,-1.0/7),9.9);
    
    data_bootstrap(N,data,x2,delta,delta2);
        
    for (k=0;k<=K;k++)
    {
        for (i=1;i<=n;i++)
            freq2[k][i]=0;
    }
    
    j=0;
    
    
    for (i=1;i<=N;i++)
    {
        if (x2[i]>x2[i-1])
        {
            j++;
            data2[j]=x2[i];
            freq2[delta2[i]][j]=1;
        }
        else
        {
            data2[j]=x2[i];
            freq2[delta2[i]][j]++;
        }
    }
    
    n2=j;
    
    
    ICM(N,n2,freq2,K,ind2,F2,n_Iterations,&phi,&iterations);
    
    for (k=1;k<=K;k++)
    {
        j=0;
        for (i=1;i<=n2;i++)
        {
            if (F2[k][i]>F2[k][i-1])
            {
                p2[k][j]=F2[k][i]-F2[k][i-1];
                jumploc2[k][j]=data2[i];
                j++;
            }
        }
        njumps2[k]=j;
    }
    
    
    
    for (k=1;k<=K;k++)
    {
        for (i=1;i<=ngrid;i++)
            dens2[k][i]=dens_estimate(A,B,k,njumps2,jumploc2,p2,h2,grid[i]);
    }
    
    for (k=1;k<=K;k++)
    {
        for (i=1;i<=ngrid;i++)
            Fsmooth2[k][i]=bdf(A,B,k,njumps2,jumploc2,p2,h1,grid[i]);
    }
    
    
    for (i=0;i<=98;i++)
    {
        for (k=1;k<=K;k++)
            //f3[iter][k][i] = Fsmooth2[k][10*(i+1)]-Fsmooth[k][10*(i+1)];
            f3[iter][k][i] = (Fsmooth2[k][10*(i+1)]-Fsmooth[k][10*(i+1)])/sqrt(varF(N,n2,freq2,k,K,F2,delta2,A,B,data2,h1,grid[10*(i+1)]));
    }
}


void bootstrap_hazard(int iter, double A, double B, int N, int n, int K, int ngrid, double grid[], double data[],double x2[], double data2[], int delta[], int delta2[], int **freq2, double **p2,double **jumploc2, int njumps2[], double **F2, int **ind2, double **dens2, double **Fsmooth2, double ***f3, double **hazard)
{
    int i,j,k,n2,iterations,n_Iterations=1000;
    double h1,h2,phi,c;
    
    c=B;
    h1=c*pow(N,-1.0/5);
    h2=fmin(c*pow(N,-1.0/7),9.9);
    
    data_bootstrap(N,data,x2,delta,delta2);
    
    for (k=0;k<=K;k++)
    {
        for (i=1;i<=n;i++)
            freq2[k][i]=0;
    }
    
    j=0;
    
    
    for (i=1;i<=N;i++)
    {
        if (x2[i]>x2[i-1])
        {
            j++;
            data2[j]=x2[i];
            freq2[delta2[i]][j]=1;
        }
        else
        {
            data2[j]=x2[i];
            freq2[delta2[i]][j]++;
        }
    }
    
    n2=j;
    
    
    ICM(N,n2,freq2,K,ind2,F2,n_Iterations,&phi,&iterations);
    
    for (k=1;k<=K;k++)
    {
        j=0;
        for (i=1;i<=n2;i++)
        {
            if (F2[k][i]>F2[k][i-1])
            {
                p2[k][j]=F2[k][i]-F2[k][i-1];
                jumploc2[k][j]=data2[i];
                j++;
            }
        }
        njumps2[k]=j;
    }
    
    
    
    for (k=1;k<=K;k++)
    {
        for (i=1;i<=ngrid;i++)
            dens2[k][i]=dens_estimate(A,B,k,njumps2,jumploc2,p2,h2,grid[i]);
    }
    
    for (k=1;k<=K;k++)
    {
        for (i=1;i<=ngrid;i++)
            Fsmooth2[k][i]=bdf(A,B,k,njumps2,jumploc2,p2,h1,grid[i]);
    }
    
    
    for (i=0;i<=98;i++)
    {
        for (k=1;k<=K;k++)
            f3[iter][k][i] = dens2[k][10*(i+1)]/(1-Fsmooth2[k][10*(i+1)])-hazard[k][10*(i+1)];
    }
    
}

