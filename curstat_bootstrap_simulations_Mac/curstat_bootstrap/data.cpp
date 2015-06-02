//
//  data.c
//  comp
//
//  Created by Piet on 20-03-14.
//  Copyright (c) 2014 Piet. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <chrono>
#include <random>
#include "icm.h"


int CompareTime(const void *a, const void *b)
{
    if ((*(SampleTime *) a).t < (*(SampleTime *) b).t)
        return -1;
    if ((*(SampleTime *) a).t > (*(SampleTime *) b).t)
        return 1;
    return 0;
}

void curstatgen(int n, double tt[], int delta[])
{
    int	i;
    SampleTime *obs;
    double x;
    
    obs = new SampleTime[n];
    
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    for (i = 1; i <= n;i++)
    {
        x= dis(generator);
        x=-log(1-(1-exp(-2))*x);
        tt[i]=2*dis(generator);
        if (x<=tt[i]) delta[i]=1;
            else delta[i]=0;
        
    }
    
    for (i=0;i<n;i++)
    {
        obs[i].t=tt[i+1];
        obs[i].delta=delta[i+1];
    }
    
    qsort(obs,n,sizeof(SampleTime),CompareTime);
    
    for (i=1;i<=n;i++)
    {
        tt[i]=obs[i-1].t;
        delta[i]=obs[i-1].delta;
    }
    
    delete[] obs;
}



void data_bootstrap(int n, int *m, double x[], double x2[], double data2[], int **freq, int delta[], int delta2[])
{
    int	i,j,k;
    SampleTime *obs;
    
    obs = new SampleTime[n+1];
    
    for (k=1;k<=n;k++)
        freq[0][k]=freq[1][k]=0;
    
    for (k=1;k<=n;k++)
    {
        j=1+rand()%n;
        
        x2[k]=x[j];
        delta2[k]=delta[j];
    }
    
    for (i=0;i<n;i++)
    {
        obs[i].t=x2[i+1];
        obs[i].delta=delta2[i+1];
    }
    
    qsort(obs,n,sizeof(SampleTime),CompareTime);
    
    for (i=1;i<=n;i++)
    {
        x2[i]=obs[i-1].t;
        delta2[i]=obs[i-1].delta;
    }

    
    x2[0]=0;
    
    j=0;
    
    for (i=1;i<=n;i++)
    {
        if (x2[i]>x2[i-1])
        {
            j++;
            data2[j]=x2[i];
            freq[delta2[i]][j]=1;
        }
        else
        {
            data2[j]=x2[i];
            freq[delta2[i]][j]++;
        }
    }
    
    *m=j;


    delete[] obs;
}

void data_binom(int n, int delta[], double F[])
{
    int	i;
    double x;
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis(0,1);
    
    for (i=1;i<=n;i++)
    {
        x= dis(gen);
        
        if (x<=F[i])
            delta[i]=1;
        else
            delta[i]=0;
    }
    
}
