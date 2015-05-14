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


void data_bootstrap(int n, double x[], double x2[], int delta[], int delta2[])
{
    int	i,j,k;
    SampleTime *obs;
    
    obs = new SampleTime[n+1];
    
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
    
    delete[] obs;
}
