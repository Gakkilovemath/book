//
//  dens.cpp
//  curstat_bootstrap
//
//  Created by Piet on 15/05/15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include "icm.h"

double KK(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y = (16.0 + 35*x - 35*pow(x,3) + 21*pow(x,5) - 5*pow(x,7))/32.0;
    else
    {
        if (x>1)
            y=1;
        else
            y=0;
        
    }
    
    return y;
}

double bdf(double B, int m, double t[], double p[], double u, double h)
{
    int			k;
    double		t1,t2,t3,sum;
    
    
    sum=0;
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        t2=(u+t[k])/h;
        t3=(2*B-u-t[k])/h;
        sum+= (KK(t1)+KK(t2)-KK(t3))*p[k];
        //sum+= KK(t1)*p[k];
    }
    return fmin(1-1.e-10,fmax(sum,1.0e-10));
}
