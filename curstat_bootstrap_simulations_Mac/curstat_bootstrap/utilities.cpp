//
//  icm.cpp
//  comprisk
//
//  Created by Piet on 07-02-15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#include <iostream>
#include "icm.h"
#include "math.h"

#define SQR(x) ((x)*(x))


void cumsum(int m, double cs[], int delta[])
{
    int	i;
    
    cs[1]= delta[1];
    for (i=2;i<=m;i++)
        cs[i] = cs[i-1] + delta[i];
}


void convexmin(int n, double cumw[], double cs[], double w[], double y[])
{
    int	i, j, m;
    
    y[1] = cs[1]/w[1];
    for (i=1;i<= n;i++)
    {
        y[i] = (cs[i] - cs[i - 1])/w[i];
        if (y[i-1]>y[i])
        {
            j = i;
            while ((y[j-1] > y[i]) && (j>1))
            {
                j=j-1;
                y[i] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
                for (m=j;m<i;m++)	y[m] = y[i];
            }
        }
    }
    
    for (i=1;i<= n;i++)
    {
        if (y[i]<0)
            y[i]=0;
        if (y[i]>1)
            y[i]=1;
    }
}

double varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u)
{
    int			i;
    double		t1,t2,t3,sum;
    
    
    sum=0;
    
    for (i=1;i<=n;i++)
    {
        t1=(u-t[i])/h;
        t2=(u+t[i]-2*A)/h;
        t3=(2*B-u-t[i])/h;
        
        sum += SQR(K(t1)-K(t2)-K(t3))*(SQR(F[i]-1)*freq[1][i]+SQR(F[i])*freq[0][i]);
    }
    
    sum = sum/(N*N);
    
    return sum;
}



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
    return fmax(0,sum);
}


weight_t weight(double x)
{
    double y1, y2;
    weight_t temp;
    
    
    y1 = -((2048*(-16 + x*(64 + 5*x*(-32 + x*(43 + 7*(-4 + x)*x)))))/
           (pow(1 + x,6)*(5359 + 5*x*(-3550 +x*(4909 + x*(-3620 + x*(1517 + 35*(-10 + x)*x)))))));
    
    temp.alpha =y1;
    
    y2 = 80640*pow(-1 + x,4)/(pow(1 + x,6)* (5359 + 5*x*(-3550 + x*(4909 + x*(-3620 + x*(1517 + 35*(-10 + x)*x))))));
    
    temp.alpha = y1;
    temp.beta = y2;
    
    return temp;
}


double K(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y=(35.0/32)*pow(1-u,3);
    else
        y=0.0;
    
    return y;
}

double dens_estimate(double A, double B,  int njumps, double jumploc[], double p[], double h, double u)
{
    int i;
    double		rho,sum,x;
    weight_t	weights;
    
    sum=0;
    
    if (u>=A+h && u<=B-h)
    {
        for (i=1;i<=njumps;i++)
        {
            x=(u-jumploc[i])/h;
            sum+= K(x)*p[i]/h;
        }
    }
    else
    {
        if (u<A+h)
        {
            rho=(u-A)/h;
            weights=weight(rho);
            for (i=1;i<=njumps;i++)
            {
                x=(u-jumploc[i])/h;
                sum += (K(x)*weights.alpha + x*K(x)*weights.beta)*p[i]/h;
            }
        }
        else
        {
            if (u>B-h)
            {
                rho = (B-u)/h;
                weights=weight(rho);
                for (i=1;i<=njumps;i++)
                {
                    x=(u-jumploc[i])/h;
                    sum += (K(x)*weights.alpha - x*K(x)*weights.beta)*p[i]/h;
                }
            }
        }
    }
    return fmax(0,sum);
}
