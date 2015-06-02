//
//  main_comp.cpp
//  comprisk
//
//  Created by Piet Groeneboom on 07-02-15.
//  Copyright (c) 2015 Piet Groeneboom. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <unistd.h>
#include <sys/param.h>
#include <mach-o/dyld.h>
#include <math.h>
#include "icm.h"

double  max(int n, double x[]);
double  min(int n, double x[]);


int compare(const void *a, const void *b)
{
    double x = *(double*)a;
    double y = *(double*)b;
    
    if (x < y)
        return -1;
    if (x > y)
        return 1;
    return 0;
}

int     GetFileSize(FILE *file)
{
    int    cause, ndata;
    double time;
    
    ndata = 0;
    rewind(file);
    while ( fscanf(file,"%lf %d", &time, &cause) != EOF)
        ndata++;
    return ndata;
}

int     GetFileSize2(int K, FILE *file)
{
    int n,state,ch;
    
    state=n=ch=0;
    
    rewind(file);
    
    while ((ch=getc(file)) != EOF)
    {
        if (ch == ' '|| ch=='\t')
            state = 0;
        else
        {
            if (state == 0)
            {
                state = 1;
                n++;
            }
        }
    }
    
    return n/(K+2);
}




int     GetNumberElements_first_line(FILE *file)
{
    int n,state,ch;
    
    state=n=ch=0;
    
    rewind(file);
    
    while ('\n' != (ch=getc(file)))
    {
        if (ch == ' '|| ch=='\t')
            state = 0;
        else
        {
            if (state == 0)
            {
                state = 1;
                n++;
            }
        }
    }
    
    return n;
}


void main_comp(double *A1, double *B1, int N, int n, int K, int ngrid, double data[], double data1[], int delta[], int **freq, double *grid, double **F, double **Fsmooth, double **hazard, const char *appdir)
{
    int             i,j,k,**ind,iterations,n_Iterations=1000;
    int             *njumps;
    double          phi,**fdens;
    double          step,**p,**jumploc;
    double          A,B,c,h;
    SampleTime		*obs;
    
    
    delta[0]=0;
    data[0]=data1[0]=0;
        
    obs = new SampleTime[N];
    
    for (i=0; i<N; i++)
    {
        obs[i].t = data[i+1];
        obs[i].delta = delta[i+1];
    }
    
    qsort(obs,N,sizeof(SampleTime),CompareTime);
    
    for (i=0; i<N; i++)
    {
        data[i+1]= obs[i].t ;
        delta[i+1] = obs[i].delta;
    }
    
    ind = new int *[K+2];
    
    
    for (i=0;i<K+2;i++)
        ind[i]= new int[N+1];
    
    for (k=1;k<=K+1;k++)
        F[k][0]=0;
    
    njumps = new int[K+2];
    fdens= new double *[K+2];
    p= new double *[K+2];
    jumploc= new double *[K+2];
    
    
    for (i=0;i<K+2;i++)
    {
        fdens[i]= new double[ngrid+1];
        p[i]= new double[n+1];
        jumploc[i]= new double[n+1];
    }
    
    
    ICM(N,n,freq,K,ind,F,n_Iterations,&phi,&iterations);
    
    ofstream file2_("MLE.txt");
    
    if (file2_.is_open())
    {
        for (i=1;i<=n;i++)
        {
            if (F[K+1][i]>F[K+1][i-1])
            {
                file2_ << setprecision(11) << setw(20) << data1[i];
                for (k=1;k<=K+1;k++)
                    file2_ << setprecision(11) <<  setw(20) << F[k][i];
                file2_ << "\n";
            }
        }
        file2_.close();
    }

    
    for (k=1;k<=K+1;k++)
    {
        j=0;
        for (i=1;i<=n;i++)
        {
            if (F[k][i]>F[k][i-1])
            {
                p[k][j]=F[k][i]-F[k][i-1];
                jumploc[k][j]=data1[i];
                j++;
            }
        }
        njumps[k]=j;
    }
    
    A=min(N,data);
    B=max(N,data);
    *A1=A;
    *B1=B;
    step=(B-A)/ngrid;
    
    c=B;
    
    for (i=0;i<=ngrid;i++)
        grid[i]=A+i*step;
    
    h = fmin(c*pow(N,-1.0/7),9.9);
    
    
    for (i=0;i<=ngrid;i++)
    {
        for (k=1;k<K+1;k++)
            fdens[k][i]=dens_estimate(A,B,k,njumps,jumploc,p,h,grid[i]);
    }
    
    h = B*pow(N,-1.0/5);
    
    
    for (i=0;i<=ngrid;i++)
    {
        for (k=1;k<=K;k++)
            Fsmooth[k][i]=bdf(A,B,k,njumps,jumploc,p,h,grid[i]);
    }
    
    for (i=0;i<=ngrid;i++)
    {
        for (k=1;k<=K;k++)
            hazard[k][i]=fdens[k][i]/(1-Fsmooth[k][i]);
    }
    
    ofstream file_("SMLE.txt");

    if (file_.is_open())
    {
        for (i=1;i<=ngrid;i++)
        {
            file_ << setprecision(10) << setw(20) << grid[i];
            for (k=1;k<K+1;k++)
                file_ << setprecision(10) <<  setw(20) << Fsmooth[k][i];
            file_ << "\n";
        }
        file_.close();
    }
    
    ofstream file1_("data.txt");
    
    if (file1_.is_open())
    {
        for (i=1;i<=N;i++)
        {
            file1_ << setprecision(11) << setw(20) << data[i];
            file1_ <<  setw(10) << delta[i];
            file1_ << "\n";
        }
        file1_.close();
    }
    
    
    
    ofstream file3_("hazard.txt");
    
    if (file3_.is_open())
    {
        for (i=1;i<=ngrid;i++)
        {
            file3_ << setprecision(10) << setw(20) << grid[i];
            for (k=1;k<K+1;k++)
                file3_ << setprecision(10) <<  setw(20) << hazard[k][i];
            file3_ << "\n";
        }
        file3_.close();
    }
    
    // free memory
    
    for (i = 0; i < K+2; i++)
         delete[] ind[i], delete[] p[i], delete[] jumploc[i], delete[] fdens[i];
    
    delete[] ind, delete[] p, delete[] jumploc, delete[] njumps, delete[] fdens;
    
    delete[] obs;
}

void main_comp_SMLE_bootstrap(int iter, double A, double B, int N, int n, int K, int ngrid, double data[], int delta[], double grid[], double **Fsmooth, double ***f3, const char *appdir)
{
    int             i,k;
    int             *delta2,**freq2,**ind2,*njumps2;
    double          *data2,*x2,**F2,**dens2,**Fsmooth2;
    double          **p2,**jumploc2;
    
    delta2 = new int[N+1];
    data2 = new double[N+1];
    x2 = new double[N+1];
    
    delta2[0]=0;
    data2[0]=0;
    x2[0]=0;

    freq2 = new int *[K+1];
    
    for (i=0;i<=K;i++)
    {
        freq2[i]= new int[N+1];
    }
    
    for (k=0;k<K+1;k++)
    {
        for (i=1;i<=N;i++)
            freq2[k][i]=0;
    }


    F2= new double *[K+2];
    ind2 = new int *[K+2];
        
    
    for (i=0;i<K+2;i++)
    {
        F2[i]= new double[N+1];
        ind2[i]= new int[N+1];
    }
    
    for (k=1;k<=K+1;k++)
        F2[k][0]=0;
    
    Fsmooth2= new double *[K+2];
    njumps2 = new int[K+2];
    dens2= new double *[K+2];
    p2= new double *[K+2];
    jumploc2= new double *[K+2];
    
    for (i=0;i<K+2;i++)
    {
        Fsmooth2[i]= new double[ngrid+1];
        dens2[i]= new double[ngrid+1];
        p2[i]= new double[n+1];
        jumploc2[i]= new double[n+1];
    }
    
    
    for (k=1;k<=K+1;k++)
        Fsmooth2[k][0]=0;
    
    bootstrap_SMLE(iter,A,B,N,n,K,ngrid,grid,data,x2,data2,delta,delta2,freq2,p2,jumploc2,njumps2,F2,ind2,dens2,Fsmooth2,f3,Fsmooth);

	
    // free memory
    
    for (i = 0; i < K+2; i++)
        delete[] F2[i], delete[] ind2[i], delete[] p2[i], delete[] jumploc2[i], delete[] Fsmooth2[i], delete[] dens2[i];
    
    for (i = 0; i < K+1; i++)
        delete[] freq2[i];
    
    delete[] delta2, delete[] data2, delete[] x2;
    
    delete[] F2, delete[] ind2, delete[] p2, delete[] jumploc2, delete[] Fsmooth2, delete[] dens2,
    delete[] njumps2, delete[] freq2;
}

void main_comp_hazard_bootstrap(int iter, double A, double B, int N, int n, int K, int ngrid, double data[], int delta[], double grid[], double **hazard, double ***f3, const char *appdir)
{
    int             i,k;
    int             *delta2,**freq2,**ind2,*njumps2;
    double          *data2,*x2,**F2,**dens2,**Fsmooth2;
    double          **p2,**jumploc2;
    
    delta2 = new int[N+1];
    data2 = new double[N+1];
    x2 = new double[N+1];
    
    delta2[0]=0;
    data2[0]=0;
    x2[0]=0;
    
    freq2 = new int *[K+1];
    
    for (i=0;i<=K;i++)
    {
        freq2[i]= new int[N+1];
    }
    
    for (k=0;k<K+1;k++)
    {
        for (i=1;i<=N;i++)
            freq2[k][i]=0;
    }
    
    
    F2= new double *[K+2];
    ind2 = new int *[K+2];
    
    
    for (i=0;i<K+2;i++)
    {
        F2[i]= new double[N+1];
        ind2[i]= new int[N+1];
    }
    
    for (k=1;k<=K+1;k++)
        F2[k][0]=0;
    
    Fsmooth2= new double *[K+2];
    njumps2 = new int[K+2];
    dens2= new double *[K+2];
    p2= new double *[K+2];
    jumploc2= new double *[K+2];
    
    for (i=0;i<K+2;i++)
    {
        Fsmooth2[i]= new double[ngrid+1];
        dens2[i]= new double[ngrid+1];
        p2[i]= new double[n+1];
        jumploc2[i]= new double[n+1];
    }
    
    
    for (k=1;k<=K+1;k++)
        Fsmooth2[k][0]=0;
    
    bootstrap_hazard(iter,A,B,N,n,K,ngrid,grid,data,x2,data2,delta,delta2,freq2,p2,jumploc2,njumps2,F2,ind2,dens2,Fsmooth2,f3,hazard);
    
    
    // free memory
    
    for (i = 0; i < K+2; i++)
        delete[] F2[i], delete[] ind2[i], delete[] p2[i], delete[] jumploc2[i], delete[] Fsmooth2[i], delete[] dens2[i];
    
    for (i = 0; i < K+1; i++)
        delete[] freq2[i];
    
    delete[] delta2, delete[] data2, delete[] x2;
    
    delete[] F2, delete[] ind2, delete[] p2, delete[] jumploc2, delete[] Fsmooth2, delete[] dens2,
    delete[] njumps2, delete[] freq2;
}



void confidence_intervals_SMLE(double A, double B, int N, int n, int **freq, int npoints, int K, int NumIt, double grid[], double **F, double data[], int delta[], double **Fsmooth, double **lowbound, double **upbound, double ***f3)
{
    int iter,i,k,below,above;
    double *f4,h;
    
    h= B*pow(N,-1.0/5);
    
    f4 = new double[NumIt];
    
    below=(int)(0.025*NumIt-1);
    above=(int)(0.975*NumIt-1);
    
    for (i=1;i<npoints;i++)
    {
        for (k=1;k<=K;k++)
        {
            for (iter=0;iter<NumIt;iter++)
                f4[iter]=f3[iter][k][i-1];
            
            qsort(f4,NumIt,sizeof(double),compare);
            
            lowbound[k][i]= Fsmooth[k][10*i]-f4[above]*sqrt(varF(N,n,freq,k,K,F,delta,A,B,data,h,grid[10*i]));
            upbound[k][i]= Fsmooth[k][10*i]-f4[below]*sqrt(varF(N,n,freq,k,K,F,delta,A,B,data,h,grid[10*i]));
            
            //lowbound[k][i]= Fsmooth[k][10*i]-f4[above];
            //upbound[k][i]= Fsmooth[k][10*i]-f4[below];
             
            if (i>1)
            {
                lowbound[k][i]=fmax(lowbound[k][i],lowbound[k][i-1]);
                upbound[k][i]=fmax(upbound[k][i],upbound[k][i-1]);
            }
        }
    }
    
    ofstream file_("CI_SMLE.txt");
    
    if (file_.is_open())
    {
        for (i=1;i<npoints;i++)
        {
            file_ << setprecision(10) << setw(20) << grid[10*i];
            for (k=1;k<K+1;k++)
                file_ << setprecision(10) <<  setw(20) << fmax(0,lowbound[k][i])
                << setprecision(10) <<  setw(20) << fmax(0,upbound[k][i]);
            file_ << "\n";
        }
        file_.close();
    }
    
    
    delete[] f4;
}

void confidence_intervals_hazard(int npoints, int K, int NumIt, double grid[], double **hazard,
                          double **lowbound, double **upbound, double ***f3)
{
    int iter,i,k,below,above;
    double *f4;
    
    f4 = new double[NumIt];
    
    below=(int)(0.025*NumIt-1);
    above=(int)(0.975*NumIt-1);
    
    for (i=1;i<npoints;i++)
    {
        for (k=1;k<=K;k++)
        {
            for (iter=0;iter<NumIt;iter++)
                f4[iter]=f3[iter][k][i-1];
            
            qsort(f4,NumIt,sizeof(double),compare);
            
            lowbound[k][i]= hazard[k][10*i]-f4[above];
            upbound[k][i]= hazard[k][10*i]-f4[below];
        }
    }
    
    ofstream file_("CI_hazard.txt");
    
    if (file_.is_open())
    {
        for (i=1;i<npoints;i++)
        {
            file_ << setprecision(10) << setw(20) << grid[10*i];
            for (k=1;k<K+1;k++)
                file_ << setprecision(10) <<  setw(20) << fmax(0,lowbound[k][i])
                << setprecision(10) <<  setw(20) << fmax(0,upbound[k][i]);
            file_ << "\n";
        }
        file_.close();
    }
    
    
    delete[] f4;
}




double min(int n, double x[])
{
    int i;
    double a=x[1];
    
    for (i=2;i<=n;i++)
    {
        if (x[i]<a)
            a=x[i];
    }
    return a;
}

double max(int n, double x[])
{
    int i;
    double a=x[1];
    
    for (i=2;i<=n;i++)
    {
        if (x[i]>a)
            a=x[i];
    }
    return a;
}


