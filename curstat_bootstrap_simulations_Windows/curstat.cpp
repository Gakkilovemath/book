//
//  main.cpp
//  trial_comprisk
//
//  Created by Piet on 07-02-15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include "main.h"

double F0(double x);

double F0(double x)
{
    return (1-exp(-x))/(1-exp(-2.0));
}



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


int main(int argc, const char * argv[])
{
    int             i,j,m,m2,n,*delta,*ind,npoints,n_Iterations,NumIt,NumIt2;
    int             *delta2,**freq,**freq2,*ind2,iter,iter2,njumps,ngrid;
    double          *data,*data2,*tt,*F,*F1,*F2,*dens,*dens2,*Fsmooth,*Fsmooth2;
    double          step,*grid,*data0,*data1,*p,*p2,*jumploc,*jumploc2,*hazard;
    double          A,B,c,h1,h2,*lowbound,*upbound;
    double          *cumw,*cs,*w,*y,*y2;
    double          **f3,*f4;
    int             below,above,*percentage;
    clock_t         StartTime, StopTime;
    double          Time_bootstrap;
    
	cout << "Piet Groeneboom 2015" << endl << "For more information please see:" << endl;
	cout << "Nonparametric Estimation under Shape Constraints, Section 9.5," << endl;
	cout << "Piet Groeneboom & Geurt Jongbloed, Cambridge University Press, 2014." << endl << endl;
	
	n_Iterations=1000;
	NumIt=1000;
    NumIt2=1000;
    ngrid=1000;
    n=1000;

	npoints = 100;
    
    below=(int)(0.025*NumIt2-1);
    above=(int)(0.975*NumIt2-1);
    
    njumps=0;
    
    delta = new int[n+1];
    data = new double[n+1];
    data0 = new double[n+1];
    data1 = new double[n+1];
    delta2 = new int[n+1];
    data2 = new double[n+1];
    
    delta[0]=0;
    delta2[0]=0;
    data[0]=data1[0]=data2[0]=0;
    
    
    curstatgen(n,data,delta);
 
    F= new double[n+1];
    F1= new double[n+1];
    F2= new double[n+1];
    tt=new double[n+1];
    ind= new int[n+1];
    ind2= new int[n+1];
    cumw= new double[n+1];
    cs= new double[n+1];
    w= new double[n+1];
    y= new double[n+1];
    y2= new double[n+1];
    
    Fsmooth= new double[ngrid+1];
    Fsmooth2= new double[ngrid+1];
    dens= new double[ngrid+1];
    dens2= new double[ngrid+1];
    p= new double[n+1];
    p2= new double[n+1];
    jumploc= new double[n+1];
    jumploc2= new double[n+1];
    hazard= new double[ngrid+1];
    
    freq = new int*[2];
    for (i=0;i<=1;i++)
        freq[i] = new int[n+1];
    
    freq2 = new int*[1];
    for (i=0;i<=1;i++)
        freq2[i] = new int[n+1];
    
    f4= new double[NumIt2];
    
    f3  = new double*[NumIt2+1];
    
    for (iter2=0;iter2<NumIt2+1;iter2++)
        f3[iter2] = new double[npoints];
    
    lowbound=new double[npoints];
    upbound=new double[npoints];
    percentage = new int[npoints];
    
    for (i=0;i<npoints;i++)
        percentage[i]=0;


    F[0]=F1[0]=F2[0]=0;
    cumw[0]=cs[0]=0;
    
    y[0]=y2[0]=0;
    
    for (i=1;i<=n;i++)
        w[i]=1.0;
    
    A=0;
    B=2.0;
    step=B/ngrid;
    
    c=B;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]=i*step;
    
    StartTime = clock();
    
    for (iter=1;iter<=NumIt;iter++)
    {
        curstatgen(n,data,delta);
        
        for (i=1;i<=n;i++)
        {
            freq[1][i]=delta[i];
            freq[0][i]=1-delta[i];
        }
        
        for (i=1;i<=n;i++)
        {
            w[i]=1.0;
            cumw[i]=i*1.0;
        }
        
        cumsum(n,cs,delta);
        convexmin(n,cumw,cs,w,y);
        
        j=0;
        
        for (i=1;i<=n;i++)
        {
            if (y[i]>y[i-1])
            {
                j++;
                p[j]=y[i]-y[i-1];
                F[j]=y[i];
                data0[j]=data[i];
            }
        }
        
        njumps=j;
        
        // bandwidth for SMLE
        h1 = c*pow(n,-1.0/5);
        
        // bandwidth for density
        h2 = c*pow(n,-1.0/7);
        
        for (i=1;i<=n;i++)
            F1[i]=bdf(B,njumps,data0,p,data[i],h1);
        
        for (i=0;i<=ngrid;i++)
            Fsmooth[i]=bdf(B,njumps,data0,p,grid[i],h1);
        
        for (i=0;i<=ngrid;i++)
            dens[i]=dens_estimate(0.0,B,njumps,data0,p,h2,grid[i]);
        
        for (i=0;i<=ngrid;i++)
            hazard[i]=dens[i]/(1-Fsmooth[i]);
        
        for (iter2=1;iter2<=NumIt2;iter2++)
        {
            data_bootstrap(n,&m,data,tt,data1,freq2,delta,delta2);
            
            //data_binom(n,delta2,F1);
            
            cs[0]=cumw[0]=0;
            
            for (i=1;i<=m;i++)
            {
                w[i]=(double)(freq2[0][i]+freq2[1][i]);
                cs[i]=cs[i-1]+(double)freq2[1][i];
                cumw[i]=cumw[i-1]+w[i];
            }
            
            convexmin(m,cumw,cs,w,y2);
            
            /*for (i=1;i<=n;i++)
                cumw[i]=i*1.0;
            
            cs[0]=cumw[0]=0;
            
            cumsum(n,cs,delta2);
            
            convexmin(n,cumw,cs,w,y2);
            
            m=n;
            
            for (i=1;i<=n;i++)
            {
                data1[i]=data[i];
                freq2[0][i]=1-delta2[i];
                freq2[1][i]=delta2[i];
            }*/
            
            
            j=0;
            
            for (i=1;i<=m;i++)
            {
                if (y2[i]>y2[i-1])
                {
                    j++;
                    data2[j]=data1[i];
                    p2[j]=y2[i]-y2[i-1];
                }
            }
            
            m2=j;
            
            
            for (i=0;i<=ngrid;i++)
                Fsmooth2[i]=bdf(B,m2,data2,p2,grid[i],h1);
            
            for (i=1;i<=99;i++)
                f3[iter2][i]=(Fsmooth2[10*i]-Fsmooth[10*i])/sqrt(varF(n,m,freq2,y2,0.0,B,data1,h1,grid[10*i]));
            
            //for (i=1;i<=99;i++)
                //f3[iter2][i]=Fsmooth2[10*i]-Fsmooth[10*i];
        }
        
        for (i=1;i<npoints;i++)
        {
            for (iter2=0;iter2<NumIt2;iter2++)
                f4[iter2]=f3[iter2+1][i];
            
            qsort(f4,NumIt2-1,sizeof(double),compare);
            
            lowbound[1]= Fsmooth[10]-f4[above]*sqrt(varF(n,n,freq,y,0.0,B,data,h1,grid[10]));
            upbound[1]= Fsmooth[10]-f4[below]*sqrt(varF(n,n,freq,y,0.0,B,data,h1,grid[10]));
            
            if (i>1)
            {
                lowbound[i]= fmax(lowbound[i-1],Fsmooth[10*i]-f4[above]*sqrt(varF(n,n,freq,y,0.0,B,data,h1,grid[10*i])));
                upbound[i]= fmax(upbound[i-1],Fsmooth[10*i]-f4[below]*sqrt(varF(n,n,freq,y,0.0,B,data,h1,grid[10*i])));
            }
            
            //lowbound[i]= Fsmooth[10*i]-f4[above];
            //upbound[i]= Fsmooth[10*i]-f4[below];
            
            if (F0(grid[i*10])<lowbound[i] || F0(grid[i*10])>upbound[i])
                percentage[i]++;
        }
        
         printf("%6d   %12.6f   %12.6f    %12.6f    %12.6f    %6d\n",iter,F0(grid[500]),Fsmooth[500],lowbound[50],upbound[50],percentage[50]);
    }
    
    
    StopTime  = clock();
    Time_bootstrap   = (StopTime - StartTime) / (double) CLOCKS_PER_SEC;
    
    printf("\nThe computations took %10.5f seconds\n",Time_bootstrap);
       
    ofstream file0_("MLE.txt");
    
    if (file0_.is_open())
    {
        file0_ << setprecision(10) << setw(20) << grid[0];
        file0_ << setprecision(11) <<  setw(20) << F[0];
        file0_ << "\n";
        for (i=1;i<=njumps;i++)
        {
            file0_ << setprecision(11) << setw(20) << data0[i];
            file0_ << setprecision(11) <<  setw(20) << F[i];
            file0_ << "\n";
        }
        file0_ << setprecision(10) << setw(20) << grid[ngrid];
        file0_ << setprecision(11) <<  setw(20) << F[njumps];
        file0_ << "\n";
        file0_.close();
    }
    
    ofstream file_("CI_SMLE.txt");
    
    if (file_.is_open())
    {
        file_ << setprecision(10) << setw(20) << grid[10];
        file_ << setprecision(11) <<  setw(20) << fmax(0,lowbound[1])
                    << setprecision(11) <<  setw(20) << fmax(0,upbound[1]);
        file_ << "\n";
        for (i=2;i<npoints;i++)
        {
            file_ << setprecision(10) << setw(20) << grid[10*i];
            file_ << setprecision(11) <<  setw(20) << lowbound[i]
                << setprecision(11) <<  setw(20) << upbound[i];
            file_ << "\n";
        }
        file_.close();
    }
    
    
    ofstream file1_("SMLE.txt");
    
    if (file1_.is_open())
    {
        file1_ << "\n";
        for (i=1;i<=ngrid;i++)
        {
            file1_ << setprecision(10) << setw(20) << grid[i];
            file1_ << setprecision(11) <<  setw(20) << Fsmooth[i];
            file1_ << "\n";
        }
        file1_.close();
    }
    
    ofstream file2_("data.txt");
    
    if (file2_.is_open())
    {
        for (i=1;i<=n;i++)
        {
            file2_ << setprecision(11) << setw(20) << data[i];
            file2_ <<  setw(10) << delta[i];
            file2_ << "\n";
        }
        file2_.close();
    }
    
    ofstream file3_("hazard.txt");
    
    if (file3_.is_open())
    {
        for (i=1;i<=ngrid;i++)
        {
            file3_ << setprecision(10) << setw(20) << grid[i];
            file3_ << setprecision(11) <<  setw(20) << hazard[i];
            file3_ << "\n";
        }
        file3_.close();
    }
    
	printf("\nPictures of 95 percent confidence intervals can be made with CI_SMLE.R\n");
    printf("\nFreeing memory\n");
	
    // free memory
    
    delete[] delta, delete[] data, delete[] data0, delete[] data1, delete[] delta2, delete[] data2,
    delete[] F, delete[] F1, delete[] F2, delete[] tt, delete[] ind, delete[] ind2, delete[] cumw,
    delete[] cs, delete[] w, delete[] y, delete[] y2;
    
    delete[] Fsmooth, delete[] Fsmooth2, delete[] dens, delete[] p, delete[] p2,
    delete[] jumploc, delete[] jumploc2, delete[] hazard;
    
    delete[] lowbound, delete[] upbound, delete[] percentage;

    for (i=0;i<=1;i++)
        delete[] freq[i];
    
    delete[] freq;
    
    for (iter2=0;iter2<NumIt2+1;iter2++)
            delete[] f3[iter2];
    
    delete[] f3;

	printf("\nProgram Complete\n");

	getchar();

	return 0;
}




