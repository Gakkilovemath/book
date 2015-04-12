//
//  icm.cpp
//  trial_comprisk
//
//  Created by Piet on 07-02-15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <Rcpp.h>

#define SQR(x) ((x)*(x))

using namespace std;
using namespace Rcpp;

void	ICM(int ndata, int ndata1, int **freq, int K, int **ind, double **F,
            int n_It, double *phi1, int *iteration, double tol);
void 	ICM_iteration(int ndata, int ndata1, int K, int nlast, int first[], int m[], int **freq, int **ind, double **y,
                      double **F, double **F_new, double **cumw, double **cs, double **v, double **w,
                      double **nabla, double *lambda1, double *alpha1);
double 	compute_lambda(int ndata, int ndata1, int K, int **freq, double **F);
void    compute_Fplus(int ndata, int K, double **F);
void    compute_nabla(int ndata, int ndata1, int K, int first[], int m[], int **freq, double **F, double lambda, double **nabla);
double  f_alpha(int ndata, int ndata1, int K, int **freq, double alpha, double **F, double **F_new, double lambda);
double  f_alpha_prime(int ndata, int ndata1, int K, int **freq, double alpha, double **F, double **F_new, double lambda);
int		fenchelviol(int K, int m[], int **ind, double **F, double **nabla, double tol, double *partsum, double *inprod);
double  phi_ICM(int ndata, int K, int nlast, int **freq, double **F);
void    Compute_F(int K, int ndata, int **freq, int first[], double **y, double **F);
void    cumsum(int K, int m[], double **v, double **cs);
void    convexmin(int K, int m[], double **cumw, double **cs, double **w, double **y);
double golden(int ndata, int ndata1, int K, int **freq, double **F, double **F_new, double a1, double b1,
              double (*f)(int,int,int,int**,double,double**,double**,double), double lambda);

double  KK(double x);
double  bdf(double A, double B, int k, int njumps[], double **jumploc, double **p, double u, double h);
double  max(int n, double x[]);
double  min(int n, double x[]);


// [[Rcpp::export]]

List ComputeMLE(DataFrame input)
{
    double  a,b,*data1,*data2,**F,phi,*grid,**p,**jumploc,step,bandwidth,tol=1.0e-10;
    int     i,j,k,K,m,ndata,ndata2,iterations,*delta1,**ind,**freq,*njumps;
    int     n_Iterations=1000, ngrid=1000;
    
    Rcpp::DataFrame DF = Rcpp::DataFrame(input);
    Rcpp::NumericVector data = DF["V1"];
    Rcpp::IntegerVector delta = DF["V2"];
    
    // determine the number of rows of the data frame
    
    ndata = (int)data.size();
    
    // copy the data vector for use of the C++ procedures
    
    data1 = new double[ndata+1];
    delta1 = new int[ndata+1];
    
    for (i=1;i<=ndata;i++)
    {
        data1[i]=(double)data[i-1];
        delta1[i]=(int)delta[i-1];
    }
    
    
    // determine maximum (b) and minimum (a) for interval on which the SMLE is computed
    // The SMLE uses the values of the MLE which is computed first
    // step is the distance between successive points on the grid
       
    a=min(ndata,data1);
    b=max(ndata,data1);
    step=(b-a)/ngrid;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i] = a+i*step;
    
    // Determine the number of risks
    
    j=0;
    
    for (i=1;i<=ndata;i++)
    {
        if (delta1[i]>j)
            j=delta1[i];
    }
    
    K = j;
    

    // create new data vector to take the ties into account
    
    j=0;
    
    // assuming all data are >0 we put data[0]=0
    
    data1[0]=0;
    
    // create new array data2 for forming data without ties
    
    data2 = new double[ndata+1];
    
    for (i=1;i<=ndata;i++)
    {
        if (data1[i]>data1[i-1])
        {
            j++;
            data2[j]=data1[i];
        }
        else
            data2[j]=data1[i];
    }
    
    ndata2=j;
    
    // the matrix (double pointer) freq will contain the frequencies of risks at the observations
    
    freq = new int *[K+1];
    
    for (i=0;i<K+1;i++)
        freq[i]= new int[ndata2+1];
    
    for (k=0;k<K+1;k++)
    {
        for (i=1;i<=ndata2;i++)
            freq[k][i]=0;
    }
    
    j=0;
    for (i=1;i<=ndata;i++)
    {
        if (data1[i]>data1[i-1])
        {
            j++;
            freq[delta1[i]][j]=1;
        }
        else
            freq[delta1[i]][j] +=1;
    }


    
    // reserve memory for the matrix (double pointer) F which will give the MLE
    // ind is an integer matrix which is going to give the points at which the components of the MLE will be evaluated
    
    F= new double *[K+2];
    ind = new int *[K+2];
    
    for (i=0;i<K+2;i++)
    {
        F[i]= new double[ndata2+1];
        ind[i]= new int[ndata2+1];
    }
    
    // initialize the first values of F in order to later build F recursively
    for (i=1;i<=K+1;i++)
        F[i][0]=0;

    
    // Compute the MLE using the  iterative convex minorant algorithm
    ICM(ndata,ndata2,freq,K,ind,F,n_Iterations,&phi,&iterations,tol);
    
    // Determine the number of points of jump of the sum function F_+
    j=0;
    for (i=1;i<=ndata2;i++)
    {
        if (F[K+1][i]>F[K+1][i-1])
            j++;
    }
    
    m=j;
    
    // reserve memory for the size of the jumps of the components (p) and the locations of these jump (jumploc)
    
    p= new double *[K+2];
    jumploc = new double *[K+2];
    
    for (i=0;i<K+2;i++)
    {
        p[i]= new double[m+1];
        jumploc[i]= new double[m+1];
    }
    
    njumps = new int[m+1];

    // form the output matrix which will contain the MLE
    
    NumericMatrix out1 = NumericMatrix(m,K+2);
    
    j=0;
    for (i=1;i<=ndata2;i++)
    {
        if (F[K+1][i]>F[K+1][i-1])
        {
            j++;
            out1(j-1,0)=data2[i];
            for (k=1;k<=K+1;k++)
                out1(j-1,k) = F[k][i];
        }
    }
    
    // determine sizes and locations of jumps of the subdistribution functions
    
    for (k=1;k<=K+1;k++)
    {
        j=0;
        for (i=1;i<=ndata2;i++)
        {
            if (F[k][i]>F[k][i-1])
            {
                j++;
                p[k][j-1]=F[k][i]-F[k][i-1];
                jumploc[k][j-1]=data2[i];
            }
        }
        njumps[k]=j;
    }
    
    // determination of the bandwidth of the SMLE using a simple rule of thumb
    
    bandwidth = (b-a)*pow(ndata,-1.0/5);
    
    // form the output matrix which will contain the SMLE
    
    NumericMatrix out2 = NumericMatrix(ngrid+1,K+2);
    
    // computation of the SMLE
    
    for (i=0;i<=ngrid;i++)
    {
        out2(i,0)=grid[i];
        for (k=1;k<=K+1;k++)
            out2(i,k)=bdf(a,b,k,njumps,jumploc,p,grid[i],bandwidth);
    }
    
    // make the list for the output, containing the MLE and SMLE matrices
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("SMLE")=out2,Rcpp::Named("log_likelihood")=-phi);
    
    // free memory
    for (i = 0; i < K+2; i++)
        delete[] F[i], delete[] ind[i], delete[] p[i], delete[] jumploc[i];
    
    for (i = 0; i < K+1; i++)
        delete[] freq[i];
    
    delete[] F, delete[] ind, delete[] p, delete[] jumploc, delete[] freq, delete[] njumps, delete[] delta1, delete[] data1, delete[] data2;
    
    return out;
}

void	ICM(int N, int n, int **freq, int K, int **ind, double **F,
            int n_It, double *phi1, int *iteration, double tol)
{
    int		*first,nlast,i,j,k,iteration1,*m,*n1,m_total;
    double	**cumw,**cs,**v,**w,**nabla,**F_new;
    double	**y;
    double	sum,alpha,lambda,inprod,partialsum,phi;
    
    y= new double *[K+2];
    F_new=new double *[K+2];
    cumw=new double *[K+2];
    cs=new double *[K+2];
    v=new double *[K+2];
    w=new double *[K+2];
    nabla=new double *[K+2];
    
    first= new int[K+2];
    m = new int[K+2];
    n1 = new int[K+2];
    
    for (i=0;i<K+2;i++)
    {
        y[i]= new double[n+1];
        F_new[i]= new double[n+1];
        cumw[i]= new double[n+1];
        cs[i]= new double[n+1];
        v[i]= new double[n+1];
        w[i]= new double[n+1];
        nabla[i]= new double[n+1];
    }
    
    j=n;
    
    while (freq[0][j]==0)
        j--;
    
    nlast=j;
    
    for (k=1;k<=K;k++)
    {
        j=1;
        while (freq[k][j]==0)
            j++;
        first[k]=j;
    }
    
    
    for (k=0;k<=K;k++)
        m[k]=n1[k]=0;
    
    for (k=1;k<=K;k++)
    {
        j=0;
        for (i=first[k];i<=n;i++)
        {
            if (freq[k][i]>0)
            {
                j++;
                ind[k][j]=i;
                n1[k] += freq[k][i];
            }
        }
        m[k]=j;
    }
    
    for (k=1;k<=K;k++)
    {
        for (i=0;i<first[k];i++)
            F[k][i]=F_new[k][i]=0;
    }
    
    
    
    m_total=0;
    
    for (k=1;k<=K;k++)
        m_total += n1[k];
    
    
    for (k=1;k<=K;k++)
    {
        for (i=first[k];i<=n;i++)
        {
            if (freq[k][i]>0)
                F[k][i]=F[k][i-1]+freq[k][i]*1.0/(m_total+K);
            else
                F[k][i]=F[k][i-1];
        }
    }
    
    for (i=1;i<=n;i++)
    {
        sum=0;
        for (k=1;k<=K;k++)
            sum += F[k][i];
        F[K+1][i]=sum;
    }
    
    
    for (k=1;k<=K;k++)
    {
        j=0;
        for (i=first[k];i<=n;i++)
        {
            if (freq[k][i]>0)
            {
                j++;
                y[k][j]=F[k][i];
            }
        }
    }
    
    if (nlast<n)
        lambda=compute_lambda(N,n,K,freq,F);
    else
        lambda=0;
    
    for (k=1;k<=K;k++)
    {
        for (i=first[k];i<=n;i++)
            F_new[k][i]=F[k][i];
    }
    
    compute_Fplus(n,K,F_new);
    
    compute_nabla(N,n,K,first,m,freq,F,lambda,nabla);
    
    iteration1=0;
    
    //Rcout << endl << endl;
    
    //Rcout << setw(10) << "iteration" << setw(10) << "alpha" << setw(10) << "F_+" << setw(15) << "phi"
                //<< setw(15) << "lambda" << setw(15) << "inprod" << setw(20) << "partial sum" << endl << endl;
    
    while (iteration1<=n_It && fenchelviol(K,m,ind,y,nabla,tol,&partialsum,&inprod))
    {
        iteration1++;
        
        ICM_iteration(N,n,K,nlast,first,m,freq,ind,y,F,F_new,cumw,cs,v,w,nabla,&lambda,&alpha);
        
        phi=phi_ICM(n,K,nlast,freq,F);
        
        //Rcout << setw(5) << iteration1 << setw(14) << alpha << setw(14) << F[K+1][n] << setprecision(10) <<setw(15)  << phi << setw(10) << lambda << setprecision(10) << setw(20) << inprod  << setprecision(10)  << setw(20) << partialsum << endl;
    }
    
    //Rcout << endl << endl;
    
    *phi1=phi_ICM(n,K,nlast,freq,F);
    *iteration=iteration1;
    
    for (i=0;i<K+2;i++){
        delete[] y[i], delete[] F_new[i], delete[] cumw[i], delete[] cs[i], delete[] v[i], delete[] w[i], delete[] nabla[i];
    }
    
    delete[] y, delete[] cumw, delete[] cs, delete[] v, delete[] w, delete[] nabla;
    delete[] m, delete[] n1, delete[] first;
    
}


void 	ICM_iteration(int N, int n, int K, int nlast, int first[], int m[], int **freq, int **ind, double **y,
                      double **F, double **F_new, double **cumw, double **cs, double **v, double **w,
                      double **nabla, double *lambda1, double *alpha1)
{
    int 	i,j,k;
    double 	a,lambda,alpha;
    
    lambda=*lambda1;
    
    for (k=1;k<=K;k++)
    {
        for (i=1;i<=m[k];i++)
            v[k][i]=w[k][i]=0;
    }
    
    
    for (k=1;k<=K;k++)
    {
        j=0;
        for (i=first[k];i<=n;i++)
        {
            if (freq[k][i]>0)
            {
                j++;
                v[k][j] = freq[k][i]/(N*F[k][i]);
                w[k][j] = freq[k][i]/(N*SQR(F[k][i]));
            }
            if (freq[0][i]>0)
            {
                v[k][j] -= freq[0][i]/(N*(1-F[K+1][i]));
                w[k][j] += freq[0][i]/(N*SQR(1-F[K+1][i]));
            }
        }
        
        if (nlast<n)
            v[k][m[k]] -=lambda;
        
        for (i=1;i<=m[k];i++)
            v[k][i] += w[k][i]*F[k][ind[k][i]];
    }
    
    cumsum(K,m,v,cs);
    
    convexmin(K,m,cumw,cs,w,y);
    
    Compute_F(K,n,freq,first,y,F_new);
    
    if (F_new[K+1][nlast]>=1)
    {
        a=(1-F[K+1][nlast])/(F_new[K+1][nlast]-F[K+1][nlast]);
        for (k=1;k<=K;k++)
        {
            for (i=first[k];i<=n;i++)
                F_new[k][i]=F[k][i]+a*(F_new[k][i]-F[k][i]);
        }
        compute_Fplus(n,K,F_new);
        
        for (k=1;k<=K;k++)
        {
            j=0;
            for (i=first[k];i<=n;i++)
            {
                if (freq[k][i]>0)
                {
                    j++;
                    y[k][j]=F_new[k][i];
                }
            }
        }
    }
    
    if (f_alpha_prime(N,n,K,freq,1.0,F,F_new,lambda)<=0) alpha=1;
    else
        alpha=golden(N,n,K,freq,F,F_new,0.1,1,f_alpha,lambda);
    
    for (k=1;k<=K;k++)
    {
        for (i=first[k];i<=n;i++)
            F[k][i] += alpha*(F_new[k][i]-F[k][i]);
    }
    
    compute_Fplus(n,K,F);
    
    if (nlast<n)
        lambda=compute_lambda(N,n,K,freq,F);
    else
        lambda=0;
    
    compute_nabla(N,n,K,first,m,freq,F,lambda,nabla);
    
    *lambda1=lambda;
    *alpha1=alpha;
}

double 	compute_lambda(int N, int n, int K, int **freq, double **F)
{
    int i;
    double sum=0;
    
    sum=0;
    for (i=1;i<=n;i++)
    {
        if (freq[0][i]>0)
            sum -= freq[0][i]/(N*(1-F[K+1][i]));
    }
    
    sum += 1.0;
    
    if (sum<0)
        sum=0;
    
    return sum;
}

void compute_Fplus(int n, int K, double **F)
{
    int i,k;
    double sum;
    
    for (i=1;i<=n;i++)
    {
        sum=0;
        for (k=1;k<=K;k++)
            sum += F[k][i];
        F[K+1][i]=sum;
    }
}


void compute_nabla(int N, int n, int K, int first[], int m[], int **freq, double **F, double lambda, double **nabla)
{
    int i,j,k;
    
    for (k=1;k<=K;k++)
    {
        for (i=1;i<=m[k];i++)
            nabla[k][i]=0;
    }
    
    for (k=1;k<=K;k++)
    {
        j=0;
        for (i=first[k];i<=n;i++)
        {
            if (freq[k][i]>0)
            {
                j++;
                nabla[k][j] -= freq[k][i]/(N*F[k][i]);
            }
            if (freq[0][i]>0)
                nabla[k][j] += freq[0][i]/(N*(1-F[K+1][i]));
        }
        
        if (freq[0][n]==0)
            nabla[k][m[k]] += lambda;
    }
}

double f_alpha(int N, int n, int K, int **freq, double alpha, double **F, double **F_new, double lambda)
{
    int i,k;
    double sum,a,tol=1.0e-30;
    
    sum=0;
    
    for (i=1;i<=n;i++)
    {
        for (k=1;k<=K;k++)
        {
            if (freq[k][i]>0)
            {
                a = (1-alpha)*F[k][i]+alpha*F_new[k][i];
                if (a>tol)
                    sum -= freq[k][i]*log(a);
                else
                    sum -= freq[k][i]*log(tol);
            }
        }
        if (freq[0][i]>0)
        {
            a = (1-alpha)*(1-F[K+1][i])+alpha*(1-F_new[K+1][i]);
            if (a>tol)
                sum -= freq[0][i]*log(a);
            else
                sum -= freq[0][i]*log(tol);
        }
    }
    if (freq[0][n]==0)
        return sum+N*lambda*alpha*(F_new[K+1][n]-F[K+1][n]);
    else
        return sum;
}

double f_alpha_prime(int N, int n, int K, int **freq, double alpha, double **F, double **F_new, double lambda)
{
    int		i,k;
    double	a,sum,tol=1.0e-30;
    
    sum=0;
    
    
    for (i=1;i<=n;i++)
    {
        for (k=1;k<=K;k++)
        {
            if (freq[k][i]>0)
            {
                a = (1-alpha)*F[k][i]+alpha*F_new[k][i];
                if (a>tol)
                    sum -= freq[k][i]*(F_new[k][i]-F[k][i])/a;
                else
                    sum -= freq[k][i]*(F_new[k][i]-F[k][i])/tol;
            }
        }
        if (freq[0][i]>0)
        {
            a = (1-alpha)*(1-F[K+1][i])+alpha*(1-F_new[K+1][i]);
            if (a>tol)
                sum += freq[0][i]*(F_new[K+1][i]-F[K+1][i])/a;
            else
                sum += freq[0][i]*(F_new[K+1][i]-F[K+1][i])/tol;
        }
    }
    if (freq[0][n]==0)
        return sum + N*lambda*(F_new[K+1][n]-F[K+1][n]);
    else
        return sum;
}


int		fenchelviol(int K, int m[], int **ind, double **y, double **nabla, double tol, double *partsum, double *inprod)
{
    double	sum,sum2;
    int	i,k;
    int	fenchelvioltemp=0;
    
    sum2=0;
    
    for (k=1;k<=K;k++)
    {
        sum=0;
        for (i=1;i<=m[k];i++)
        {
            sum -= nabla[k][i];
            if (sum<sum2) sum2 = sum;
        }
    }
    
    sum=0;
    
    for (k=1;k<=K;k++)
    {
        for (i=1;i<=m[k];i++)
            //sum += nabla[k][i]*F[k][ind[k][i]];
            sum += nabla[k][i]*y[k][i];
    }
    
    sum=fabs(sum);
    
    *inprod = sum;
    *partsum = sum2;
    
    if (sum > tol || sum2 < -tol ) fenchelvioltemp = 1;
    
    return fenchelvioltemp;
}


double phi_ICM(int n, int K, int nlast, int **freq, double **F)
{
    int i,k;
    double sum,tol=1.0e-30;
    
    sum=0;
    
    for (i=1;i<=n;i++)
    {
        for (k=1;k<=K;k++)
        {
            if (freq[k][i]>0)
            {
                if (F[k][i]>tol)
                    sum -= freq[k][i]*log(F[k][i]);
                else
                    return 1.0/tol;
            }
        }
        if (freq[0][i]>0)
        {
            if (1-F[K+1][i]>tol)
                sum -= freq[0][i]*log(1-F[K+1][i]);
            else
                return 1.0/tol;
        }
    }
    
    return sum;    
}


void Compute_F(int K, int n, int **freq, int first[], double **y, double **F)
{
    int k,i,j;
    
    for (k=1;k<=K;k++)
    {
        j=0;
        for (i=first[k];i<=n;i++)
        {
            if (freq[k][i]>0)
            {
                j++;
                F[k][i]=y[k][j];
            }
            else
                F[k][i]=F[k][i-1];
        }
    }
    
    compute_Fplus(n,K,F);
}


void cumsum(int K, int m[], double **v, double **cs)
{
    int	i,k;
    
    for (k=1;k<=K;k++)
    {
        cs[k][1]= v[k][1];
        for (i=2;i<=m[k];i++)
            cs[k][i] = cs[k][i-1] + v[k][i];
    }
}


void convexmin(int K, int m[], double **cumw, double **cs, double **w, double **y)
{
    int	i,j,k,n,m1;
    
    for (k=1;k<=K;k++)
    {
        m1=m[k];
        cs[k][0] = 0;
        cumw[k][0] = 0;
        
        for (i = 1; i<= m1;i++)
            cumw[k][i] = cumw[k][i-1] + w[k][i];
        
        y[k][1] = cs[k][1]/w[k][1];
        
        for (i=2;i<= m1;i++)
        {
            y[k][i] = (cs[k][i] - cs[k][i-1])/w[k][i];
            if (y[k][i-1]>y[k][i])
            {
                j = i;
                while ((y[k][j-1] > y[k][i]) && (j>1))
                {
                    j=j-1;
                    y[k][i] = (cs[k][i]-cs[k][j-1])/(cumw[k][i]-cumw[k][j-1]);
                    for (n=j;n<i;n++)	y[k][n] = y[k][i];
                } 
            }
        }
        
        for (i=1;i<=m1;i++)
        {
            if (y[k][i]<0)
                y[k][i]=0;
            if (y[k][i]>1)
                y[k][i]=1;
        }
    }
}


double golden(int N, int n, int K, int **freq, double **F, double **F_new, double a1, double b1, double (*f)(int,int,int,int**,double,double**,double**,double), double lambda)
{
    double a,b,eps=1.0e-10;
    
    a=a1;
    b=b1;
    
    double k = (sqrt(5.0) - 1.0) / 2;
    double xL = b - k*(b - a);
    double xR = a + k*(b - a);
    
    while (b-a > eps)
    {
        if ((*f)(N,n,K,freq,xL,F,F_new,lambda)<(*f)(N,n,K,freq,xR,F,F_new,lambda))
        {
            b = xR;
            xR = xL;
            xL = b - k*(b - a);
        }
        else
        {
            a = xL;
            xL = xR;
            xR = a + k * (b - a);
        }
    }
    return (a+b)/2;
    
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


double bdf(double A, double B, int k, int njumps[], double **jumploc, double **p, double u, double h)
{
    int			j;
    double		t1,t2,t3,sum;
    
    sum=0;
    
    for (j=0;j<njumps[k];j++)
    {
        t1=(u-jumploc[k][j])/h;
        t2=(u+jumploc[k][j]-2*A)/h;
        t3=(2*B-u-jumploc[k][j])/h;
        sum+= (KK(t1)+KK(t2)-KK(t3))*p[k][j];
    }
    return sum;
}


