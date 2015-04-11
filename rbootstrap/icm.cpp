//  icm.cpp
//  trial_comprisk
//
//  Created by Piet on 07-02-15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <Rcpp.h>


#define SQR(x) ((x)*(x))

using namespace Rcpp;

typedef struct{double alpha, beta;} weight_t;
typedef struct
{
    double t;
    int delta;
} SampleTime;

void	ICM(int N, int n, int **freq, int K, int **ind, double **F,
            int n_It, double *phi1, int *iteration, double tol);
void 	ICM_iteration(int N, int n, int K, int nlast, int first[], int m[], int **freq, int **ind, double **y,
                      double **F, double **F_new, double **cumw, double **cs, double **v, double **w,
                      double **nabla, double *lambda1, double *alpha1);
double 	compute_lambda(int N, int n, int K, int **freq, double **F);
void    compute_Fplus(int N, int K, double **F);
void    compute_nabla(int N, int n, int K, int first[], int m[], int **freq, double **F, double lambda, double **nabla);
double  f_alpha(int N, int n, int K, int **freq, double alpha, double **F, double **F_new, double lambda);
double  f_alpha_prime(int N, int n, int K, int **freq, double alpha, double **F, double **F_new, double lambda);
int		  fenchelviol(int K, int m[], int **ind, double **F, double **nabla, double tol, double *partsum, double *inprod);
double  phi_ICM(int N, int K, int nlast, int **freq, double **F);
void    Compute_F(int K, int N, int **freq, int first[], double **y, double **F);
void    cumsum(int K, int m[], double **v, double **cs);
void    convexmin(int K, int m[], double **cumw, double **cs, double **w, double **y);
double  golden(int N, int n, int K, int **freq, double **F, double **F_new, double a1, double b1,
              double (*f)(int,int,int,int**,double,double**,double**,double), double lambda);
double  max(int n, double x[]);
double  min(int n, double x[]);
int     compare(const void *a, const void *b);
int     CompareTime(const void *a, const void *b);

weight_t weight(double x);
double  K(double x);
double  KK(double x);
double  dens_estimate(double A, double B, int k, int njumps[], double **jumploc, double **p, double h, double u);
double  bdf(double A, double B, int k, int njumps[], double **jumploc, double **p, double h, double u);
double  max(int n, double x[]);
double  min(int n, double x[]);
void    data_bootstrap(int n, double x[], double x2[], int delta[], int delta2[]);
void    bootstrap(double A, double B, int NumIt, int N, int K, int ngrid, double grid[], double data[],
                  double x2[], double data2[], int delta[], int delta2[], int **freq2, double **p2,
                  double **jumploc2, int njumps2[], double **F2, int **ind2, double **dens2, double **Fsmooth2,
                  double ***f3, double **hazard);

// [[Rcpp::export]]

List ComputeMLE(DataFrame input)
{
    double  a,b,*data0,*data1,**F,phi,*grid,**p,**jumploc,step,h1,h2,tol=1.0e-10;
    double  **Fsmooth,**fdens,**hazard,*x2,*data2,***f3,**lowbound,**upbound,*f4;
    int     i,j,k,K,m,N,n,iterations,*delta0,**ind,**freq,*njumps,*delta2;
    int     percentile1,percentile2,iter,n_Iterations=1000,ngrid=1000,NumIt=1000;
    SampleTime *obs;
    
    
    Rcpp::DataFrame DF = Rcpp::DataFrame(input);
    Rcpp::NumericVector data = DF["V1"];
    Rcpp::IntegerVector delta = DF["V2"];
    
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "Piet Groeneboom 2015, academic use only." << std::endl;
    Rcpp::Rcout << "For more information please see:" << std::endl;
    Rcpp::Rcout << "Nonparametric Estimation under Shape Constraints, pp. 10-11 and 367-369," << std::endl;
    Rcpp::Rcout << "Piet Groeneboom & Geurt Jongbloed, Cambridge University Press, 2014." << std::endl << std::endl;
    
    
    // determine the number of rows of the data frame
    
    N = (int)data.size();
    
    Rcpp::Rcout << "Sample size:" << std::setw(10) << N << std::endl;
    
    obs = new SampleTime[N];
    
    for (i=0;i<N;i++)
    {
        obs[i].t= (double)data[i];
        obs[i].delta= (int)delta[i];
    }

    qsort(obs,N,sizeof(SampleTime),CompareTime);
  
    // copy the data vector for use of the C++ procedures
    
    data0 = new double[N+1];
    x2 = new double[N+1];
    data2 = new double[N+1];
    delta0 = new int[N+1];
    delta2 = new int[N+1];
    
    x2[0]=0;
    
    for (i=1;i<=N;i++)
    {
        data0[i]=obs[i-1].t;
        delta0[i]=obs[i-1].delta;
    }
    
    delete[] obs;
    
    // determine maximum (b) and minimum (a) for interval on which the SMLE is computed
    // The SMLE uses the values of the MLE which is computed first
    // step is the distance between successive points on the grid
    
    a=min(N,data0);
    b=max(N,data0);
    
    step=(b-a)/ngrid;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i] = a+i*step;
    
    // Determine the number of risks
    
    j=0;
    
    for (i=1;i<=N;i++)
    {
        if (delta0[i]>j)
            j=delta0[i];
    }
    
    K = j;
    
    // percentile points for 95% bootstrap confidence intervals
    
    percentile1=25;
    percentile2=975;
    
    // assuming all data are >0 we put data[0]=0
    
    data0[0]=0;
    
    // create new data vector to take the ties into account
    
    j=0;
    
    data1 = new double[N+1];
    
    for (i=1;i<=N;i++)
    {
        if (data0[i]>data0[i-1])
        {
            j++;
            data1[j]=data0[i];
        }
        else
            data1[j]=data0[i];
    }
    
    n=j;
    
    Rcpp::Rcout << "Number of unique observations after correction for ties:" << std::setw(5) << n << std::endl << std::endl;
    
    // the matrix (double pointer) freq will contain the frequencies of risks at the observations
    
    freq = new int *[K+1];
    
    for (i=0;i<=K;i++)
        // freq[i]= new int[n+1];
        freq[i]= new int[N+1];
    
    for (k=0;k<=K;k++)
    {
        // for (i=1;i<=n;i++)
        for (i=1;i<=N;i++)
            freq[k][i]=0;
    }
    
    j=0;
    for (i=1;i<=N;i++)
    {
        if (data0[i]>data0[i-1])
        {
            j++;
            freq[delta0[i]][j]=1;
        }
        else
            freq[delta0[i]][j]++;
    }
    
    
    
    // reserve memory for the matrix (double pointer) F which will give the MLE
    // ind is an integer matrix which is going to give the points at which the components of the MLE will be evaluated
    
    F= new double *[K+2];
    ind = new int *[K+2];
    
    for (i=0;i<K+2;i++)
    {
        F[i]= new double[N+1];
        ind[i]= new int[N+1];
    }
    
    // initialize the first values of F in order to later build F recursively
    for (i=1;i<=K+1;i++)
        F[i][0]=0;
    
    
    // Compute the MLE using the  iterative convex minorant algorithm
    ICM(N,n,freq,K,ind,F,n_Iterations,&phi,&iterations,tol);
    
    // Determine the number of points of jump of the sum function F_+
    j=0;
    for (i=1;i<=n;i++)
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
        p[i]= new double[N+1];
        jumploc[i]= new double[N+1];
    }
    
    njumps = new int[N+1];
    
    Rcpp::Rcout << "Determining the output matrix which will contain the MLE" << std::endl;
    
    NumericMatrix out1 = NumericMatrix(m,K+2);
    
    j=0;
    for (i=1;i<=n;i++)
    {
        if (F[K+1][i]>F[K+1][i-1])
        {
            j++;
            out1(j-1,0)=data1[i];
            for (k=1;k<=K+1;k++)
                out1(j-1,k) = F[k][i];
        }
    }
    
    // determine sizes and locations of jumps of the subdistribution functions
    
    for (k=1;k<=K+1;k++)
    {
        j=0;
        for (i=1;i<=n;i++)
        {
            if (F[k][i]>F[k][i-1])
            {
                j++;
                p[k][j-1]=F[k][i]-F[k][i-1];
                jumploc[k][j-1]=data1[i];
            }
        }
        njumps[k]=j;
    }
    
    // determination of the bandwidth of the SMLE and density using a simple rule of thumb
    
    h1 = b*pow(N,-1.0/5);
    h2 = fmin(9.9,b*pow(N,-1.0/6));
    
    // matrices for SMLE, density and hazard
    
    Fsmooth = new double *[K+2];
    fdens = new double *[K+2];
    hazard = new double *[K+2];
    
    for (i=0;i<K+2;i++)
    {
        Fsmooth[i]= new double[ngrid+1];
        fdens[i]= new double[ngrid+1];
        hazard[i]= new double[ngrid+1];
    }
    
    f3  = new double**[NumIt+1];
    
    for (iter=0;iter<=NumIt;iter++)
        f3[iter] = new double *[K+1];
    
    for (iter=0;iter<=NumIt;iter++)
        for (k=0;k<=K;k++)
            f3[iter][k] = new double[100];
    
    
    // Compute SMLE
    
    for (i=0;i<=ngrid;i++)
    {
        for (k=1;k<=K;k++)
        {
            Fsmooth[k][i]=bdf(a,b,k,njumps,jumploc,p,h1,grid[i]);
            fdens[k][i] = dens_estimate(a,b,k,njumps,jumploc,p,h2,grid[i]);
            hazard[k][i] = fdens[k][i]/(1-Fsmooth[k][i]);
        }
    }
    
    Rcpp::Rcout << "Starting bootstrap" << std::endl;
    
    bootstrap(a,b,NumIt,N,K,ngrid,grid,data0,x2,data2,delta0,delta2,freq,p,jumploc,njumps,
              F,ind,fdens,Fsmooth,f3,hazard);
    
    Rcpp::Rcout << std::endl << std::endl;
    
    lowbound = new double *[K+1];
    upbound = new double *[K+1];
    
    for (k=0;k<K+1;k++)
    {
        lowbound[k] = new double[100];
        upbound[k] = new double[100];
    }
    
    Rcpp::Rcout << "Sorting results of bootstrap" << std::endl << std::endl;
    
    f4 = new double[NumIt];
    
    for (i=1;i<=99;i++)
    {
        for (k=1;k<=K;k++)
        {
            for (iter=0;iter<NumIt;iter++)
                f4[iter]=f3[iter][k][i];
            
            qsort(f4,NumIt,sizeof(double),compare);
            
            // lower and upper bounds for bootstrap confidence intervals
            
            lowbound[k][i]= fmax(0,hazard[k][10*i]-f4[percentile2-1]);
            upbound[k][i]= fmax(0,hazard[k][10*i]-f4[percentile1-1]);
        }
    }
    
    Rcpp::Rcout << "Forming matrix with hazards" << std::endl;
    
    NumericMatrix out2 = NumericMatrix(ngrid+1,K+1);
    
    // computation of the SMLE
    
    for (i=0;i<=ngrid;i++)
    {
        out2(i,0)=grid[i];
        for (k=1;k<=K;k++)
            out2(i,k) = hazard[k][i];
    }
    
    Rcpp::Rcout << "Forming matrix with confidence intervals" << std::endl;
    
    NumericMatrix out3 = NumericMatrix(99,2*K+1);
    
    
    for (i=0;i<99;i++)
    {
        out3(i,0)=grid[10*(i+1)];
        for (k=1;k<=K;k++)
        {
            out3(i,2*k-1)=lowbound[k][i+1];
            out3(i,2*k)=upbound[k][i+1];
        }
    }
    
    Rcpp::Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("hazard")=out2,Rcpp::Named("CI")=out3,Rcpp::Named("log_likelihood")=-phi);
    
    //List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("hazard")=out2,Rcpp::Named("log_likelihood")=-phi);
    
    Rcpp::Rcout << "Freeing memory" << std::endl;
    
    // free memory

    for (i = 0; i < K+2; i++)
        delete[] F[i], delete[] ind[i], delete[] p[i], delete[] jumploc[i], delete[] Fsmooth[i], delete[] fdens[i], delete[] hazard[i];
    
    for (i = 0; i < K+1; i++)
        delete[] freq[i], delete[] lowbound[i], delete[] upbound[i];
    
    delete[] F, delete[] ind, delete[] p, delete[] jumploc, delete[] freq, delete[] njumps, delete[] delta0, delete[] delta2, delete[] data0, delete[] data1, delete[] data2, delete[] x2,
    delete[] grid, delete[] Fsmooth, delete[] fdens, delete[] hazard, delete[] f4, delete[] lowbound,
    delete[] upbound;
    
    for (iter = 0;iter < NumIt;iter++)
    {
        for (i =0;i< K+1;i++)
            delete[] f3[iter][i];
    }
    
    for (iter = 0;iter < NumIt;iter++)
        delete[] f3[iter];
    
    delete[] f3;
    
    return out;
}



void	ICM(int N, int n, int **freq, int K, int **ind, double **F,
            int n_It, double *phi1, int *iteration, double tol)
{
    int		*first,nlast,i,j,k,iteration1,*m,*n1,m_total;
    double	**cumw,**cs,**v,**w,**nabla,**F_new;
    double	**y;
    double	alpha,lambda,inprod,partialsum;
    
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
    
    for (k=1;k<=K;k++)
        y[k][0]=0;
    
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
        double sum=0;
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
    
    while (iteration1<=n_It && fenchelviol(K,m,ind,y,nabla,tol,&partialsum,&inprod))
    {
        iteration1++;
        
        ICM_iteration(N,n,K,nlast,first,m,freq,ind,y,F,F_new,cumw,cs,v,w,nabla,&lambda,&alpha);
    }
    
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
    double 	lambda,alpha;
    
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
        double a=(1-F[K+1][nlast])/(F_new[K+1][nlast]-F[K+1][nlast]);
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
    
    for (i=1;i<=n;i++)
    {
        double sum=0;
        for (k=1;k<=K;k++)
            sum += F[k][i];
        F[K+1][i]=sum;
    }
}


void compute_nabla(int N, int n, int K, int first[], int m[], int **freq, double **F, double lambda, double **nabla)
{

    int i,k;
    
    for (k=1;k<=K;k++)
    {
        for (i=1;i<=m[k];i++)
            nabla[k][i]=0;
    }
    
    for (k=1;k<=K;k++)
    {
        int j=0;
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

    for (int k=1;k<=K;k++)
    {
        int j=0;
        for (int i=first[k];i<=n;i++)
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
    int	i,j,n;
    
    for (int k=1;k<=K;k++)
    {
        int m1=m[k];
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


double bdf(double A, double B, int k, int njumps[], double **jumploc, double **p, double h, double u)
{

    double sum=0;
    
    for (int j=0;j<njumps[k];j++)
    {
        double t1=(u-jumploc[k][j])/h;
        double t2=(u+jumploc[k][j]-2*A)/h;
        double t3=(2*B-u-jumploc[k][j])/h;
        sum+= (KK(t1)+KK(t2)-KK(t3))*p[k][j];
    }
    return fmax(0,sum);
}

double dens_estimate(double A, double B, int k, int njumps[], double **jumploc, double **p, double h, double u)
{
    int i;
    double		rho,sum,x;
    weight_t	weights;
    
    sum=0;
    
    if (u>=A+h && u<=B-h)
    {
        for (i=0;i<njumps[k];i++)
        {
            x=(u-jumploc[k][i])/h;
            sum+= K(x)*p[k][i]/h;
        }
    }
    else
    {
        if (u<A+h)
        {
            rho=(u-A)/h;
            weights=weight(rho);
            for (i=0;i<njumps[k];i++)
            {
                x=(u-jumploc[k][i])/h;
                sum += (K(x)*weights.alpha + x*K(x)*weights.beta)*p[k][i]/h;
            }
        }
        else
        {
            if (u>B-h)
            {
                rho = (B-u)/h;
                weights=weight(rho);
                for (i=0;i<njumps[k];i++)
                {
                    x=(u-jumploc[k][i])/h;
                    sum += (K(x)*weights.alpha - x*K(x)*weights.beta)*p[k][i]/h;
                }
            }
        }
    }
    return fmax(0,sum);
}


void bootstrap(double A, double B, int NumIt, int N, int K, int ngrid, double grid[], double data[],
               double x2[], double data2[], int delta[], int delta2[], int **freq2, double **p2,
               double **jumploc2, int njumps2[], double **F2, int **ind2, double **dens2, double **Fsmooth2, double ***f3, double **hazard)
{
    int iter,i,k,iterations,n_Iterations=1000;
    double h1,h2,phi,c,tol=1.0e-10;
    
    c=B;
    h1=c*pow(N,-1.0/5);
    h2=fmin(c*pow(N,-1.0/6),9.9);
    
    srand((unsigned int)time(NULL));
    
    for (iter=0;iter<NumIt;iter++)
    {
        Rcpp::Rcout << std::setw(10) << iter+1 << "\r";
        data_bootstrap(N,data,x2,delta,delta2);
        
        for (k=0;k<=K;k++)
        {
            for (i=1;i<=N;i++)
                freq2[k][i]=0;
        }
        
        int j=0;
        
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
        
        int n2=j;
        
        
        ICM(N,n2,freq2,K,ind2,F2,n_Iterations,&phi,&iterations,tol);
        
        j=0;
        
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
        
        
        for (i=1;i<=99;i++)
        {
            for (k=1;k<=K;k++)
                f3[iter][k][i] = dens2[k][10*i]/(1-Fsmooth2[k][10*i])-hazard[k][10*i];
        }
        
    }
    
}

// the compare function for ordering results of bootstrap

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

// the compare function for ordering 2-dimensional array in data_bootstrap

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
    int	i;
    SampleTime *obs;
    
    obs = new SampleTime[n+1];
    
    for (i=1;i<=n;i++)
    {
        int j=1+rand()%n;
        
        x2[i]=x[j];
        delta2[i]=delta[j];
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







