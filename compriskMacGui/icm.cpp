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

void	ICM(int N, int n, int **freq, int K, int **ind, double **F,
                              int n_It, double *phi1, int *iteration)
{
	int		*first,nlast,i,j,k,iteration1,*m,*n1,m_total;
	double	**cumw,**cs,**v,**w,**nabla,**F_new;
	double	**y;
	double	sum,alpha,lambda,inprod,partialsum,eta=1.0e-10;
	
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
        y[i]= new double[N+1];
        F_new[i]= new double[N+1];
        cumw[i]= new double[N+1];
        cs[i]= new double[N+1];
        v[i]= new double[N+1];
        w[i]= new double[N+1];
        nabla[i]= new double[N+1];
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
    
	while (iteration1<=n_It && fenchelviol(K,m,ind,y,nabla,eta,&partialsum,&inprod))
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

/* The lines of input file are formatted as:
 time, cause
 where double "time" is the censoring time
 and "cause"=1,...K indicates the cause of failure
 and "cause"=0 means no failure has occurred yet.
 */


int	    CheckFileFormat(FILE *file, int nr_causes)
{
    int    cause;
    double time;
    
    rewind(file);
    return fscanf(file,"%lf %d", &time, &cause)==2;
}


void free_memory(int K, int NumIt, int **freq, double **F, double **Fsmooth, double **hazard, double **lowbound,
                 double **upbound, double data[], double data1[], int delta[], double grid[], double ***f3)
{
    int i,iter;
    for (i = 0; i < K+2; i++)
    {
        delete[] Fsmooth[i];
        delete hazard[i];
        delete F[i];
    }
    
    for (i = 0; i < K+1; i++)
        delete[] lowbound[i], delete upbound[i];
    
    for (i = 0; i < K+2; i++)
        delete[] freq[i];
    
    delete[] data;
    delete[] data1;
    delete[] delta;
    delete[] freq;
    delete[] Fsmooth;
    delete[] hazard;
    
    delete[] lowbound;
    delete[] upbound;
    delete[] grid;
    delete[] F;
    
    for (iter = 0;iter < NumIt;iter++)
    {
        for (i=0;i<K+1;i++)
            delete[] f3[iter][i];
    }
    
    for (iter = 0;iter < NumIt;iter++)
    {
        delete[] f3[iter];
    }
    
    delete[] f3;

}

double varF(int N, int n, int **freq, int k, int K1, double **F, int delta[], double A, double B, double t[], double h, double u)
{
    int			i,j,ind,count;
    double		t1,t2,t3,sum;
    
    
    sum=0;
    
    for (i=1;i<=n;i++)
    {
        t1=(u-t[i])/h;
        t2=(u+t[i]-2*A)/h;
        t3=(2*B-u-t[i])/h;
        
        if (delta[k]>0)
            ind =1;
        else
            ind=0;
        
        count=0;
        for (j=0;j<=K1;j++)
            count += freq[j][i];
        
        sum += SQR(K(t1)-K(t2)-K(t3))*(SQR(F[k][i]-1)*freq[k][i]+SQR(F[k][i])*(count-freq[k][i]));
    }
    
    sum = sum/(N*N);
    
    return sum;
}

int compute_K(int N, int delta[])
{
    int i,j=0;
    
    for (i=1;i<=N;i++)
    {
        if (delta[i]>j)
            j=delta[i];
    }
    
    return j;
}


int compute_n(int N, double data[])
{
    int i,j;
    
    j=0;
    
    for (i=1;i<=N;i++)
    {
        if (data[i]>data[i-1])
            j++;
    }
    
    return j;
}

/*void    reserve_memory1(int K,int N, int n, int ngrid, double **hazard, double *grid)
 {
 int k;
 
 grid= new double[ngrid+1];
 
 hazard= new double *[K+2];
 for (k=0;k<K+2;k++)
 hazard[k]= new double[ngrid+1];
 
 }
 
 void reserve_memory2(int K, int NumIt, int npoints, double **lowbound, double **upbound, double ***f3)
 {
 int k,iter;
 
 lowbound = new double *[K+1];
 upbound = new double *[K+1];
 
 for (k=0;k<K+1;k++)
 {
 lowbound[k] = new double[npoints];
 upbound[k] = new double[npoints];
 }
 
 f3  = new double**[NumIt];
 
 for (iter=0;iter<NumIt;iter++)
 f3[iter] = new double *[K+1];
 
 for (iter=0;iter<NumIt;iter++)
 for (k=0;k<K+1;k++)
 f3[iter][k] = new double[npoints];
 
 }*/


