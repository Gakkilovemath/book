//
//  icm.h
//  trial_comprisk
//
//  Created by Piet on 07-02-15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

typedef struct
{
    double t;
    int delta;
}
SampleTime;


typedef struct{double alpha, beta;}
weight_t;

void    curstatgen(int n, double tt[], int delta[]);
void    main_bootstrap(int N, int K, double data[], int delta[],
                    double **hazard, double **lowbound, double **upbound, int ngrid, double grid[]);
int     CompareTime(const void *a, const void *b);
int     compare(const void *a, const void *b);
void	ICM(int ndata, int ndata1, int **freq, int K, int **ind, double **F,
                              int n_It, double *phi1, int *iteration);
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
void    cumsum(int m, double cs[], int delta[]);
void    convexmin(int n, double cumw[], double cs[], double w[], double y[]);
double  golden(int ndata, int ndata1, int K, int **freq, double **F, double **F_new, double a1, double b1,
              double (*f)(int,int,int,int**,double,double**,double**,double), double lambda);
double  dens_estimate(double A, double B,  int njumps, double jumploc[], double p[], double h, double u);
double  bdf(double A, double B, int k, int njumps[], double **jumploc, double **p, double h, double u);
void    data_bootstrap(int n, int *m, double x[], double x2[], double data2[], int **freq, int delta[], int delta2[]);
weight_t weight(double x);
double  K(double x);
double  KK(double x);
double  bdf(double B, int m, double t[], double p[], double u, double h);
void    data_binom(int n, int delta[], double F[]);
double  varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u);


