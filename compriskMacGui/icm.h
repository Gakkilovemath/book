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

int     GetFileSize(FILE *file);
int     GetFileSize2(int K, FILE *file);
int     GetNumberElements_first_line(FILE *file);
void main_comp(double *A1, double *B1, int N, int n, int K, int ngrid, double data[], double data1[], int delta[],
               int **freq, double *grid, double **F, double **Fsmooth, double **hazard, const char *appdir);
void    main_comp_SMLE_bootstrap(int iter, double A, double B, int N, int n, int K, int ngrid, double data[], int delta[], double grid[], double **Fsmooth, double ***f3, const char *appdir);
void    main_comp_hazard_bootstrap(int iter, double A, double B, int N, int n, int K, int ngrid,
                            double data[], int delta[], double grid[], double **hazard, double ***f3, const char *appdir);
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
void    cumsum(int K, int m[], double **v, double **cs);
void    convexmin(int K, int m[], double **cumw, double **cs, double **w, double **y);
double  golden(int ndata, int ndata1, int K, int **freq, double **F, double **F_new, double a1, double b1,
              double (*f)(int,int,int,int**,double,double**,double**,double), double lambda);
void    data_bootstrap(int n, double x[], double x2[], int delta[], int delta2[]);
void    confidence_intervals_SMLE(double A, double B, int N, int n, int **freq, int npoints, int K, int NumIt,
                          double grid[], double **F, double data[], int delta[], double **Fsmooth, double **lowbound, double **upbound, double ***f3);
void    confidence_intervals_hazard(int npoints, int K, int NumIt, double grid[], double **hazard,
                                 double **lowbound, double **upbound, double ***f3);



double  dens_estimate(double A, double B, int k, int njumps[], double **jumploc, double **p, double h, double u);
double  bdf(double A, double B, int k, int njumps[], double **jumploc, double **p, double h, double u);
void    bootstrap_SMLE(int iter, double A, double B, int N, int n, int K, int ngrid, double grid[], double data[],double x2[], double data2[], int delta[], int delta2[], int **freq2, double **p2,double **jumploc2, int njumps2[], double **F2, int **ind2, double **dens2, double **Fsmooth2, double ***f3, double **hazard);
void    bootstrap_hazard(int iter, double A, double B, int N, int n, int K, int ngrid, double grid[],
                         double data[],double x2[], double data2[], int delta[], int delta2[], int **freq2, double **p2,double **jumploc2, int njumps2[], double **F2, int **ind2, double **dens2, double **Fsmooth2, double ***f3, double **hazard);
weight_t weight(double x);
double  K(double x);
double  KK(double x);
int	    CheckFileFormat(FILE *file, int nr_causes);
double  varF(int N, int n, int **freq, int k, int K, double **F, int delta[], double A, double B, double t[], double h, double u);
int compute_n(int N, double data[]);
int compute_K(int N, int delta[]);

/*void   reserve_memory1(int K,int N, int n, int ngrid, double **hazard, double *grid);
void   reserve_memory2(int K, int NumIt, int npoints, double **lowbound, double **upbound, double ***f3);*/

void    free_memory(int K, int NumIt, int **freq, double **F, double **Fsmooth, double **hazard, double **lowbound, double **upbound, double data[], double data1[], int delta[], double grid[], double ***f3);


