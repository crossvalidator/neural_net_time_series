#include "stdio.h"
#include  "math.h"
#include "stdlib.h"
#define ALPHA 1
#define BETA 0.5
#define GAMA 2

double tune_variable[N]={0.0,0.0,0.0,0.0};

/******************************************************************/
void downh(double *magnitude, double *phase)
{
}
/******************************************************************************/ 
double *Downhill(double *magnitude, double *phase)
{
  	double fminv=0.0,final=0.0;
	double temperature=T;
	int i=0;
	double x[N],f[M];
	double a[N]={0.0,0.0,0.0,0.0001};
	double b[N]={1.0e+1,1.0e+2,1.0,0.5};
	double coea[N]={1, 10, 0.1, 0.04999};//1, 10, 0.1, 0.09
	double *optx;
	optx = calloc(N, sizeof(double));

	downh(magnitude,phase);
	for (i=0;i<N;i++)
		printf("optx[%d]=%f\n",i,optx[i]);
	printf("fminv=%f\n",fminv);
	return optx;
   
}
