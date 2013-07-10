//#include "ancfis.h"
#include "stdio.h"
#include  "math.h"
#include "stdlib.h"
#define N 4
#define M 800 //2000
#define K 2
#define T 100 //500//100
#define NUMSTEPS 2 //10
#define ALPHA 0.99
#define WEIGHT 0.95
double rr=1.0;//1.0
int itern=0;
double tune_variable[N]={0.0,0.0,0.0,0.0};
double rn(double *ran_point)
{
int m;
double s,u,v,p;
s=65536.0;
u=2053.0;
v=13849.0;
m=(int)(*ran_point/s);
*ran_point=*ran_point-m*s;
*ran_point=u*(*ran_point)+v;
m=(int)(*ran_point/s);
*ran_point=*ran_point-m*s;
p=*ran_point/s;
return p;
}
 double f_value(double x[], double *magnitude, double *phase )
{
	int j;
	double final=0.0;
	for(j = 0 ; j < In_vect_n*In_vect_n ; j++) 
	final += pow((magnitude[j]-fabs(x[3]*sin(x[0]*phase[j]+x[1])+x[2])), 2.0);
	return final;
}
void fun_restrict(double x[], double c[], double d[], double w[])
{
	   c[0]=0.0; c[1]=0.0;
       d[0]=x[2];  d[1]=1.0;
       w[0]=x[3];  w[1]=fabs(x[3])+fabs(x[2]);
       /*c[0]=0.0;
       d[0]=1.0;
       w[0]=fabs(x[3])+fabs(x[2]);*/	   
}

double MonteCarlo(double a[], double b[], double x[],double f[],double optx[], double coea[], double *fmin, double *t,double *magnitude, double *phase)
{
        extern double f_value();
		extern void fun_restrict();
		double rn();
        int i=0,j=0,r=0,g=0;
		double c[K],d[K],w[K],cx[N],cy[N][M],tempoptx[N];
		double temp=0,temp1=0, tempminv=0.0,randnumber=0.0, annealingTemperature=0.0,funminv=0.0;
		tempminv=*fmin;
		funminv=*fmin;
		for(i=0;i<N;i++)
			tempoptx[i]=optx[i];

		for (i=0;i<N;i++)
		{
		  cx[i]=rn(&rr);
		  for (j=0; j<M;j++)
		  {
			  cy[i][j]=1-2*cx[i]*cx[i];
			  cx[i]=cy[i][j];
		  }
		}
		for(j=0;j<M;j++)
		{
			for(i=0;i<N;i++)
		       x[i]=tempoptx[i]+coea[i]*cy[i][j];
			g=0;r=0;
			while((r<N) &&(g==0))
			{
				if((a[r]<=x[r])&&(b[r]>=x[r]))
					r=r+1;
				else
					g=1;
			}
			if(g==0)
			{
				fun_restrict(x,c,d,w);
				r=0;
				while((r<K) && (g==0))
				{
					if((c[r]<=w[r])&&(d[r]>=w[r]))
						r=r+1;
					else
						g=1;
				}
			}
			if(g==0)
			{
				for(i=0;i<N;i++)
				if (tune_variable[i]<fabs(x[i]-tempoptx[i])){
			       tune_variable[i]=fabs(x[i]-tempoptx[i]);
				}
				f[j]=f_value(x, magnitude, phase);
                if(tempminv>f[j])
					
				{ for(i=0;i<N;i++)
				  tempoptx[i]=x[i];
				  tempminv=f[j];
				  if(funminv>f[j])
				  {
					  for(i=0;i<N;i++)

					  optx[i]=x[i];
					 
					funminv=f[j];
				  }
				  
				}
				else
				{ randnumber=rn(&rr);
				 annealingTemperature=*t;
				 if(randnumber<exp(-(f[j]-tempminv)/annealingTemperature))
				 {
					 for(i=0;i<N;i++)
				       tempoptx[i]=x[i];
          			 tempminv=f[j];
				 }
				}
			}
			else
				f[j]=1.0E+010;
		}	
	   
*fmin=funminv;
return(*fmin);
}	
/******************************************************************/
void chaos(double a[], double b[], double x[],double f[],double optx[], double coea[], double *minv, double *t,double *magnitude, double *phase)
	{
		extern double f_value();
		extern void fun_restrict();
		extern double MonteCarlo();
		double rn();
		int i=0,j=0,l=0,p=0,r=0,g=0;
		double c[K],d[K],w[K],cx[N],cy[N][M],xx[N][M],temp=0,temp1=0, function_minimum;
		for (i=0;i<N;i++)
		{
		  cx[i]=0.1;
		  for (j=0; j<M;j++)
		  {
			  cy[i][j]=4*cx[i]*(1-cx[i]);
			  cx[i]=cy[i][j];
		  }
		}
		for(j=0;j<M;j++)
		{
			for(i=0;i<N;i++)
		       {
				xx[i][j]=a[i]+(b[i]-a[i])*cy[i][j];
				x[i]=xx[i][j];
			}
			g=0;
			if(g==0)
			{
				fun_restrict(x,c,d,w);
				r=0;
				while((r<K) && (g==0))
				{
					if((c[r]<=w[r])&&(d[r]>=w[r]))
						r=r+1;
					else
						g=1;
				}
			}
			if(g==0)
			{
				f[j]=f_value(x, magnitude, phase);
			}
			else
				f[j]=1.0E+010;

		}
		for(j=0;j<M;j++)
		{
			p=j;
			for( l=j+1;l<M;l++)
			  if(f[l]<f[p]) p=l;
			  if(p!=j)
			  {
			    temp=f[p];
			    f[p]=f[j];
			    f[j]=temp;
			    for (i=0;i<N;i++)
			    {
				temp1=xx[i][p];
			        xx[i][p]=xx[i][j];
				xx[i][j]=temp1;
			    }
			  }
		}
		 for (i=0;i<N;i++)
		    optx[i]=xx[i][0];
		 
		 *minv=f[0];
	     while (*t>0.01)//0.01 sd: this is Tmin
		 {
             for (j=0;j<NUMSTEPS;j++) // Appears to be Lmax
			  function_minimum=MonteCarlo(a,b,x,f,optx,coea,minv,t,magnitude,phase);
			 *t*=0.98; // sd: this is parameter Beta.
			 /*for (i=0;i<N;i++)
			  coea[i]=coea[i]*0.99;*/
			 for (i=0;i<N;i++){
			   coea[i]=coea[i]*(1-ALPHA)+ALPHA*WEIGHT*tune_variable[i]; 			  
			   tune_variable[i]=0.0;}
			 //printf("t=%f\n",*t);
		 }
}
/******************************************************************************/ 
double *chaotic_simulatedAnnealling(double *magnitude, double *phase)
{
  	extern void chaos();
	double fminv=0.0,final=0.0;
	double temperature=T;
	int i=0;
	double x[N],f[M];
	double a[N]={0.0,0.0,0.0,0.0001};
	double b[N]={1.0e+1,1.0e+2,1.0,0.5};
	double coea[N]={1, 10, 0.1, 0.04999};//1, 10, 0.1, 0.09
	double *optx;
	optx = calloc(N, sizeof(double));

	chaos(a,b,x,f,optx,coea,&fminv, &temperature,magnitude,phase);
	for (i=0;i<N;i++)
		printf("optx[%d]=%f\n",i,optx[i]);
	printf("fminv=%f\n",fminv);
	return optx;
   
}
