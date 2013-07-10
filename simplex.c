//simplex.c. 

/* Global variables -------------------------------------------------------- */
int NDIM;                   /* Dimension of approximation space        */
double* MX[25];             /* NDIM-dimensional vector of values       */
double MdX[25];             /* NDIM-dimensional vector of differences  */

/* Improved simplex approximation method ----------------------------------- */
void aps()
{
double MP[26][25];          /* Main matrix of simplex vertices         */
double Pb[25];              /* The point Pb.                           */
double Pr[25];              /* The point Pr.                           */
double Prr[25];             /* The point Prr or Prr'.                  */
double Y[26];               /* Results Y[i]=FUNK(MP[i])                */
double RTol;                /* Real toleration                         */
double Al,Bt,Gm,Ypr,Yprr,Yavr,Rmp;
double xa,xb,xc,xd,lMin;
int i,j,MPTS;
int ITR0,ITR1,NITR,MITR;
int iLo,iHi,iNHi;
FILE *fe;

MPTS=NDIM+1;Al=1.0;Bt=0.5;Gm=2.0;Rmp=MPTS;ITR0=0;
printf("Downhill simplex approximation method for %d dimensions.\n",NDIM);
if((fe=fopen("zparam00.dat","rt"))==0) {ITR1=0;MITR=100;NITR=40;} else {
  fscanf(fe,"%d %d\n",&ITR1,&NITR);MITR=100;
  for(i=0;i<MPTS;i++) fscanf(fe,"%le %le",MX[i],MdX+i);
  fclose(fe);}
while(ITR1<NITR) {ITR0=0;
 for(i=0;i<NDIM;i++) {for(j=0;j<=NDIM;j++) MP[j][i]=*(MX[i]);MP[i][i]+=MdX[i]*(0.9+0.2*rand()/RAND_MAX);}
 /* Rotarion of the matrix MP - not implemented */
 for(j=0;j<=NDIM;j++) {for(i=0;i<NDIM;i++) *(MX[i])=MP[j][i];Y[j]=func();}
 while(ITR0<MITR) {
  Yavr=0;iLo=0;if(Y[0]>Y[1]) {iHi=0;iNHi=1;} else {iHi=1;iNHi=0;}
  for(i=0;i<MPTS;i++) {Yavr+=Y[i];if(Y[i]<Y[iLo]) iLo=i;
   if(Y[i]>Y[iHi]) {iNHi=iHi;iHi=i;} else if(Y[i]>Y[iNHi]) if(i!=iHi) iNHi=i;}
  ITR0++;RTol=2.0*fabs(Y[iHi]-Y[iLo])/(fabs(Y[iHi])+fabs(Y[iLo]));Yavr/=Rmp;
  printf(">%3d%4d RTol=%e Ymin=%f Ymax=%f Yavr=%f\n",ITR1,ITR0,RTol,Y[iLo],Y[iHi],Yavr);
  for(j=0;j<NDIM;j++) Pb[j]=0.0;
  for(i=0;i<MPTS;i++) if(i!=iHi) for(j=0;j<NDIM;j++) Pb[j]+=MP[i][j];
  for(j=0;j<NDIM;j++) {Pb[j]/=NDIM;Pr[j]=(1.0+Al)*Pb[j]-Al*MP[iHi][j];}
  for(j=0;j<NDIM;j++) *(MX[j])=Pr[j];Ypr=func();
  if(Ypr<=Y[iLo]) {for(j=0;j<NDIM;j++) Prr[j]=Gm*Pr[j]+(1.0-Gm)*Pb[j];
    for(j=0;j<NDIM;j++) *(MX[j])=Prr[j];Yprr=func();
    if(Ypr>Yprr) {for(j=0;j<NDIM;j++) MP[iHi][j]=Prr[j];Y[iHi]=Yprr;}
    else {for(j=0;j<NDIM;j++) MP[iHi][j]=Pr[j];Y[iHi]=Ypr;}}
  else {
    if(Ypr>=Y[iNHi]) {
      if(Ypr<Y[iHi]) {for(j=0;j<NDIM;j++) MP[iHi][j]=Pr[j];Y[iHi]=Ypr;}
      for(j=0;j<NDIM;j++) Prr[j]=Bt*MP[iHi][j]+(1.0-Bt)*Pb[j];
      for(j=0;j<NDIM;j++) *(MX[j])=Prr[j];Yprr=func();
      if(Yprr<Y[iHi]) {for(j=0;j<NDIM;j++) MP[iHi][j]=Prr[j];Y[iHi]=Yprr;}
      else {
/*      for(i=0;i<MPTS;i++) if(i!=iLo) {
          for(j=0;j<NDIM;j++) {Pr[j]=0.5*(MP[i][j]+MP[iLo][j]);MP[i][j]=Pr[j];}
          for(j=0;j<NDIM;j++) *(MX[j])=Pr[j];Y[i]=func();}}}                     */
        for(j=0;j<NDIM;j++) Pr[j]=0.5*(MP[iHi][j]+MP[iLo][j]);
        for(j=0;j<NDIM;j++) *(MX[j])=Pr[j];Ypr=func();
        if(Ypr<Y[iHi]) {for(j=0;j<NDIM;j++) MP[iHi][j]=Pr[j];Y[iHi]=Ypr;}
        else {
          for(j=0;j<NDIM;j++) Prr[j]=-MP[iHi][j]+2.0*MP[iLo][j];
          for(j=0;j<NDIM;j++) *(MX[j])=Prr[j];Yprr=func();
          if(Yprr<Y[iHi]) {for(j=0;j<NDIM;j++) MP[iHi][j]=Prr[j];Y[iHi]=Yprr;}
          else {
            xa=3*Y[iHi]-8*Ypr+6*Y[iLo]-Yprr;
            xb=Y[iHi]-2*Y[iLo]+Yprr;
            xc=-0.5*Y[iHi]+8*Ypr/3-2*Y[iLo]+Yprr/6;
            xd=xb*xb-4*xa*xc;
            if(xd>0) {
              lMin=0.5*(-xb-sqrt(xd))/xa;
              for(j=0;j<NDIM;j++) Pr[j]=lMin*MP[iHi][j]+(1-lMin)*MP[iLo][j];
              for(j=0;j<NDIM;j++) *(MX[j])=Pr[j];Ypr=func();}
            if(Ypr<Y[iHi]) {for(j=0;j<NDIM;j++) MP[iHi][j]=Pr[j];Y[iHi]=Ypr;}
            else {for(j=0;j<NDIM;j++) MP[iHi][j]=MP[iLo][j];Y[iHi]=Y[iLo];}
          }}}}
    else {for(j=0;j<NDIM;j++) MP[iHi][j]=Pr[j];Y[iHi]=Ypr;}}}
 iLo=0;for(i=0;i<MPTS;i++) if(Y[i]<Y[iLo]) iLo=i;
 for(i=0;i<NDIM;i++) *(MX[i])=MP[iLo][i];ITR1++;
 fe=fopen("zraport0.dat","at");fprintf(fe,"%5d %f  ",ITR1,Y[iLo]);
 for(i=0;i<NDIM;i++) {fprintf(fe,"%16.8e ",*MX[i]);}
 fprintf(fe,"\n");fclose(fe);
 printf("=======> Ymin=%f  (%d)\n",Y[iLo],iLo);
 /* for(i=0;i<NDIM;i++) {MdX[i]*=0.95;} */
 fe=fopen("zparam00.dat","wt");fprintf(fe,"%d %d\n",ITR1,NITR);
 for(i=0;i<NDIM;i++) {fprintf(fe,"%16.8e %16.8e\n",*MX[i],MdX[i]);}fclose(fe);}
printf("Approximation with simplex method finished.\n");
}
