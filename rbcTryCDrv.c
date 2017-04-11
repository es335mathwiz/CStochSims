

/*Mathematica Creation Date{2017, 4, 11, 15, 56, 59.365725}*/
/*rbc example model*/
#include "./lagLead.h"
#include <math.h>
#include "useSparseAMA.h"
#define aDummy(t)     (stateVector[(t-(-1))*4+0])
#define cc(t)     (stateVector[(t-(-1))*4+1])
#define kk(t)     (stateVector[(t-(-1))*4+2])
#define theta(t)     (stateVector[(t-(-1))*4+3])
#define linPt$aDummy(t)     (linearizationPoint[(t-(-1))*4+0])
#define linPt$cc(t)     (linearizationPoint[(t-(-1))*4+1])
#define linPt$kk(t)     (linearizationPoint[(t-(-1))*4+2])
#define linPt$theta(t)     (linearizationPoint[(t-(-1))*4+3])

#define modelShock(n) (shockVec[n])

  







void rbcExampleDerivative(double *stateVector,double *parameters,
double * shockVec,
double * aMat,int * jaMat,int *iaMat,double * homotopyAlpha,double * linearizationPoint
)
{int i;
double bMat[11];
int ibMat[4+1];
int jbMat[11];
double cMat[11];
int icMat[4+1];
int jcMat[11];
int aOne=1;int ierr;int maxNumberHElements;
int hrows=4;
int hcols=3*4;
int okay[50000];
if(*homotopyAlpha>=1.0) {
double okay10;
double okay14;
double okay4;
double okay5;
double okay6;
double okay8;




okay4=cc(tPaOne);

okay5=1/okay4;

okay6=kk(0);

okay10=pow(okay6,-0.64);

okay8=theta(0);

okay14=kk(tMaOne);

aMat[0]=-pow(cc(0),-2.);

aMat[1]=0.21888*okay5*okay8*pow(okay6,-1.6400000000000001);

aMat[2]=-0.34199999999999997*okay10*okay5;

aMat[3]=0.34199999999999997*okay10*okay8*pow(okay4,-2.);

aMat[4]=-0.36*okay8*pow(okay14,-0.64);

aMat[5]=1.;

aMat[6]=1.;

aMat[7]=-pow(okay14,0.36);

aMat[8]=-0.95*pow(theta(tMaOne),-0.050000000000000044);

aMat[9]=1.;

aMat[10]=1.;

iaMat[0]=6.;

iaMat[1]=7.;

iaMat[2]=8.;

iaMat[3]=10.;

iaMat[4]=3.;

iaMat[5]=6.;

iaMat[6]=7.;

iaMat[7]=8.;

iaMat[8]=4.;

iaMat[9]=8.;

iaMat[10]=5.;

jaMat[0]=1.;

jaMat[1]=5.;

jaMat[2]=9.;

jaMat[3]=11.;

jaMat[4]=12.;

} else {
double okay1;
double okay5;
double okay7;


okay1=linPt$cc(0);

okay5=linPt$kk(0);

okay7=linPt$theta(0);

aMat[0]=1/okay1-(0.34199999999999997*okay7*pow(okay5,-0.64))/linPt$cc(t\
PaOne);

aMat[1]=okay1+okay5-okay7*pow(linPt$kk(tMaOne),0.36);

aMat[2]=okay7-pow(linPt$theta(tMaOne),0.95);

aMat[3]=linPt$aDummy(0);

iaMat[0]=1.;

iaMat[1]=2.;

iaMat[2]=3.;

iaMat[3]=4.;

iaMat[4]=5.;

jaMat[0]=1.;

jaMat[1]=1.;

jaMat[2]=1.;

jaMat[3]=1.;

/*initialize cMat to zero sparse matrix*/
for(i=0;i<=4+1;i++){icMat[i]=1;}
for(i=0;i<11;i++){cMat[i]=0;};
if(*homotopyAlpha>0.0) {
double okay1;
double okay10;
double okay14;
double okay16;
double okay3;
double okay8;


okay1=cc(0);

okay14=kk(0);

okay3=linPt$cc(0);

okay8=linPt$kk(0);

okay10=linPt$theta(0);

okay16=theta(0);

aMat[0]=1/okay1-1/okay3-(0.34199999999999997*okay16*pow(okay14,-0.64))/\
cc(tPaOne)+(0.34199999999999997*okay10*pow(okay8,-0.64))/linPt$cc(tPaO\
ne);

aMat[1]=okay1+okay14-okay3-okay8-okay16*pow(kk(tMaOne),0.36)+okay10*pow\
(linPt$kk(tMaOne),0.36);

aMat[2]=-okay10+okay16+pow(linPt$theta(tMaOne),0.95)-pow(theta(tMaOne),\
0.95);

aMat[3]=aDummy(0)-linPt$aDummy(0);

iaMat[0]=1.;

iaMat[1]=2.;

iaMat[2]=3.;

iaMat[3]=4.;

iaMat[4]=5.;

jaMat[0]=1.;

jaMat[1]=1.;

jaMat[2]=1.;

jaMat[3]=1.;

for(i=0;i<11;i++){cMat[i]=cMat[i]*(*homotopyAlpha);};
}
maxNumberHElements=11;
aplb_(&hrows,&hcols,&aOne,bMat,jbMat,ibMat,cMat,jcMat,icMat,
aMat,jaMat,iaMat,&maxNumberHElements,okay,&ierr);
}
}
