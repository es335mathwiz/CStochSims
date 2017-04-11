/*Mathematica Creation Date {2017, 4, 11, 15, 20, 55.270366}*/
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

void rbcExample(double *stateVector,double *parameters,
double * shockVec,
double * aMat,int * jaMat,int *iaMat,double * homotopyAlpha,double * linearizationPoint
)
{
int i;
double bMat[4];
//int ibMat[4+1];
//int jbMat[4];
if(*homotopyAlpha>=1.0) {
double okay1;
double okay5;
double okay7;


okay1=cc(0);

okay5=kk(0);

okay7=theta(0);

aMat[0]=1/okay1-(0.34199999999999997*okay7*pow(okay5,-0.64))/cc(tPaOne)\
;

aMat[1]=okay1+okay5-okay7*pow(kk(tMaOne),0.36);

aMat[2]=okay7-pow(theta(tMaOne),0.95);

aMat[3]=aDummy(0);

for(i=0;i<4-0;i++){aMat[i]=aMat[i]+shockVec[i];};
iaMat[0]=1.;

iaMat[1]=2.;

iaMat[2]=3.;

iaMat[3]=4.;

iaMat[4]=5.;

jaMat[0]=1.;

jaMat[1]=1.;

jaMat[2]=1.;

jaMat[3]=1.;

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

for(i=0;i<4-0;i++){aMat[i]=aMat[i]+shockVec[i];};
iaMat[0]=1.;

iaMat[1]=2.;

iaMat[2]=3.;

iaMat[3]=4.;

iaMat[4]=5.;

jaMat[0]=1.;

jaMat[1]=1.;

jaMat[2]=1.;

jaMat[3]=1.;

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

for(i=0;i<4;i++){aMat[i]=aMat[i]+(*homotopyAlpha*bMat[i]);};
}
}
}
