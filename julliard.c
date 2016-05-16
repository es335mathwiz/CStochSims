/*Mathematica Creation Date{2016, 2, 25, 15, 0, 22.090232}*/
/*julliard example from paper describing stack*/
#include "../lagLead.h"
#include <math.h>
#define ey(t)     (stateVector[(t-(-1))*5+0])
#define pdot(t)     (stateVector[(t-(-1))*5+1])
#define rr(t)     (stateVector[(t-(-1))*5+2])
#define rs(t)     (stateVector[(t-(-1))*5+3])
#define y(t)     (stateVector[(t-(-1))*5+4])
#define linPt$ey(t)     (linearizationPoint[(t-(-1))*5+0])
#define linPt$pdot(t)     (linearizationPoint[(t-(-1))*5+1])
#define linPt$rr(t)     (linearizationPoint[(t-(-1))*5+2])
#define linPt$rs(t)     (linearizationPoint[(t-(-1))*5+3])
#define linPt$y(t)     (linearizationPoint[(t-(-1))*5+4])

#define modelShock(n) (shockVec[n])

  







void julMod(double *stateVector,double *parameters,
double * shockVec,
double * aMat,int * jaMat,int *iaMat,double * homotopyAlpha,double * linearizationPoint
)
{
int i;
double bMat[5];
int ibMat[6];
int jbMat[5];
if(*homotopyAlpha>=1.0) {


aMat[0]=1.;
aMat[1]=1.;
aMat[2]=1.;
aMat[3]=1.;
aMat[4]=1.;

/*for(i=0;i<5;i++){aMat[i]=aMat[i]+shockVec[i];};*/
iaMat[0]=1.;
iaMat[1]=1.;
iaMat[2]=1.;
iaMat[3]=1.;
iaMat[4]=1.;
jaMat[0]=1.;
jaMat[1]=2.;
jaMat[2]=3.;
jaMat[3]=4.;
jaMat[4]=5.;
jaMat[5]=6.;
} else {


aMat[0]=1.;
aMat[1]=1.;
aMat[2]=1.;
aMat[3]=1.;
aMat[4]=1.;

/*for(i=0;i<5;i++){aMat[i]=aMat[i]+shockVec[i];};*/
iaMat[0]=1.;
iaMat[1]=1.;
iaMat[2]=1.;
iaMat[3]=1.;
iaMat[4]=1.;
jaMat[0]=1.;
jaMat[1]=2.;
jaMat[2]=3.;
jaMat[3]=4.;
jaMat[4]=5.;
jaMat[5]=6.;
if(*homotopyAlpha>0.0) {


aMat[0]=1.;
aMat[1]=1.;
aMat[2]=1.;
aMat[3]=1.;
aMat[4]=1.;

iaMat[0]=1.;
iaMat[1]=1.;
iaMat[2]=1.;
iaMat[3]=1.;
iaMat[4]=1.;
jaMat[0]=1.;
jaMat[1]=2.;
jaMat[2]=3.;
jaMat[3]=4.;
jaMat[4]=5.;
jaMat[5]=6.;
for(i=0;i<5;i++){aMat[i]=aMat[i]+(*homotopyAlpha*bMat[i]);};
}
}
}

