

/*Mathematica Creation Date{2017, 4, 5, 13, 12, 18.919406}*/
/*rbc example model*/
#include "useSparseAMA.h"
#include "/a/msulx1/lcl/msu/home/m1gsa00/git/CStochSims/lagLead.h"
#include <math.h>
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

jaMat[0]=1.;

jaMat[1]=5.;

jaMat[2]=9.;

jaMat[3]=11.;

jaMat[4]=12.;

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

} else {
double okay1;
double okay10;
double okay11;
double okay15;
double okay2;
double okay21;
double okay24;
double okay25;
double okay26;
double okay30;
double okay35;
double okay39;
double okay8;
double okay9;


okay2=linPt$cc(0);

okay11=linPt$kk(0);

okay8=linPt$cc(tPaOne);

okay15=linPt$theta(0);

okay9=1/okay8;

okay21=pow(okay11,-0.64);

okay1=cc(0);

okay10=kk(0);

okay30=linPt$kk(tMaOne);

okay35=pow(okay30,0.36);

okay24=-okay15;

okay25=theta(0);

okay26=okay24+okay25;

okay39=linPt$theta(tMaOne);

aMat[0]=1/okay2-0.34199999999999997*okay15*okay21*okay9-0.3419999999999\
9997*okay21*okay26*okay9+0.21888*(okay10-okay11)*okay15*okay9*pow(okay\
11,-1.6400000000000001)-(okay1-okay2)*pow(okay2,-2.)+0.341999999999999\
97*okay15*okay21*(-okay8+cc(tPaOne))*pow(okay8,-2.);

aMat[1]=okay1+okay10-okay15*okay35-okay26*okay35-0.36*okay15*(-okay30+k\
k(tMaOne))*pow(okay30,-0.64);

aMat[2]=okay25-pow(okay39,0.95)-0.95*pow(okay39,-0.050000000000000044)*\
(-okay39+theta(tMaOne));

aMat[3]=aDummy(0);

jaMat[0]=1.;

jaMat[1]=2.;

jaMat[2]=3.;

jaMat[3]=4.;

jaMat[4]=5.;

iaMat[0]=1.;

iaMat[1]=1.;

iaMat[2]=1.;

iaMat[3]=1.;

/*initialize cMat to zero sparse matrix*/
for(i=0;i<=4+1;i++){icMat[i]=1;}
for(i=0;i<11;i++){cMat[i]=0;};
if(*homotopyAlpha>0.0) {
double okay1;
double okay10;
double okay11;
double okay12;
double okay13;
double okay17;
double okay19;
double okay23;
double okay28;
double okay3;
double okay30;
double okay31;
double okay34;
double okay35;
double okay40;
double okay46;
double okay48;


okay1=cc(0);

okay3=linPt$cc(0);

okay13=linPt$kk(0);

okay10=linPt$cc(tPaOne);

okay17=linPt$theta(0);

okay11=1/okay10;

okay23=pow(okay13,-0.64);

okay19=cc(tPaOne);

okay12=kk(0);

okay28=theta(0);

okay35=linPt$kk(tMaOne);

okay34=kk(tMaOne);

okay40=pow(okay35,0.36);

okay30=-okay17;

okay31=okay28+okay30;

okay46=linPt$theta(tMaOne);

okay48=theta(tMaOne);

aMat[0]=1/okay1+0.34199999999999997*okay11*okay17*okay23-1/okay3+0.3419\
9999999999997*okay11*okay23*okay31-0.34199999999999997*okay17*(-okay10\
+okay19)*okay23*pow(okay10,-2.)-(0.34199999999999997*okay28*pow(okay12\
,-0.64))/okay19-0.21888*okay11*(okay12-okay13)*okay17*pow(okay13,-1.64\
00000000000001)+(okay1-okay3)*pow(okay3,-2.);

aMat[1]=okay17*okay40+okay31*okay40-okay28*pow(okay34,0.36)+0.36*okay17\
*(okay34-okay35)*pow(okay35,-0.64);

aMat[2]=0.95*(-okay46+okay48)*pow(okay46,-0.050000000000000044)+pow(oka\
y46,0.95)-pow(okay48,0.95);

aMat[3]=0.;

jaMat[0]=1.;

jaMat[1]=2.;

jaMat[2]=3.;

jaMat[3]=4.;

jaMat[4]=5.;

iaMat[0]=1.;

iaMat[1]=1.;

iaMat[2]=1.;

iaMat[3]=1.;

for(i=0;i<11;i++){cMat[i]=cMat[i]*(*homotopyAlpha);};
}
maxNumberHElements=11;
aplb_(&hrows,&hcols,&aOne,bMat,jbMat,ibMat,cMat,jcMat,icMat,
aMat,jaMat,iaMat,&maxNumberHElements,okay,&ierr);
}
}
