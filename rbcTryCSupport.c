


/*Mathematica Creation Date{2017, 4, 5, 13, 12, 18.919406}*/
/*rbc example model*/
#include "lagLead.h"
#include <math.h>
//staticf double maxarg1,maxarg2;
#include <math.h>

double FMAX(double a,double b)
{
  return(a > b ? a : b);
}
double FMIN(double a,double b)
{
  return(a < b ? a : b);
}
double FABS(double a)
{
  return(a > 0 ? a: -a);
}

double doRightSmaller(double a,double b)
{
  return(a < b ? 0 : 1);
}
double doSign(double a)
{
  /*  return(fabs(a) >0.01?(a > 0 ? 1 : -1):2*a) ;*/
  return(a > 0 ? 1 : -1) ;

}
#define modelShock(n) (0)  


#define aDummy(t)     (stateVector[(t-(-1))*4+0])
#define cc(t)     (stateVector[(t-(-1))*4+1])
#define kk(t)     (stateVector[(t-(-1))*4+2])
#define theta(t)     (stateVector[(t-(-1))*4+3])
#define linPt$aDummy(t)     (linearizationPoint[(t-(-1))*4+0])
#define linPt$cc(t)     (linearizationPoint[(t-(-1))*4+1])
#define linPt$kk(t)     (linearizationPoint[(t-(-1))*4+2])
#define linPt$theta(t)     (linearizationPoint[(t-(-1))*4+3])

  




void rbcExamplePeriodicPointGuesser
(double * parameters,int period,
	double guessVector[3][4])
{
//int i,j;
//double svalue;
int timeOffset;
for(timeOffset=0;
	timeOffset<period+ 3 - 1;
			timeOffset++)
	{
guessVector[timeOffset][0]=0.;
guessVector[timeOffset][1]=0.35984508755628597;
guessVector[timeOffset][2]=0.18703194520402708;
guessVector[timeOffset][3]=1.;
}
}

void rbcExampleModelDimensions(int * numberOfEquations, int * lags,
int * leads, int * numberOfParameters,
int * numberOfDataValues, int * numberOfShocks,int * numberExogenous)
{
*numberOfEquations=4;
*lags=1;
*leads=1;
*numberOfParameters=0;
*numberOfDataValues=500;
*numberOfShocks=500;
*numberExogenous=0;
}
void rbcExampleUpsilon(double *parameters,
double * aMat,int * jaMat,int *iaMat
)
{
aMat[0]=1.;
iaMat[0]=1.;
iaMat[1]=2.;
jaMat[0]=1.;
}
void rbcExampleExogH(double *parameters,double *stateVector,
double * aMat,int * jaMat,int *iaMat
)
{
aMat[0]=1.;
iaMat[0]=1.;
iaMat[1]=2.;
jaMat[0]=1.;
}
void rbcExampleSelectZ(double * aMat,int * jaMat,int *iaMat
)
{
//aMat;
iaMat[0]=1.;
//jaMat;
}
