



/*Mathematica Creation Date{2017, 4, 5, 13, 12, 18.919406}*/
/*rbc example model*/
#include <stdlib.h>
#include "runItExternalDefs.h"



#define PATHLENGTH 1000

int numberOfEquations=4;
char * namesArray[] =  {"aDummy", "cc", "kk", "theta"};
char * paramNamesArray[] = {};
int numberOfParameters=0;
int * parameters[]={};
int numDATA=500;
int numSHOCKS=500;
double * theData;

main(int argc, char * argv[])
{
#include "runItInvariantLocalDefs.h"
#include "runrbcTryCLocalDefs.h"
printf(" runIt.mc, 2016 m1gsa00 \n");

rbcExampleDataVals=(double *)calloc(numberOfEquations*numDATA,sizeof(double));
for(i=0;i<numDATA;i++){rbcExampleData(i,rbcExampleDataVals+(i*numberOfEquations));}

rbcExampleShockVals=(double *)calloc(numberOfEquations*numSHOCKS,sizeof(double));
for(i=0;i<numSHOCKS;i++){rbcExampleShocks(i,rbcExampleShockVals+(i*numberOfEquations));}


processCommandLine(argc,argv,namesArray,*numberOfEquations,
paramNamesArray,numberOfParameters,parameters,
	rbcExampleDataVals,numDATA,numSHOCKS,
	pathLength,replications,t0,stochasticPathLength,
intControlParameters,doubleControlParameters,flnm);

/*
rbcExamplePeriodicPointGuesser(parameters,1,rbcExampleFP);

FPnewt(numberOfEquations,lags,leads,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
failedQ);

*/
}

#include "runItOther.h"

/*
printf("generating perm vec
");
 generateDraws(1,(stochasticPathLength),(*replications),numSHOCKS,julliardPermVec);
printf("done generating perm vec
");
*/


/*


altComputeAsymptoticQMatrix(
numberOfEquations,lags,leads,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,pathLength,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
AMqMatrix,AMqMatrixj,AMqMatrixi,
failedQ
);
*/

/*
if(failedQ[0])
{  printf("problems computing  Q matrix
");return(1);}
else {printf("computed Q matrix
");}
*//*
printf("computed Q matrix
");
for(i=0;i< *pathLength;i++){
rbcExamplePeriodicPointGuesser(parameters,1,
julliardPathQ+(i *julNEQS));}
*totalTime=dtime(userSystemTime);
printf("after computing Q matrix
totalTime=%f,userSystemTime=%f,systemTime=%f
",
*totalTime,*userSystemTime,*(userSystemTime+1));
printf("using q matrix
");*/
/*

stochSim(numberOfEquations,lags,leads,pathLength,
rbcExample,rbcExampleDerivative,parameters,
replications,t0,tf,julliardPermVec,
julliardShocks,numberOfShocks,
julliardData,numberOfData,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,AMqMatrix,AMqMatrixj,AMqMatrixi,
julliardFP,
julliardPathQ,
failedQ);

*/

