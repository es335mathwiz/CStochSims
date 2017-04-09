



/*Mathematica Creation Date{2017, 4, 5, 13, 12, 18.919406}*/
/*rbc example model*/
#include <stdlib.h>
#include <stdio.h>
#include "stackC.h"
#include "stochProto.h"

// not needed #include "runItExternalDefs.h"




//char * namesArray[] =  {"aDummy", "cc", "kk", "theta"};
//char * paramNamesArray[] = {};
unsigned int Parameters=0;
//int * parameters[]={};
unsigned int NEQS=4;
unsigned int NLAGS=1;
unsigned int NLEADS=1;
unsigned int numDATA=500;
unsigned int numSHOCKS=500;
double * theData;

#include "runItInvariantLocalDefs.h"
#include "runrbcTryCLocalDefs.h"

int main(int argc, char * argv[])
{
printf(" runIt.mc, 2016 m1gsa00 \n");

rbcExampleDataVals=(double *)calloc((*numberOfEquations)*numDATA,sizeof(double));
for(i=0;i<numDATA;i++){rbcExampleData(i,rbcExampleDataVals+(i*(*numberOfEquations)));}

rbcExampleShockVals=(double *)calloc((*numberOfEquations)*numSHOCKS,sizeof(double));
for(i=0;i<numSHOCKS;i++){rbcExampleShocks(i,rbcExampleShockVals+(i*(*numberOfEquations)));}


processCommandLine(argc,argv,namesArray,*numberOfEquations,
paramNamesArray,numberOfParameters,parameters,
	rbcExampleDataVals,numDATA,numSHOCKS,
	pathLength,replications,t0,stochasticPathLength,
intControlParameters,doubleControlParameters,flnm);
unsigned int exogRows[0];
unsigned int exogCols[0];
unsigned int exogenizeQ[1]={0};


unsigned int * rbcExampleFailedQ;
rbcExampleFailedQ=(unsigned int *)calloc(*replications,sizeof(unsigned int));
unsigned int maxNumberElements=100;
double ** fmats, ** smats;
unsigned int ** fmatsj, ** smatsj,** fmatsi, **smatsi;
allocFPNewt(*numberOfEquations,NLAGS,NLEADS,*pathLength,maxNumberElements,
&rbcExampleFP,&rbcExampleIntercept,
&fmats,&fmatsj,&fmatsi,
&smats,&smatsj,&smatsi);

rbcExamplePeriodicPointGuesser(parameters,1,rbcExampleFP);
FPnewt(numberOfEquations,&NLAGS,&NLEADS,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,rbcExampleIntercept,exogRows,exogCols,exogenizeQ,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
&maxNumberElements,
rbcExampleFailedQ,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo);

FILE * outFile;
outFile=fopen(flnm,"w");

printf("saving values for variable in file named %s\n",flnm);
fprintf(outFile,"RunParams={%d,%d,%d,%d,%d,%d,%d};\n",
    NEQS,NLAGS,NLEADS,
     *pathLength,*t0,*stochasticPathLength,*replications);


freeFPNewt(NLAGS,*pathLength,
&rbcExampleFP,&rbcExampleIntercept,
&fmats,&fmatsj,&fmatsi,
&smats,&smatsj,&smatsi);

/*
fPrintMathInt(outFile,*replications,rbcExampleFailedQ,"rbcExampleFailedQ");
fPrintMathInt(outFile,*replications * (*stochasticPathLength),
      rbcExamplePermVec,"rbcExamplePermVec");
fPrintMathDbl(outFile,(*replications * NEQS*(*stochasticPathLength+NLAGS)),
      rbcExamplePathQ,"Results");
fPrintMathDbl(outFile,(NEQS*(numDATA)),rbcExampleDataVals,"dataArray");
fPrintMathDbl(outFile,(NEQS*(numSHOCKS)),rbcExampleShockVals,"shocksArray");
*/

     fclose(outFile);

return(0);

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

