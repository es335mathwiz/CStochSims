%	$Id: stochRun.w,v 1.6 2000/12/06 14:53:34 m1gsa00 Exp m1gsa00 $	
\documentclass[html]{article}
\usepackage{moreverb}
\begin{document}

\title{``C'' Program Calling stochSim()
$Author: m1gsa00 $\\
%\begin{minipage}{\textwidth}\verbatim{$Date: 2000/12/06 14:53:34 $}\end{minipage}
%\footnote{\begin{minipage}{\textwidth}\verbatim{$Header: /mq/home4/m1gsa00/aim/stochSims/RCS/stochRun.w,v 1.6 2000/12/06 14:53:34 m1gsa00 Exp m1gsa00 $}\end{minipagez}}
}
\author{Gary Anderson}
\maketitle
 
\centerline{{\bf Abstract}}
\begin{quote}
  This paper provides an example calling the prototype stochSim() routine.
\end{quote}



\newpage
\tableofcontents
\newpage


\section{How to Use This Document}
\label{sec:intro}

The numbers just inside the angle bracket refers to the page where the
scrap is located.

\section{Compiling and Running}
\label{sec:comprun}


\begin{description}

\item[{\bf creating model specific ``C'' routines}] Section \ref{sec:julmod}
presents a Mathematica program characterizing the information required for
stochastic simulations.
Section \ref{sec:prepstep}
presents a Mathematica program that uses the model information to generate
all the model specific ``C''  code.


\item[{\bf compiling and linking the program}] Section \ref{sec:makefile}
presents a  makefile  for creating the program from it's constituent parts.
To use it issue the command 

gmake stochRun

The make file contains the following command to link the 
program {\bf stochRun.o} with the nonlinear model {\bf julliard.o}
and the other {\bf stochProto} files
creating the program stochRun:
\begin{verbatim}
stochRun:	stochRun.o stackC.o  myNewt.o \
			ma50ad.o stochProto.o julliard.o
		f77 -o stochRun -O4 stochRun.o myNewt.o\
		stackC.o julliard.o ma50ad.o stochProto.o\
		$(fastSPARSELIB) $(LINKFLAGS) $(fastLAPACKLIB)
\end{verbatim}

The file julliard.o represents the only model specific program.

\item[{\bf command line switches}] \ 

stochRun -a {\em t0} -s { \em stochasticPathLength} -r {\em replications }-l {\em pathLength}  -f {\em outputFilename}
\begin{description}
\item[{\bf t0}] The index of the first data value. This determines the initial
conditions for the set of simulations.
\item[{\bf stochasticPathLength}] The length of the stochastic path. Must be at least 1.
\item[{\bf replications}] Number of draws.
\item[{\bf pathLength}] Length of ``stack'' in perfect foresight solution. Must be
greater than 0 in current implementation.
\item[{\bf outputFilename}] Output file will contain Mathematica lists of
  \begin{enumerate}
\item a vector describing model input parameters
  \begin{enumerate}
  \item     number of model equations
  \item number of lags,
  \item number of leads,
  \item certainty equivalence computation path length,
  \item index into data for initial conditions,
  \item stochastic path length,
  \item the number of draws or  replications
  \end{enumerate}
  \item the {\bf failedQ} vector,
  \item  the shockPermutationIndices vector,
  \item  simulation results
\item the model data used in the program. 
\item the model shocks used in the program.
  \end{enumerate}
\end{description}
\end{description}



\section{main}
\label{sec:main}

The example program runs stocastic simulations using the
5 equation model presented by Julliard in presenting Stack\footnote{
The numbers just inside the $<$  $>$ provide
the page numbers of the component's definition.}

@o stochRun.c -d
@{
#define genericPeriodicPointGuesser julModPeriodicPointGuesser
#define genericDerivative julModDerivative
#define generic julMod
#define genericData julModData
#define genericShocks julModShocks
#define genericNLAGS 1
#define genericNLEADS 5
#define genericNEQS 5
#define SHOCKS 30
#define DATA 50



@<defines and includes@>


/*unsigned int  dtime(double * userSystemTime);*/

#include "stochProto.h"

int main(int argc, char * argv[])
{
@<main variable declarations@>
@<main storage allocations determined at compile time@>
@<process command line@>
@<main scalar variable initializations@>
@<main storage allocations determined at run time@>
@<obtain fixed point for terminal constraint@>
@<obtain shocks and data@>
@<generate shock indices@>
@<carryout stochastic sims@>
@<report results@>
@<main storage deallocations@> 
return(0);
}
@| main
@}


\subsection{The Model}
\label{sec:julliard}
Section \ref{sec:def} in the Appendix  presents the Mathematica input
describing the model.
The stack code requires two routines.
The first computes $f(x)$ the second computes $\frac{\partial f}{\partial x}$.
They both must return the result in Compressed Sparse Row format(CSR) format\cite{saad94}.
A Mathematica program generated the code characterizing the Julliard model.
The programs reside in the file 
julliard.c.


@d defines and includes
@{
/*#include <string>*/
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include "useSparseAMA.h"
#include "stackC.h"
void generic(double * xvec,double * pvec,
double * alhs,
unsigned int * jalhs,
unsigned int * ialhs
);
void genericDerivative(double * xvec,double * pvec,
double * alhs,
unsigned int * jalhs,
unsigned int * ialhs);

void genericData(unsigned int t,double * vectorOfVals);
void genericShocks(unsigned int t,double * vectorOfVals);
void genericPeriodicPointGuesser
(double * parameters,unsigned int period,
	double *);

void generateDraws(unsigned int t0Index,unsigned int tfIndex,unsigned int replications,unsigned int shocksAvailable,
unsigned int * iarray);
void altComputeAsymptoticQMatrix(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
void (* func)(),void (* dfunc)(),double * params,
double canadaFP[],unsigned int * pthLngth,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
double * qMat,unsigned int * qMatj,unsigned int * qMati,
unsigned int * ierr
);
void stochSim(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,
unsigned int * replications,
unsigned int * t0,unsigned int * tf,unsigned int * permVecs,
double * shockTable,unsigned int * shocksAvailable,
double * dataTable,unsigned int * dataAvailable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,
double x[],
unsigned int *failedQ);
void fPrintMathDbl(FILE * file,unsigned int length,double * matrix,char *  matrixName);

void fPrintMathInt(FILE * file,unsigned int length,unsigned int * matrix,char *  matrixName);



@}

The user may specify the maximum number of non zero elements
for the sparse matrices,
and the maximum length of the simulation path and the maximum number
or replications.
The driver program uses these to allocate  space for the computations.

@d defines and includes
@{
#define MAXELEMENTS 20000
//#define PATHLENGTH 25
//#define REPLICATIONS 5000


@}



\subsection{Main Variable Declarations}
\label{sec:mainVariables}

The algorithms access the dimensions of the 
nonlinear system through ``FORTRAN'' style
integer variables.
@d main variable declarations
@{
char * paramNamesArray[] = {};

double  parameters[0]={};
char * namesArray[] =  {"aDummy", "cc", "kk", "theta"};
unsigned int numberOfParameters=0;
unsigned int numberOfData=DATA;
unsigned int numberOfShocks=SHOCKS;
unsigned int numberOfEquations=genericNEQS;
unsigned int lags[1]={genericNLAGS};
unsigned int leads[1]={genericNLEADS};
unsigned int pathLength[1]={PATHLENGTH};
unsigned int t0[1]={0};
unsigned int tf[1]={0};
unsigned int replications[1]={1};
double totalTime[1];
unsigned int numDATA=500;
unsigned int numSHOCKS=500;
double userSystemTime[2];
/*unsigned int shockIndex[1];*/
/*void * calloc(unsigned num,unsigned int amt);*/

double *julModDataVals=(double *)calloc(numberOfEquations*numDATA,sizeof(double));
unsigned int i;
for(i=0;i<numDATA;i++){julModData(i,julModDataVals+(i*numberOfEquations));}

double *julModShockVals=(double *)calloc(numberOfEquations*numSHOCKS,sizeof(double));
for(i=0;i<numSHOCKS;i++){julModShocks(i,julModShockVals+(i*numberOfEquations));}


@| numberOfEquations lags leads
@}


@d main scalar variable initializations
@{
/*unsigned int  dtime(double * userSystemTime);*/

/**totalTime=dtime(userSystemTime);*/
printf("initializing variables\n totalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));

@}

The main routine explicitly frees these calloc's.








\subsection{Using the Q Matrix to Compute Stochastic Simulations}
\label{sec:qmatrix}


The following code obtains the asymptotic Q matrix, 
initializes the first path to the fixed point value, 
and calls the
stochSim routine.


@d carryout stochastic sims
@{



altComputeAsymptoticQMatrix(
&numberOfEquations,lags,leads,
generic,genericDerivative,genericParam,
genericFP,pathLength,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
AMqMatrix,AMqMatrixj,AMqMatrixi,
chk
);
/*
if(chk[0])
{  printf("problems computing  Q matrix\n");return(1);}
else {printf("computed Q matrix\n");}
*/
printf("computed Q matrix\n");
for(i=0;i<*pathLength;i++){
genericPeriodicPointGuesser(genericParam,1,
genericPathQ+(i *genericNEQS));}
/**totalTime=dtime(userSystemTime);*/
printf("after computing Q matrix\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));
printf("using q matrix\n");
cPrintSparse(5,AMqMatrix,AMqMatrixj,AMqMatrixi);


stochSim(&numberOfEquations,lags,leads,pathLength,
generic,genericDerivative,genericParam,
replications,t0,tf,genericPermVec,
genericShocksArray,&numberOfShocks,
genericDataArray,&numberOfData,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,AMqMatrix,AMqMatrixj,AMqMatrixi,
genericFP,
genericPathQ,
failedQ);
@|setgm sval
@}

@d obtain shocks and data
@{


for(i=0;i<DATA;i++){genericData(i,genericDataArray+(i*genericNEQS));}
for(i=0;i<SHOCKS;i++){genericShocks(i,genericShocksArray+(i*genericNEQS));}
@}
@d generate shock indices
@{

printf("generating perm vec\n");
 generateDraws(1,(*stochasticPathLength),(*replications),SHOCKS,genericPermVec);
printf("done generating perm vec\n");



@}

@o stochRunHide.c
@{
void fPrintMathDbl(FILE * file,unsigned int length,double * matrix,char *  matrixName)
{
unsigned int i;
fprintf(file,"%s={",matrixName);
for(i=0;(i<length-1);i++){
fprintf(file,"%30.20f,",matrix[i]);}
fprintf(file,"%30.20f};\n",matrix[length-1]);
}
void fPrintMathInt(FILE * file,unsigned int length,unsigned int * matrix,char *  matrixName)
{
unsigned int i;
fprintf(file,"%s={",matrixName);
for(i=0;(i<length-1);i++){
fprintf(file,"%d,",matrix[i]);}
fprintf(file,"%d};\n",matrix[length-1]);
}
@}



@d report results
@{
FILE * outFile;
outFile=fopen(flnm,"w");

printf("saving values for variable in file named %s\n",flnm);
fprintf(outFile,"genericRunParams={%d,%d,%d,%d,%d,%d,%d};\n",
    genericNEQS,genericNLAGS,genericNLEADS,
     *pathLength,*t0,*stochasticPathLength,*replications);
fPrintMathInt(outFile,*replications,failedQ,"genericFailedQ");
fPrintMathInt(outFile,*replications * (*stochasticPathLength),
      genericPermVec,"genericPermVec");
fPrintMathDbl(outFile,(*replications * genericNEQS*(*stochasticPathLength+genericNLAGS)),
      genericPathQ,"genericResults");
fPrintMathDbl(outFile,(genericNEQS*(DATA)),genericDataArray,"genericDataArray");
fPrintMathDbl(outFile,(genericNEQS*(SHOCKS)),genericShocksArray,"genericShocksArray");
     fclose(outFile);


@}


@d main variable declarations
@{
double atof();
@}



@d process command line
@{
/*vbl=0;*//*hack so that if shock irrelevant variable if no other variables shocked*/
*pathLength=1;
*replications=1;
*t0=genericNLAGS+1;
*stochasticPathLength=1;
printf("default values:(pathLength=%u,replications=%u,t0=%u,stochasticPathLength=%u)\n",*pathLength,*replications,*t0,*stochasticPathLength);

processCommandLine(argc,argv,namesArray,numberOfEquations,
paramNamesArray,numberOfParameters,parameters,
	julModDataVals,numDATA,numSHOCKS,
	pathLength,replications,t0,stochasticPathLength,
intControlParameters,doubleControlParameters,flnm);


@}

@d find variable name match
@{
i=0;
while( (i <genericNEQS) &&(strcmp(argv[2],genericNamesArray[i])))i++;
if(i==genericNEQS){
vbl=0;/*shock something that's irrelevant*/
printf("i don't know the variable %s: ignoring this variable value pair\n",
argv[2]);} else {vbl = i;}
@}

The limits the size of sparse matrices to the amount specified on the argument list.
@d main variable declarations
@{
unsigned int maxNumberElements[1]={MAXELEMENTS};
@}


The routines will need $L(\tau+\theta+1)$ doubles to hold the fixed point during calculation.
@d main variable declarations
@{

unsigned int stochasticPathLength[1]={1};
unsigned int * genericPermVec;
double * genericShocksArray;
double * genericDataArray;
double * genericFP;
double genericParam[2]={0.5,0.6};
double * genericPathQ;
double **fmats;unsigned int  **fmatsj;unsigned int  **fmatsi;
double **smats;unsigned int  **smatsj;unsigned int  **smatsi;
static char flnm[50] = "stochOut.m";


@<define names array@> 
@}
@d define names array
@{
/*char * genericNamesArray[] =  
{"ey","pdot","rr","rs","y"};*/

@}
@d main storage allocations determined at run time
@{
failedQ=(unsigned int *)calloc(*replications,sizeof(unsigned int));
for(i=0;i<*replications;i++)failedQ[i]=0;
*tf=(*t0)+(*stochasticPathLength)-1;

genericPermVec=(unsigned int *)calloc(
     (*stochasticPathLength)*(*replications),sizeof(unsigned int));
genericPathQ=(double *)calloc(
    *replications*
    genericNEQS*(genericNLAGS+genericNLEADS+(*pathLength)+(*stochasticPathLength)),
sizeof(double));
double ** ptrToPtrToDouble = NULL;
unsigned int ** ptrToPtrToInt = NULL;
fmats =(double **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToDouble));
fmatsj =(unsigned int **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToInt));
fmatsi =(unsigned int **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToInt));
smats =(double **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToDouble));
smatsj =(unsigned int **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToInt));
smatsi =(unsigned int **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToInt));
for(i=0;i<(*pathLength)+genericNLAGS+1;i++){
fmats[i] =(double *)calloc(MAXELEMENTS,sizeof(double));
fmatsj[i] =(unsigned int *)calloc(MAXELEMENTS,sizeof(unsigned int));
fmatsi[i] =(unsigned int *)calloc(
     genericNEQS*(genericNLAGS+genericNLEADS)+1,sizeof(unsigned int));
smats[i] =(double *)calloc(MAXELEMENTS,sizeof(double));
smatsj[i] =(unsigned int *)calloc(MAXELEMENTS,sizeof(unsigned int));
smatsi[i] =(unsigned int *)calloc(
     genericNEQS*(genericNLAGS+genericNLEADS)+1,sizeof(unsigned int));
}


@}

@d main storage allocations determined at compile time
@{
genericShocksArray=(double *)calloc(
     genericNEQS*(SHOCKS),sizeof(double));
genericDataArray=(double *)calloc(
     genericNEQS*(DATA),sizeof(double));
genericFP=(double *)calloc(
     genericNEQS*(genericNLAGS+genericNLEADS+1),sizeof(double));
AMqMatrix=(double *)
   calloc(MAXELEMENTS,sizeof(double));
AMqMatrixj=(unsigned int *)
   calloc(MAXELEMENTS,sizeof(unsigned int));
AMqMatrixi=(unsigned int *)
   calloc((genericNEQS*(genericNLEADS+genericNLAGS)),
        sizeof(unsigned int));
@}


@d main storage deallocations
@{
free(failedQ);
free(AMqMatrix);
free(AMqMatrixj);
free(AMqMatrixi);
free(genericShocksArray);
free(genericDataArray);
free(genericPermVec);
free(genericFP);
free(genericPathQ);
for(i=0;i<(*pathLength)+genericNLAGS+1;i++){
free(smats[i]);
free(smatsj[i]);
free(smatsi[i]);
free(fmats[i]);
free(fmatsj[i]);
free(fmatsi[i]);
}
free(fmats);
free(fmatsj);
free(fmatsi);
free(smats);
free(smatsj);
free(smatsi);
@}


chk returns 0 for success 1 if the routine may not have converged.

@d main variable declarations
@{
unsigned int * failedQ;
unsigned int chk[1]={0};
@}
@d obtain fixed point for terminal constraint
@{
printf("$Id: stochRun.w,v 1.6 2000/12/06 14:53:34 m1gsa00 Exp m1gsa00 $\n");

genericPeriodicPointGuesser(genericParam,1,genericFP);

unsigned int * inIntControl={0};
double * inDoubleControl={0};
unsigned int * outIntControl={0};
double * outDoubleControl={0};
double * linearizationPoint={0};


FPnewt(&numberOfEquations,lags,leads,
generic,genericDerivative,genericParam,
genericFP,linearizationPoint,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
chk ,inIntControl,inDoubleControl,outIntControl,outDoubleControl);

/*
if(chk[0])
{  printf("problems computing  FP solution\n");return(1);}
else {printf("computed FP solution\n");}
*/

printf("computed FP solution\n");


@}


@d main variable declarations
@{
/*unsigned int * hColumns;*/
@}

@d main variable declarations
@{
double * AMqMatrix;
unsigned int * AMqMatrixj;
unsigned int * AMqMatrixi;
/*double * asymptoticLinearization;*/
@|
asymptoticLinearization AMqMatrix
@}

@d main storage allocations determined at compile time
@{
/**totalTime=dtime(userSystemTime);*/
printf("Hello World!!  after compile time determined storage allocations\n totalTime=%f,userSystemTime=%f,systemTime=%f\n",
     *totalTime,*userSystemTime,*(userSystemTime+1));

@}
Placing times at the end of each scrap.


@d obtain fixed point for terminal constraint
@{
/**totalTime=dtime(userSystemTime);*/
printf("after fixed point computation\n totalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));
@}



@d carryout stochastic sims
@{
/**totalTime=dtime(userSystemTime);*/
printf("after using Q matrix\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));

@}



@o stochSimsUnitTests.c -d
@{

#include <stdio.h>
#include <string.h>
#include "CUnit/Basic.h"
#include<stdlib.h>
#include "useSparseAMA.h"
#define genericPeriodicPointGuesser rbcExamplePeriodicPointGuesser
#define genericDerivative rbcExampleDerivative
#define generic rbcExample
#define genericData rbcExampleData
#define genericShocks rbcExampleShocks
#define genericNLAGS 1
#define genericNLEADS 1
#define genericNEQS 4
#define SHOCKS 30
#define DATA 50


int init_genericsuite1(void)
{
return(0);
}
int clean_genericsuites(void)
{
return(0);
}



/* The main() function for setting up and running the tests.
 * Returns a CUE_SUCCESS on successful running, another
 * CUnit error code on failure.
 */
@<generic defines and includes@>

int main(int argc, char * argv[])
{
@<generic main variable declarations@>
@<generic main storage allocations determined at compile time@>
@<process command line@>
@<generic main scalar variable initializations@>

@<generic main storage allocations determined at run time@>
@<generic obtain fixed point for terminal constraint@>
@<generic obtain shocks and data@>
@<generic generate shock indices@>
@<generic carryout stochastic sims@>
@<generic report results@>
@<generic main storage deallocations@> 



   CU_pSuite pSuite = NULL;

   /* initialize the CUnit test registry */
   if (CUE_SUCCESS != CU_initialize_registry())
      return CU_get_error();



  /* add another suite to the registry */
   pSuite = CU_add_suite("Genericsuite_1", init_genericsuite1, clean_genericsuites);
   if (NULL == pSuite) {
     CU_cleanup_registry();
    return CU_get_error();
  }


   /* add the tests to the suite 
   if ((NULL == CU_add_test(pSuite, "test of oneEquationZeroLead()", oneEquationZeroLead)))
   {
      CU_cleanup_registry();
      return CU_get_error();
   }

*/


   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}

@}
@d generic obtain fixed point for terminal constraint
@{
printf("$Id: stochRun.w,v 1.6 2000/12/06 14:53:34 m1gsa00 Exp m1gsa00 $\n");

genericPeriodicPointGuesser(genericParam,1,genericFP);

unsigned int * inIntControl={0};
double * inDoubleControl={0};
unsigned int * outIntControl={0};
double * outDoubleControl={0};
double * linearizationPoint={0};

FPnewt(numberOfEquations,lags,leads,
generic,genericDerivative,genericParam,
genericFP,linearizationPoint,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
chk,inIntControl,inDoubleControl,outIntControl,outDoubleControl);

/*
if(chk[0])
{  printf("problems computing  FP solution\n");return(1);}
else {printf("computed FP solution\n");}
*/

printf("computed FP solution\n");


@}
@d generic obtain shocks and data
@{


for(i=0;i<DATA;i++){genericData(i,genericDataArray+(i*genericNEQS));}
for(i=0;i<SHOCKS;i++){genericShocks(i,genericShocksArray+(i*genericNEQS));}
@}

@d generic obtain fixed point for terminal constraint
@{
/**totalTime=dtime(userSystemTime);*/
printf("after fixed point computation\n totalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));
@}

@d generic generate shock indices
@{

printf("generating perm vec\n");
 generateDraws(1,(*stochasticPathLength),(*replications),SHOCKS,genericPermVec);
printf("done generating perm vec\n");



@}
@d generic carryout stochastic sims
@{



altComputeAsymptoticQMatrix(
numberOfEquations,lags,leads,
generic,genericDerivative,genericParam,
genericFP,pathLength,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
AMqMatrix,AMqMatrixj,AMqMatrixi,
chk
);
/*
if(chk[0])
{  printf("problems computing  Q matrix\n");return(1);}
else {printf("computed Q matrix\n");}
*/
printf("computed Q matrix\n");
for(i=0;i<*pathLength;i++){
genericPeriodicPointGuesser(genericParam,1,
genericPathQ+(i *genericNEQS));}
/**totalTime=dtime(userSystemTime);*/
printf("after computing Q matrix\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));
printf("using q matrix\n");
cPrintSparse(5,AMqMatrix,AMqMatrixj,AMqMatrixi);


stochSim(numberOfEquations,lags,leads,pathLength,
generic,genericDerivative,genericParam,
replications,t0,tf,genericPermVec,
genericShocksArray,numberOfShocks,
genericDataArray,numberOfData,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,AMqMatrix,AMqMatrixj,AMqMatrixi,
genericFP,
genericPathQ,
failedQ);
@|setgm sval
@}
@d generic carryout stochastic sims
@{
/**totalTime=dtime(userSystemTime);*/
printf("after using Q matrix\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));

@}
@d generic main storage deallocations
@{
free(failedQ);
free(AMqMatrix);
free(AMqMatrixj);
free(AMqMatrixi);
free(genericShocksArray);
free(genericDataArray);
free(genericPermVec);
free(genericFP);
free(genericPathQ);
for(i=0;i<(*pathLength)+genericNLAGS+1;i++){
free(smats[i]);
free(smatsj[i]);
free(smatsi[i]);
free(fmats[i]);
free(fmatsj[i]);
free(fmatsi[i]);
}
free(fmats);
free(fmatsj);
free(fmatsi);
free(smats);
free(smatsj);
free(smatsi);
@}


@d generic report results
@{
FILE outFile=fopen(flnm,"w");
printf("saving values for variable in file named %s\n",flnm);
fprintf(outFile,"genericRunParams={%d,%d,%d,%d,%d,%d,%d};\n",
    genericNEQS,genericNLAGS,genericNLEADS,
     *pathLength,*t0,*stochasticPathLength,*replications);
fPrintMathInt(outFile,*replications,failedQ,"genericFailedQ");
fPrintMathInt(outFile,*replications * (*stochasticPathLength),
      genericPermVec,"genericPermVec");
fPrintMathDbl(outFile,(*replications * genericNEQS*(*stochasticPathLength+genericNLAGS)),
      genericPathQ,"genericResults");
fPrintMathDbl(outFile,(genericNEQS*(DATA)),genericDataArray,"genericDataArray");
fPrintMathDbl(outFile,(genericNEQS*(SHOCKS)),genericShocksArray,"genericShocksArray");
     fclose(outFile);


@}


@d generic main storage allocations determined at run time
@{
failedQ=(unsigned int *)calloc(*replications,sizeof(unsigned int));
for(i=0;i<*replications;i++)failedQ[i]=0;
*tf=(*t0)+(*stochasticPathLength)-1;

genericPermVec=(unsigned int *)calloc(
     (*stochasticPathLength)*(*replications),sizeof(unsigned int));
genericPathQ=(double *)calloc(
    *replications*
    genericNEQS*(genericNLAGS+genericNLEADS+(*pathLength)+(*stochasticPathLength)),
sizeof(double));
double ** ptrToPtrToDouble = NULL;
unsigned int ** ptrToPtrToInt = NULL;
fmats =(double **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToDouble));
fmatsj =(unsigned int **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToInt));
fmatsi =(unsigned int **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToInt));
smats =(double **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToDouble));
smatsj =(unsigned int **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToInt));
smatsi =(unsigned int **)calloc((*pathLength)+genericNLAGS+1,sizeof(ptrToPtrToInt));
for(i=0;i<(*pathLength)+genericNLAGS+1;i++){
fmats[i] =(double *)calloc(MAXELEMENTS,sizeof(double));
fmatsj[i] =(unsigned int *)calloc(MAXELEMENTS,sizeof(unsigned int));
fmatsi[i] =(unsigned int *)calloc(
     genericNEQS*(genericNLAGS+genericNLEADS)+1,sizeof(unsigned int));
smats[i] =(double *)calloc(MAXELEMENTS,sizeof(double));
smatsj[i] =(unsigned int *)calloc(MAXELEMENTS,sizeof(unsigned int));
smatsi[i] =(unsigned int *)calloc(
     genericNEQS*(genericNLAGS+genericNLEADS)+1,sizeof(unsigned int));
}


@}


@d generic main scalar variable initializations
@{
/*unsigned int  dtime(double * userSystemTime);*/

/**totalTime=dtime(userSystemTime);*/
printf("initializing variables\n totalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));

@}

@d generic defines and includes
@{
/*#include <string>*/
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include "useSparseAMA.h"
#include "stackC.h"
void generic(double * xvec,double * pvec,
double * alhs,
unsigned int * jalhs,
unsigned int * ialhs
);
void genericDerivative(double * xvec,double * pvec,
double * alhs,
unsigned int * jalhs,
unsigned int * ialhs);
void genericData(unsigned int t,double * vectorOfVals);
void genericShocks(unsigned int t,double * vectorOfVals);
void genericPeriodicPointGuesser
(double * parameters,unsigned int period,
	double *);

void generateDraws(unsigned int t0Index,unsigned int tfIndex,unsigned int replications,unsigned int shocksAvailable,
unsigned int * iarray);
void altComputeAsymptoticQMatrix(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
void (* func)(),void (* dfunc)(),double * params,
double canadaFP[],unsigned int * pthLngth,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
double * qMat,unsigned int * qMatj,unsigned int * qMati,
unsigned int * ierr
);
void stochSim(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,
unsigned int * replications,
unsigned int * t0,unsigned int * tf,unsigned int * permVecs,
double * shockTable,unsigned int * shocksAvailable,
double * dataTable,unsigned int * dataAvailable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,
double x[],
unsigned int *failedQ);
void fPrintMathDbl(FILE * file,unsigned int length,double * matrix,char *  matrixName);

void fPrintMathInt(FILE * file,unsigned int length,unsigned int * matrix,char *  matrixName);



@}

@d generic defines and includes
@{
#define MAXELEMENTS 20000
#define PATHLENGTH 25
#define REPLICATIONS 5000


@}




@d generic main storage allocations determined at compile time
@{
/**totalTime=dtime(userSystemTime);*/
printf("Hello World!!  after compile time determined storage allocations\n totalTime=%f,userSystemTime=%f,systemTime=%f\n",
     *totalTime,*userSystemTime,*(userSystemTime+1));

@}

@d generic main variable declarations
@{

unsigned int numberOfData[1]={DATA};
unsigned int numberOfShocks[1]={SHOCKS};
unsigned int numberOfEquations[1]={genericNEQS};
unsigned int lags[1]={genericNLAGS};
unsigned int leads[1]={genericNLEADS};
unsigned int pathLength[1]={PATHLENGTH};
unsigned int t0[1]={0};
unsigned int tf[1]={0};
unsigned int replications[1]={1};
double totalTime[1];
double userSystemTime[2];
/*unsigned int shockIndex[1];*/
/*void * calloc(unsigned num,unsigned int amt);*/
@| numberOfEquations lags leads
@}


@d generic main variable declarations
@{
double atof();
unsigned int pl;
/*unsigned int vbl;*/
/*double val;*/
@}
@d generic main variable declarations
@{
unsigned int maxNumberElements[1]={MAXELEMENTS};
@}

@d generic main variable declarations
@{

FILE * outFile;
unsigned int stochasticPathLength[1]={1};
unsigned int * genericPermVec;
double * genericShocksArray;
double * genericDataArray;
double * genericFP;
double genericParam[2]={0.5,0.6};
double * genericPathQ;
double **fmats;unsigned int  **fmatsj;unsigned int  **fmatsi;
double **smats;unsigned int  **smatsj;unsigned int  **smatsi;
static char flnm[50] = "stochOut.m";
unsigned int i/*,j*/;

@<define names array@> 
@}
@d generic main variable declarations
@{
unsigned int * failedQ;
unsigned int chk[1]={0};
@}

@d generic main storage allocations determined at compile time
@{
genericShocksArray=(double *)calloc(
     genericNEQS*(SHOCKS),sizeof(double));
genericDataArray=(double *)calloc(
     genericNEQS*(DATA),sizeof(double));
genericFP=(double *)calloc(
     genericNEQS*(genericNLAGS+genericNLEADS+1),sizeof(double));
AMqMatrix=(double *)
   calloc(MAXELEMENTS,sizeof(double));
AMqMatrixj=(unsigned int *)
   calloc(MAXELEMENTS,sizeof(unsigned int));
AMqMatrixi=(unsigned int *)
   calloc((genericNEQS*(genericNLEADS+genericNLAGS)),
        sizeof(unsigned int));
@}



@d generic main variable declarations
@{
double * AMqMatrix;
unsigned int * AMqMatrixj;
unsigned int * AMqMatrixi;
/*double * asymptoticLinearization;*/
@|
asymptoticLinearization AMqMatrix
@}

@o stochSimsUnitTests.c -d
@{



@}




\appendix


\section{Non Linear Model Definition} 
\label{sec:def}

\subsection{Generic Model}
\label{sec:genericmod}


\listinginput{1}{/mq/home4/m1gsa00/aim/nonLinearModels/julliard/julliardModel.m}

\subsection{Mathematica Preparation Step}
\label{sec:prepstep}


\listinginput{1}{/mq/home4/m1gsa00/aim/nonLinearModels/julliard/prepareJulliardModel.m}

\section{Makefile}
\label{sec:makefile}

@o makefile -t
@{
#identify operating system
UNAME= $(shell uname)
NUWEBFLAGS = -t
SPAMADIR=../sparseAMA

ifeq ($(UNAME),Linux)
#compilers
CC = gcc
FCFLAGS = -c -O2 -I $(SPAMADIR)/src/main/include   -I /msu/res5/software/myUsr/include/ -I/msu/res5/software/myUsr/include/ -I /msu/res1/Software/matio-1.5.1/src
FCFLAGS = -c -g -Wall -I $(SPAMADIR)/src/main/include   -I /msu/res5/software/myUsr/include -I/msu/res5/software/myUsr/include/ -I /msu/res1/Software/matio-1.5.1/src
#lapack
LAPACKLIBS=   -L /msu/res5/software/ARPACK96forCluster -larpack_linux -L/msu/res5/software/lapackGithubForCluster -llapack -lrefblas
CUNITLIBS= -L/msu/res5/software/myUsr/lib/ -l cunit
MATIOLIBS= -L/msu/res1/Software/matio-1.5.1/src/.libs/ -lmatio  -lhdf5

endif

ifeq ($(UNAME),Darwin)
#compilers
CC = gcc-6
FCFLAGS = -c -O2 -I$(SPAMADIR)/src/main/include   -I /Users/garyanderson/myUsr/include/ -I/Users/garyanderson/myUsr/include/\
-I /usr/local/Cellar/libmatio/1.5.10/include
FCFLAGS = -c -Wall -g -I $(SPAMADIR)/src/main/include   -I /Users/garyanderson/myUsr/include/ -I/Users/garyanderson/myUsr/include/\
-I /usr/local/Cellar/libmatio/1.5.10/include
#lapack
LAPACKLIBS=  -L /Users/garyanderson/ARPACK96/  -larpack_MACOS -L /Users/garyanderson/lapack-release/ -llapack -lrefblas
CUNITLIBS= -L /Users/garyanderson/myUsr/lib -l cunit
MATIOLIBS= -L/usr/local/Cellar/libmatio/1.5.10/lib -lmatio 

endif


#compilers
FC = gfortran
SPARSEAMALIB= -L../sparseAMA -lsparseAMA
STOCHSIMSLIB= -L./ -lstochSims
.w.c:
	nuweb $(NUWEBFLAGS) $*

.w.o:
	make $*.c
	make $*.o

.c.o:
	$(CC) $(FCFLAGS) -c $*.c

.f.o:
	$(FC) $(FCFLAGS) -c $*.f

.PHONY: Build

Build: stochRun runrbcTryC
#	$(FC) -o stochRun -g  stochRun.o juillard.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS) $(CUNITLIBS) $(MATIOLIBS)
#	$(FC) -o runrbcTryC -g  runrbcTryC.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS) $(CUNITLIBS) $(MATIOLIBS)
#	$(FC) -o stochSimsUnitTests -g  stochSimsUnitTests.o juillard.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS) $(CUNITLIBS) $(MATIOLIBS)

myNewt.o:			 stackC.w
		nuweb $(NUWEBFLAGS)  stackC.w
	$(CC) $(FCFLAGS) -c myNewt.c

juillard.o:	juillard.c
	$(CC) $(FCFLAGS) -c juillard.c


stochRun:	stochRun.o  juillard.o libstochSims.a
	$(FC) -o stochRun -g  stochRun.o juillard.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS)  $(CUNITLIBS) $(MATIOLIBS)

stochSimsUnitTests:	stochSimsUnitTests.o  rbcTryC.o rbcTryCDrv.o rbcTryCData.o rbcTryCShocks.o rbcTryCSupport.o libstochSims.a
	$(FC) -o stochSimsUnitTests -g  stochSimsUnitTests.o  rbcTryC.o rbcTryCDrv.o rbcTryCData.o rbcTryCShocks.o rbcTryCSupport.o  $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS)  $(CUNITLIBS) $(MATIOLIBS)

runrbcTryC:	runrbcTryC.o  rbcTryC.o rbcTryCDrv.o rbcTryCData.o rbcTryCShocks.o rbcTryCSupport.o libstochSims.a
	$(FC) -o runrbcTryC -g  runrbcTryC.o  rbcTryC.o rbcTryCDrv.o rbcTryCData.o rbcTryCShocks.o rbcTryCSupport.o  $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS)  $(CUNITLIBS) $(MATIOLIBS)

libstochSims.a:	myNewt.o \
		stackC.o stochProto.o ranlib.o
	ar -cvq libstochSims.a myNewt.o \
		stackC.o stochProto.o ranlib.o


clean: 
	rm -f *.o stochRun stochSimsUnitTests libstochSims.a

@}


\listinginput{1}{makeStochTry}



\subsection{Sample Output}
\label{sec:sample}

\begin{verbatim}
\end{verbatim}

\section{Index}
\label{sec:index}




\subsection{Files}
\label{sec:files}



@f


\subsection{Macros}
\label{sec:macros}


@m



\subsection{Names}
\label{sec:names}




@u



\bibliographystyle{authordate2}
\bibliography{files,anderson}

\end{document}
