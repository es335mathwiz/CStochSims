%	$Id: stochRun.w,v 1.6 2000/12/06 14:53:34 m1gsa00 Exp m1gsa00 $	
\documentclass[html]{article}

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
@<defines and includes@>
int  dtime(double * userSystemTime);

void  cfree(void * ptr){free(ptr);}
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
#include "sparseAMA.h"
#include "stackC.h"
#define julNLAGS 1
#define julNLEADS 5
#define julNEQS 5
#define SHOCKS 30
#define DATA 50
void julMod(double * xvec,double * pvec,
double * alhs,
int * jalhs,
int * ialhs
);
void julModDerivative(double * xvec,double * pvec,
double * alhs,
int * jalhs,
int * ialhs);
void FPnewt(int * numberOfEquations,int * lags, int * leads,
void (* func)(),void (* dfunc)(),double * params,
double x[],
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,
int *check);
void julModData(int t,double * vectorOfVals);
void julModShocks(int t,double * vectorOfVals);
void julModPeriodicPointGuesser
(double * parameters,int period,
	double *);

void generateDraws(int t0Index,int tfIndex,int replications,int shocksAvailable,
int * iarray);
void altComputeAsymptoticQMatrix(
int * numberOfEquations,int * lags, int * leads,
void (* func)(),void (* dfunc)(),double * params,
double canadaFP[],int * pthLngth,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,
double * qMat,int * qMatj,int * qMati,
int * ierr
);
void stochSim(
int * numberOfEquations,int * lags, int * leads,int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,
int * replications,
int * t0,int * tf,int * permVecs,
double * shockTable,int * shocksAvailable,
double * dataTable,int * dataAvailable,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,double * qMat,int * qMatj,int * qMati,
double * fixedPoint,
double x[],
int *failedQ);
void fPrintMathDbl(FILE * file,int length,double * matrix,char *  matrixName);

void fPrintMathInt(FILE * file,int length,int * matrix,char *  matrixName);



@}

The user may specify the maximum number of non zero elements
for the sparse matrices,
and the maximum length of the simulation path and the maximum number
or replications.
The driver program uses these to allocate  space for the computations.

@d defines and includes
@{
#define MAXELEMENTS 20000
#define PATHLENGTH 25
#define REPLICATIONS 5000


@}



\subsection{Main Variable Declarations}
\label{sec:mainVariables}

The algorithms access the dimensions of the 
nonlinear system through ``FORTRAN'' style
integer variables.
@d main variable declarations
@{

int numberOfData[1]={DATA};
int numberOfShocks[1]={SHOCKS};
int numberOfEquations[1]={julNEQS};
int lags[1]={julNLAGS};
int leads[1]={julNLEADS};
int pathLength[1]={PATHLENGTH};
int t0[1]={0};
int tf[1]={0};
int replications[1]={1};
double totalTime[1];
double userSystemTime[2];
/*int shockIndex[1];*/
/*void * calloc(unsigned num,int amt);*/
@| numberOfEquations lags leads
@}


@d main scalar variable initializations
@{
/*int  dtime(double * userSystemTime);*/

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
numberOfEquations,lags,leads,
julMod,julModDerivative,julParam,
julliardFP,pathLength,
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
julModPeriodicPointGuesser(julParam,1,
julliardPathQ+(i *julNEQS));}
/**totalTime=dtime(userSystemTime);*/
printf("after computing Q matrix\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));
printf("using q matrix\n");
cPrintSparse(5,AMqMatrix,AMqMatrixj,AMqMatrixi);


stochSim(numberOfEquations,lags,leads,pathLength,
julMod,julModDerivative,julParam,
replications,t0,tf,julliardPermVec,
julliardShocks,numberOfShocks,
julliardData,numberOfData,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,AMqMatrix,AMqMatrixj,AMqMatrixi,
julliardFP,
julliardPathQ,
failedQ);
@|setgm sval
@}

@d obtain shocks and data
@{


for(i=0;i<DATA;i++){julModData(i,julliardData+(i*julNEQS));}
for(i=0;i<SHOCKS;i++){julModShocks(i,julliardShocks+(i*julNEQS));}
@}
@d generate shock indices
@{

printf("generating perm vec\n");
 generateDraws(1,(stochasticPathLength),(*replications),SHOCKS,julliardPermVec);
printf("done generating perm vec\n");



@}




@o stochRun.c
@{
void fPrintMathDbl(FILE * file,int length,double * matrix,char *  matrixName)
{
int i;
fprintf(file,"%s={",matrixName);
for(i=0;(i<length-1);i++){
fprintf(file,"%30.20f,",matrix[i]);}
fprintf(file,"%30.20f};\n",matrix[length-1]);
}
void fPrintMathInt(FILE * file,int length,int * matrix,char *  matrixName)
{
int i;
fprintf(file,"%s={",matrixName);
for(i=0;(i<length-1);i++){
fprintf(file,"%d,",matrix[i]);}
fprintf(file,"%d};\n",matrix[length-1]);
}
@}



@d report results
@{

printf("saving values for variable in file named %s\n",flnm);
fprintf(outFile,"julModRunParams={%d,%d,%d,%d,%d,%d,%d};\n",
    julNEQS,julNLAGS,julNLEADS,
     *pathLength,*t0,stochasticPathLength,*replications);
fPrintMathInt(outFile,*replications,failedQ,"julModFailedQ");
fPrintMathInt(outFile,*replications * (stochasticPathLength),
      julliardPermVec,"julModPermVec");
fPrintMathDbl(outFile,(*replications * julNEQS*(stochasticPathLength+julNLAGS)),
      julliardPathQ,"julModResults");
fPrintMathDbl(outFile,(julNEQS*(DATA)),julliardData,"julModData");
fPrintMathDbl(outFile,(julNEQS*(SHOCKS)),julliardShocks,"julModShocks");
     fclose(outFile);


@}


@d main variable declarations
@{
double atof();
int pl;
/*int vbl;*/
/*double val;*/
@}


@d command line pathlength 
@{
   case 'l':
	 pl=atoi(argv[2]);
     printf("got %d for path length\n",pl);
     if(pl>PATHLENGTH)
	 {
       *pathLength=PATHLENGTH;
       printf("setting pathlength to maximum=%d\n",PATHLENGTH);
       } else   if(pl<1){
       *pathLength=1;
       printf("setting pathlength to 1\n");
       } else 
     {*pathLength=pl;}
	 argc--;argv++;
	 break;
@}
@d command line replications 
@{
   case 'r':
	 pl=atoi(argv[2]);
     printf("got %d for replications\n",pl);
     if(pl>REPLICATIONS)
	 {
       *replications=REPLICATIONS;
       printf("setting repetitions to maximum=%d\n",REPLICATIONS);
       } else { *replications = pl;}
	 argc--;argv++;
	 break;
@}
@d command line t0
@{
   case 'a':
	 pl=atoi(argv[2]);
     printf("got %d for t0\n",pl);
     if(pl>DATA)
	 {
       *t0=DATA;
       printf("setting initial t0 to maximum number of data elements=%d\n",*t0);
       } else { *t0 = pl;}
     if(pl<julNLAGS)
	 {
       *t0=julNLAGS+1;
       printf("setting initial t0 to one more than number of lags=%d\n",*t0);
       } else { *t0 = pl;}
	 argc--;argv++;
	 break;
@}
@d command line stochasticPathLength
@{
   case 's':
	 stochasticPathLength=atoi(argv[2]);
     printf("got %d for stochasticPathLength\n",stochasticPathLength);
     if(stochasticPathLength<1)
	 {
       stochasticPathLength=1;
       printf("setting tf to 1\n");
       }
	 argc--;argv++;
	 break;
@}
@d command line file
@{
   case 'f':
	 strcpy(flnm,argv[2]);
	 printf("got %s for filename \n",flnm);
	 argc--;argv++;
	 argc--;argv++;
	 break;
@}
@d command line variable
@{
   case 'v':
	 val=(double)atof(argv[3]);
	 printf("got %d for vbl and %f for val\n",vbl,val);
	 argc--;argv++;
	 argc--;argv++;
	 break;
@}

@d process command line
@{
/*vbl=0;*//*hack so that if shock irrelevant variable if no other variables shocked*/
*pathLength=1;
*replications=1;
*t0=julNLAGS+1;
stochasticPathLength=1;
printf("default values:(pathLength=%d,replications=%d,t0=%d,stochasticPathLength=%d)\n",*pathLength,*replications,*t0,stochasticPathLength);


while(argc>1&&argv[1][0] == '-')
{
printf("processing command line args\n");
 switch(argv[1][1]){
@<command line pathlength@>
@<command line replications@>
@<command line t0@>
@<command line stochasticPathLength@>
@<command line file@>
 default:
   printf("%s: unknown arg %s-- not processing any more args\n",
      argv[0],argv[1]);
 }
argc--;argv++;
 }
     outFile=fopen(flnm,"w");

printf("values for run:(pathLength=%d,replications=%d,t0=%d,stochasticPathLength=%d)\n",*pathLength,*replications,*t0,stochasticPathLength);

@}

@d find variable name match
@{
i=0;
while( (i <julNEQS) &&(strcmp(argv[2],julNamesArray[i])))i++;
if(i==julNEQS){
vbl=0;/*shock something that's irrelevant*/
printf("i don't know the variable %s: ignoring this variable value pair\n",
argv[2]);} else {vbl = i;}
@}

The limits the size of sparse matrices to the amount specified on the argument list.
@d main variable declarations
@{
int maxNumberElements[1]={MAXELEMENTS};
@}


The routines will need $L(\tau+\theta+1)$ doubles to hold the fixed point during calculation.
@d main variable declarations
@{

FILE * outFile;
int stochasticPathLength=1;
int * julliardPermVec;
double * julliardShocks;
double * julliardData;
double * julliardFP;
double julParam[2]={0.5,0.6};
double * julliardPathQ;
double **fmats;int  **fmatsj;int  **fmatsi;
double **smats;int  **smatsj;int  **smatsi;
static char flnm[50] = "stochOut.m";
int i/*,j*/;

@<define names array@> 
@}
@d define names array
@{
/*char * julNamesArray[] =  
{"ey","pdot","rr","rs","y"};*/

@}
@d main storage allocations determined at run time
@{
failedQ=(int *)calloc(*replications,sizeof(int));
for(i=0;i<*replications;i++)failedQ[i]=0;
*tf=(*t0)+stochasticPathLength-1;

julliardPermVec=(int *)calloc(
     (stochasticPathLength)*(*replications),sizeof(int));
julliardPathQ=(double *)calloc(
    *replications*
    julNEQS*(julNLAGS+julNLEADS+(*pathLength)+stochasticPathLength),
sizeof(double));
double ** ptrToPtrToDouble = NULL;
fmats =(double **)calloc((*pathLength)+julNLAGS+1,sizeof(ptrToPtrToDouble));
fmatsj =(int **)calloc((*pathLength)+julNLAGS+1,sizeof(ptrToPtrToDouble));
fmatsi =(int **)calloc((*pathLength)+julNLAGS+1,sizeof(ptrToPtrToDouble));
smats =(double **)calloc((*pathLength)+julNLAGS+1,sizeof(ptrToPtrToDouble));
smatsj =(int **)calloc((*pathLength)+julNLAGS+1,sizeof(ptrToPtrToDouble));
smatsi =(int **)calloc((*pathLength)+julNLAGS+1,sizeof(ptrToPtrToDouble));
for(i=0;i<(*pathLength)+julNLAGS+1;i++){
fmats[i] =(double *)calloc(MAXELEMENTS,sizeof(double));
fmatsj[i] =(int *)calloc(MAXELEMENTS,sizeof(int));
fmatsi[i] =(int *)calloc(
     julNEQS*(julNLAGS+julNLEADS)+1,sizeof(int));
smats[i] =(double *)calloc(MAXELEMENTS,sizeof(double));
smatsj[i] =(int *)calloc(MAXELEMENTS,sizeof(int));
smatsi[i] =(int *)calloc(
     julNEQS*(julNLAGS+julNLEADS)+1,sizeof(int));
}


@}

@d main storage allocations determined at compile time
@{
julliardShocks=(double *)calloc(
     julNEQS*(SHOCKS),sizeof(double));
julliardData=(double *)calloc(
     julNEQS*(DATA),sizeof(double));
julliardFP=(double *)calloc(
     julNEQS*(julNLAGS+julNLEADS+1),sizeof(double));
AMqMatrix=(double *)
   calloc(MAXELEMENTS,sizeof(double));
AMqMatrixj=(int *)
   calloc(MAXELEMENTS,sizeof(int));
AMqMatrixi=(int *)
   calloc((julNEQS*(julNLEADS+julNLAGS)),
        sizeof(int));
@}


@d main storage deallocations
@{
cfree(failedQ);
cfree(AMqMatrix);
cfree(AMqMatrixj);
cfree(AMqMatrixi);
cfree(julliardShocks);
cfree(julliardData);
cfree(julliardPermVec);
cfree(julliardFP);
cfree(julliardPathQ);
for(i=0;i<(*pathLength)+julNLAGS+1;i++){
cfree(smats[i]);
cfree(smatsj[i]);
cfree(smatsi[i]);
cfree(fmats[i]);
cfree(fmatsj[i]);
cfree(fmatsi[i]);
}
cfree(fmats);
cfree(fmatsj);
cfree(fmatsi);
cfree(smats);
cfree(smatsj);
cfree(smatsi);
@}


chk returns 0 for success 1 if the routine may not have converged.

@d main variable declarations
@{
int * failedQ;
int chk[1]={0};
@}
@d obtain fixed point for terminal constraint
@{
printf("$Id: stochRun.w,v 1.6 2000/12/06 14:53:34 m1gsa00 Exp m1gsa00 $\n");

julModPeriodicPointGuesser(julParam,1,julliardFP);

FPnewt(numberOfEquations,lags,leads,
julMod,julModDerivative,julParam,
julliardFP,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
chk);

/*
if(chk[0])
{  printf("problems computing  FP solution\n");return(1);}
else {printf("computed FP solution\n");}
*/

printf("computed FP solution\n");


@}


@d main variable declarations
@{
/*int * hColumns;*/
@}

@d main variable declarations
@{
double * AMqMatrix;
int * AMqMatrixj;
int * AMqMatrixi;
/*double * asymptoticLinearization;*/
@|
asymptoticLinearization AMqMatrix
@}

@d main storage allocations determined at compile time
@{
*totalTime=dtime(userSystemTime);
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


@d process command line
@{
/*int  dtime(double * userSystemTime);*/

/**totalTime=dtime(userSystemTime);*/
printf("after processing command lines\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));

@}

@d carryout stochastic sims
@{
/**totalTime=dtime(userSystemTime);*/
printf("after using Q matrix\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));

@}

\appendix


\section{Non Linear Model Definition} 
\label{sec:def}

\subsection{Julliard Model}
\label{sec:julmod}


\listinginput{1}{/mq/home4/m1gsa00/aim/nonLinearModels/julliard/julliardModel.m}

\subsection{Mathematica Preparation Step}
\label{sec:prepstep}


\listinginput{1}{/mq/home4/m1gsa00/aim/nonLinearModels/julliard/prepareJulliardModel.m}

\section{Makefile}
\label{sec:makefile}

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
