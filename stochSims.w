\documentclass{article}
\newcommand{\myamp}{&}
\newcommand{\mywedge}{^}
\usepackage{html}
\usepackage{notebook}
\usepackage{amsmath}
\usepackage{latexsym}
\begin{document}
\author{Gary S. Anderson}
\title{A ``C'' Implemention for Stochastic Simulations}
\maketitle
\begin{itemize}
\item reconcile lags stoch sims
\item demonstrate  stochastic divergence nl or pl improves
\item -v exog flat
\item troll comparison
\item starter logs
\item helen report
\item scaling
\item gmres
\item term condition
\item lin aim papers
\item localize sparse lapac allocations signals
\item check aim
\item obstruct
\item tests by routine and error trapping
\item silent options
\item arg lists
\end{itemize}
A {\bf modelObj}  contains info describing 
model equations data and shocks.

\begin{itemize}
\item shocks[modelObj] produces a list of shocks
\item xData[modelObj] produces a list of data
\item eqns[modelObj] produces model equations
\item lags[modelObj] produces number of lags
\item leads[modelObj] produces number of leads
\item func[modelObj] produces a Mathematica function evaluating the model equations
\item drvFunc[modelObj] produces a Mathematica function evaluating the model equations
\item qMat[modelObj] produces a  matrix of asymptotic linear constraints
\item fp[modelObj] produces the vector of values used for computing the tail of the solution path
\end{itemize}




The stochSim function has six arguments:
\begin{description}
\item[{\bf t0Index}($t_0$)] An integer indicating the offset into the data for the
initial forecast date(Length of stochastic path  given by tfIndex-t0Index).
\item[{\bf tfIndex}($t_f$)] An integer indicating the offset into the data for the
final forecast date(Length of stochastic path  given by tfIndex-t0Index).
\item[{\bf replications}($r$)] An integer indicating the number of stochastic
paths.
\item[{\bf model}($\mathcal{M}$)]An ``object'' characterizing the model
\item[{\bf horizon}($T$)]Length of solution horizon for computing perfect 
foresight solution
\item[{\bf expType}($e$)]Expectations type(tMinusOne or t) for determining how
to use perfect foresight solution in computing $x_t$
\end{description}

The output consists of list of $r$ solutions. The solutions consist of
the values for the variables along the dynamically simulated 
stochastic path and a
list of the integers indexing the shocks that generated the path.

The code requires that $t_0-\tau>0$ and $t_f\ge t_0$, $T>0$.

Setting initial  path to data before call to generatePathX

@o stochProto.m
@{
$Path=Append[$Path,"/mq/home4/m1gsa00/aim/frbus"];

Needs["stack`"];
Needs["Statistics`DiscreteDistributions`"];
Needs["Statistics`DescriptiveStatistics`"];
Needs["Statistics`ContinuousDistributions`"];
lags[model_]:=-First[lagsLeads[model]];
leads[model_]:=Last[lagsLeads[model]];
eqns[model_]:=Length[model];
func[model_]:=modelPrep[duh[model]][[1]];
drvFunc[model_]:=modelPrep[duh[model]][[2]];
xData[model_]:=modelData[model]
shocks[model_]:=modelShocks[model]

stochSim[t0Index_Integer,tfIndex_Integer,
        replications_Integer,
		model_List,horizon_Integer,expType_Symbol]:=
With[{mlags=lags[model]},
With[{laggedDataValues=xData[model][[t0Index-Range[-mlags,-1]]]},
With[{shockSeqList=
generateDraws[t0Index,tfIndex,replications,Length[shocks[model]]],
iterFunc=(generatePathX[{model,horizon,expType,laggedDataValues},#])&},
Map[iterFunc,shockSeqList]
]]]/;(t0Index>lags[model] && tfIndex>=t0Index && 
   horizon>0 && replications>0&&((expType == tMinusOne)||(expType == t)))

@}
@o stochSims.h -d
@{
#ifndef STOCHSIMS_H
#define STOCHSIMS_H


namespace stochSims {

//#define MAXELEMENTS 20000
#define PATHLENGTH 10
#define REPLICATIONS 5000


//unsigned int intControlParameters[100];
//double doubleControlParameters[100];
//unsigned int intOutputInfo[100];
//double doubleOutputInfo[100];

#define widthIntControlInfo 70
#define widthDoubleControlInfo 50
#define widthIntOutputInfo 20
#define widthDoubleOutputInfo 10
#define useFixedPoint  0
#define useCrawlingDataPoint  1
#define useWaggingTail  2
#define useFixedDataPoint  3
#define useTailSolutionPlus  4
#define useBigEasy 1


struct controlInfo
{
unsigned int useStackQ;
unsigned int useLnsrchQ;
unsigned int numberOfDebugPairs;
unsigned int useIdentityQ;
unsigned int maxitsInput;
unsigned int maxNumberStackElements;
unsigned int maxNumberAimElements;
unsigned int numberAlphas;
unsigned int numberBetas;
unsigned int useFirstDiffQ;
unsigned int ma50PivotSearch;
unsigned int useQQ;
unsigned int terminalConstr;
unsigned int numberVarsToMonitor;
unsigned int monitoredVars[10];
/*also reserve next 9 for monotored Vars*/
unsigned int type3Q;
unsigned int debugQ;
unsigned int numberVarsToShock;
unsigned int shockedVars[9];
/*also reserve next 9 for shocked Vars*/
unsigned int tMinusOneQ;
unsigned int debugPairs[10];
/*also reserve next 10 for debug pairs*/
unsigned int homotopyXGuess[10];
unsigned int homotopyEasy;
unsigned int useBigX=0;
unsigned int useShockFileQ;
unsigned int shockFileOffset;
unsigned int dataFileOffset;
unsigned int streamingQ;
unsigned int ICshockVecLength ;
unsigned int ICnumberOfEquation;
unsigned int ICnumberOfLags;
unsigned int ICnumberOfLeads;
unsigned int ICnumberOfParameters;
unsigned int ICnumberOfDataValues;
unsigned int ICnumberOfShocks;
unsigned int ICnumberExog;
unsigned int ignoreFailQ;
double tolxInput ;
double tolfInput ;
double shrinkFactorInput ;
double expandFactorInput ;
double alaminInput ;
double alfInput ;
double ma50DropTol ;
double homotopyAlpha[10];
/*also reserve next 9 for homotopyAlpha*/
double homotopyBeta[20];
/*also reserve next 9 for homotopyBeta*/
double ma50Balance ;
double ma50DropEntry ;
double ma50DropCol ;
double shockScalar ;
};
struct outputInfo {
unsigned int failedQ;
unsigned int newtonSteps;
unsigned int fEvals;
unsigned int fDrvEvals;
unsigned int shrinkSteps;
unsigned int expandSteps;
unsigned int lnsrchStepsfailedQ;
unsigned int homotopies;
unsigned int homotopyFailures;
unsigned int currentReplication;
unsigned int currentDate;
double tolF;
double failedX;
unsigned int usePreviousHomotopyQ;
};

void addOneToFailedQ(outputInfo & inInfo){}
void subOneFromFailedQ(outputInfo & inInfo){}
void resetFailedQ(outputInfo & inInfo){}

void addOneToNewtonSteps(outputInfo & inInfo){}
void resetNewtonSteps(outputInfo & inInfo){}

void addOneToFEvals(outputInfo & inInfo){}
void resetFEvals(outputInfo & inInfo){}

void addOneToFDrvEvals(outputInfo & inInfo){}
void resetFDrvEvals(outputInfo & inInfo){}

void addOneToShrinkSteps(outputInfo & inInfo){}
void resetShrinkSteps(outputInfo & inInfo){}

void addOneToExpandSteps(outputInfo & inInfo){}
void resetExpandSteps(outputInfo & inInfo){}

void addOneToLnsrchSteps(outputInfo & inInfo){}
void resetLnsrchSteps(outputInfo & inInfo){}

void addOneToHomotopies(outputInfo & inInfo){}
void resetHomotopies(outputInfo & inInfo){}

void addOneToHomotopyFailures(outputInfo & inInfo){}
void resetHomotopyFailures(outputInfo & inInfo){}

void currentReplication(outputInfo & inInfo){}
void currentDate(outputInfo & inInfo){}

void assignRealizedTolf(outputInfo & inInfo){}
void resetRealizedTolf(outputInfo & inInfo){}

void assignRealizedTolx(outputInfo & inInfo){}
void resetRealizedTolx(outputInfo & inInfo){}

void usePreviousHomotopyQ(outputInfo & inInfo){}

/*void free();*/
/*void * calloc(unsigned num,unsigned int amt);*/
/*void pathNewt(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),void (* fdjac)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),double * params,double * shockVec,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPath,double * intercept,double * linearizationPoint,
unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double x[],
unsigned int *check,double * lastDel,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * pathNewtMa50bdJob,
unsigned int * pathNewtMa50bdIq,
double * pathNewtMa50bdFact,
unsigned int * pathNewtMa50bdIrnf,
unsigned int * pathNewtMa50bdIptrl,
unsigned int * pathNewtMa50bdIptru
);*/
long ignuin_(long low,long high);
void phrtsd_(char* phrase,long* seed1,long* seed2);
void setall_(long iseed1,long iseed2);

void generateDraws(unsigned int t0Index,unsigned int tfIndex,unsigned int replications,unsigned int shocksAvailable,
unsigned int * iarray,char * str);
//#include <stdlib.h>
//void free(void * ptr);
//void * calloc(size_t amt,size_t size);
}
@}

@o mpi.h
@{
/*nada*/
#define MPI_Datatype int
@}

@o stochSims.h -d
@{
#include <stdio.h>
namespace stochSims {

/*what goes here?*/


//#include <math.h>
#include "mpi.h"
#define TRUE 1
#define FALSE 0
#define DATA_MSG_TAG 0
#define RESULT_MSG_TAG 1
#define HALT_MSG_TAG 2

#define BUFFER_SIZE 1024

#define RESULT_ELEMENT_COUNT 3

void sendDataMessage(int,int*);
void sendHaltMessage(int);
void buildResultType(double*, int*, int*, int, int, MPI_Datatype*);
void error(int, int, int, int);
@<stochSim signature@>;
}
@}



The code makes simple draws with replacement from the empirical distribution
function(EDF) associated 
with the shocks provided by the model object.\cite{davison97}.  If there
are $n$ shocks, each vector has probability $\frac{1}{n}$ of selection
for each time period(t0Index to tfIndex) in each replication.

The following function generates $r \times (t_f-t_0+1)$ integer
indexes into the shock matrix.
@o stochProto.m
@{
generateDraws[t0Index_Integer,tfIndex_Integer,
		replications_Integer,shocksAvailable_Integer]:=
        RandomArray[
        DiscreteUniformDistribution[shocksAvailable],
        {replications,tfIndex-t0Index+1}]

@}

The following code generates a single stochastic path incorporating
a given sequence of shocks from a given initial set of $x_t$. The code
generates ``dynamic'' simulation paths, in that each solution along the
stochastic path depends on the most recently computed values in the path.

The stochastic solution path depends 
on the model, the solution horizon,  the expectations
type and the data. The routine passes these arguments to a helper
routine({\bf generateNextXTMinusOne})
that actually computes the next value of $x_t$.
The {\bf Fold} function recursively calls  the generateNextXTMinusOne function
with updated versions of the ever lengthening stochastic path.


@o stochProto.m
@{


generatePathX[{model_List,horizon_Integer,expType_Symbol,xInit_List},
				shockSeq_List]:=
		With[{shock=shocks[model],
        nxtFunc=If[expType === tMinusOne,generateNextXTMinusOne,
        generateNextXT]},
		{Fold[nxtFunc,
          {model,horizon,expType,xInit},
          shockSeq][[-1]],shockSeq}]

@}

This following routine contains all the references to the stack algorithm and to
the ``AIM'' code. This routine extends the
solution path by computing the next $x_t$ given the existing
path. It appends the new $x_t$ to the path. The routine calls {\bf aimType2 }
which applies the stack algorithm and AIM algorithm for a given {\bf stack length}
to compute the perfect foresight solution that we substitute for the
expected values in the model.
The routine calls {\bf compXEtm1} to compute the $x_t$ consistent with the given
shock and the $E(x_{t+k}|I_{t-1})$.

@o stochProto.m
@{


generateNextXTMinusOne[
{model_List,horizon_Integer,expType_Symbol,xPath_List},
	shockIndex_Integer]:=
With[{shock=shocks[model][[shockIndex]],
 numEq=Length[func[model][[2]]],
 mleads=leads[model],
 mlags=lags[model]},
With[{expectations=
(aimType2[model,horizon,
   Flatten[xPath][[
       Range[-numEq*mlags,-1]]]][[
             -1,Range[numEq*(mlags+mleads+1)]]])},
{model,horizon,expType,Append[xPath,
  compXEtm1[model,expectations,shock]]}]]

@}

The {\bf aimType2} routine applies the STACK and AIM algorithms to compute perfect foresight
solutions for the model. Note that the qMat[model] provides the precomputed AIM constraint.
Also, fp[model] provides the precomputed fixed pounsigned int for the model. This version of the
routine initializes the path beyond the lagged values to the fp[model] values.

@o generateDraws.c -d
@{

#include <random>
using namespace std;


@<define assert bump@>
#include <stdlib.h>
#include <string.h>

//#include "stochSims.h"
namespace stochSims {

void allocGenerateDraws(unsigned int t0Index,unsigned int tfIndex, unsigned int replications,unsigned int ** iarray)
{
*iarray=(unsigned int *)calloc((tfIndex-t0Index+1)*replications,sizeof(int));
}
void freeGenerateDraws(unsigned int ** iarray)
{
free(*iarray);
}


void generateDraws(unsigned int t0Index,unsigned int tfIndex,unsigned int replications,unsigned int shocksAvailable,
unsigned int * iarray,char * seedString)
{

//long seed1;
//long seed2;
static  long K1=1;
unsigned int ntot,i;
long mxint;
mxint=shocksAvailable;
std::default_random_engine generator;
std::uniform_int_distribution<int> distribution(K1,mxint);


ntot=(tfIndex-t0Index+1)*replications;
if(strcmp(seedString,"sequential")){/*need to generate random numbers*/
//phrtsd_(seedString,& seed1,& seed2);
//setall_(seed1,seed2);

    for(i=0; i<ntot; i++) {
        *(iarray+i) = distribution(generator);
    }
} else {/*generate 1...min(stochpathlength,numshocks) to fill up stochpathlength*/
for(i=0;i<ntot;i++){
iarray[i]=(i%shocksAvailable)+1;
}
}
}
}
@}

@d define assert bump
@{
#define wordybump(potentialMaxValue) \
   if(potentialMaxValue>maxElementsEncountered) \
   maxElementsEncountered=potentialMaxValue;\
printf("bump stuff(%d,%d) at line %d",potentialMaxValue,maxElementsEncountered,\
__LINE__);
#include <signal.h>

#define bump(potentialMaxValue) \
   if(potentialMaxValue>maxElementsEncountered) \
   maxElementsEncountered=potentialMaxValue;

#define pathNewtAssert(expression)  \
  if(!(expression))\
		   __pathNewtAssert (expression, __FILE__, __LINE__)

#define __pathNewtAssert(expression, file, lineno)  \
  {printf("pathNewtAssert: processid=%ld\n",getpid());\
   printf ("%s:%u: failed assertion\n", file, lineno);\
   	kill(getpid(),SIGUSR2);}


@}
@d compXEtm1 signature
@{
void compXEtm1(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),void (* fdjac)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),double * params,double * shockVec,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
double * linearizationPoint,
/*unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,*/
double easyX[],double targetX[],unsigned int * exogQ,
double x[],
unsigned int *check,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * compXMa50bdJob,
/*unsigned int * compXMa50bdIq,*/
double * compXMa50bdFact,
unsigned int * compXMa50bdIrnf,
unsigned int * compXMa50bdIptrl,
unsigned int * compXMa50bdIptru
)@}

@o stochSims.h -d
@{namespace stochSims {

@<compXEtm1 signature@>;

}
@}

@o compXEtm1.c -d
@{
@<define assert bump@>
#include  "useSparseAMA.h"
#include "stackC.h"
using namespace stackC;
//#include "stochSims.h"
using namespace stochSims;
#include <stdio.h>
@<compXEtm1 signature@>
{

FILE * debFile;
unsigned int maxElementsEncountered=0;
unsigned int maxElementsSpecified;
double * lastDel;
double * fixedPath;
double * intercept;
double * qMat;
unsigned int * qMatj;
unsigned int * qMati;
unsigned int i;
double * safex;
unsigned int aOne[1];
aOne[0]=1;
intercept =(double *)calloc(*numberOfEquations*(*leads),sizeof(double));
lastDel =(double *)calloc(*numberOfEquations*(*lags+1+ *leads),sizeof(double));
safex =(double *)calloc(*numberOfEquations*(*lags+1+ *leads),sizeof(double));
fixedPath =(double *)calloc(*numberOfEquations*(*lags+1+ *leads),sizeof(double));
qMat=(double *)calloc(*numberOfEquations* *leads,sizeof(double));
qMatj=(unsigned int *)calloc(*numberOfEquations* *leads,sizeof(int));
qMati=(unsigned int *)calloc(*numberOfEquations* *leads+1,sizeof(int));
for(i=0;i<*numberOfEquations* *leads;i++)
qMat[i]=0.0;
for(i=0;i<*numberOfEquations* *leads;i++){
qMat[i]=1.0;
qMatj[i]=i+1+(*numberOfEquations* *lags);
qMati[i]=i+1;
}
qMati[i]=i+1;


maxElementsSpecified=*maxNumberElements;
for(i=0;i<*numberOfEquations*(*lags+1+ *leads);i++){
  fixedPath[i]=targetX[i];}
for(i=0;i<*numberOfEquations*(*lags+1+ *leads);i++){
  safex[i]=x[i];}


if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + 1),
x,"x$cxtm1in");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + 1),
fixedPath,"fp$cxtm1");
fprintf(debFile,"\n(**)");
fclose(debFile);
}


pathNewt(numberOfEquations,lags, leads,aOne,
vecfunc, fdjac,params,shockVec,
fmats, fmatsj, fmatsi,
smats, smatsj, smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPath,intercept, linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
x,
check,lastDel,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
compXMa50bdJob,
/*compXMa50bdIq,*/
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
);

if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + 1),
x,"x$cxtm1out");
fprintf(debFile,"\n(**)");
fclose(debFile);
}


bump(*maxNumberElements);
*maxNumberElements=maxElementsSpecified;
if(*check !=0){
printf("compXEtm1:pathNewt nonzero shock failed. try homotopy\n");
/* easyX = targetX actual homotopy on shockVec*/
for(i=0;i<*numberOfEquations*(*lags+1+ *leads);i++){
  x[i]=easyX[i];}
homotopyPathNewt(numberOfEquations,lags, leads,aOne,
vecfunc, fdjac,params,shockVec,
fmats, fmatsj, fmatsi,
smats, smatsj, smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPath,intercept, linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
easyX,targetX,exogQ,
x,
check,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
compXMa50bdJob,
/*compXMa50bdIq,*/
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
);
bump(*maxNumberElements);
*maxNumberElements=maxElementsEncountered;
if(*check !=0){
printf("compXEtm1:homotopy failed. resetting x to computed values\n");
for(i=0;i<*numberOfEquations*(*lags+1+ *leads);i++){
  x[i]=safex[i];}
}
}


if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + 1),
x,"cxtm1");
fPrintMathDbl(debFile,*numberOfEquations,
shockVec,"shock");
fprintf(debFile,"\n(**)");
fclose(debFile);
}
free(lastDel);
free(safex);
free(intercept);
free(fixedPath);
free(qMat);
free(qMatj);
free(qMati);
}
@}

@d generateNextXTMinusOne signature
@{
void generateNextXTMinusOne(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),void (* fdjac)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),double * params,
/*unsigned int * numberExog,
double * upsilonmat,unsigned int * upsilonmatj,unsigned int * upsilonmati,void (* exdfunc)(),*/
unsigned int * shockIndex,
double * shockTable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,
double * linearizationPoint,
//unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double easyX[],double targetX[],unsigned int * exogQ,
double x[],
unsigned int *check,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * pathNewtMa50bdJob,
//unsigned int * pathNewtMa50bdIq,
double * pathNewtMa50bdFact,
unsigned int * pathNewtMa50bdIrnf,
unsigned int * pathNewtMa50bdIptrl,
unsigned int * pathNewtMa50bdIptru,
unsigned int *compXMa50bdJob,
//unsigned int * compXMa50bdIq,
double * compXMa50bdFact,
unsigned int * compXMa50bdIrnf,
unsigned int * compXMa50bdIptrl,
unsigned int * compXMa50bdIptru
)@}

@o stochSims.h -d
@{
namespace stochSims {

@<generateNextXTMinusOne signature@>;
@<generateNextXT signature@>;
}
@}


@o generateNextXTMinusOne.c -d
@{
#include <stdio.h>
#include  "useSparseAMA.h"
#include  "stackC.h"
using namespace stackC;
//#include  "stochSims.h"
using namespace stochSims;
@<define assert bump@>

@<generateNextXTMinusOne signature@>
{
FILE * debFile;
unsigned int maxElementsEncountered=0;
unsigned int maxElementsSpecified;
double * lastDel;
double * shockVec;double * safex;
unsigned int i;unsigned int ii;double *zeroShockX;double * diffFromZeroShock;
lastDel= (double *) calloc(*numberOfEquations*(*lags+*leads+*pathLength),sizeof(double));
safex= (double *) calloc(*numberOfEquations*(*lags+*leads+*pathLength),sizeof(double));
shockVec= (double *) calloc(*numberOfEquations,sizeof(double));
zeroShockX= (double *) calloc(*numberOfEquations*(*lags+*leads+1),sizeof(double));
diffFromZeroShock= (double *) calloc(*numberOfEquations,sizeof(double));
maxElementsSpecified=*maxNumberElements;
for(i=0;i<*numberOfEquations*(*leads+*pathLength);i++){
x[i+*numberOfEquations*(*lags)]=targetX[i+*numberOfEquations*(*lags)];}
for(i=0;i<*numberOfEquations*(*lags+*pathLength+ *leads);i++){
  safex[i]=x[i];}

for(i=0;i<*numberOfEquations;i++)shockVec[i]=0;

if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
x,"x$nxtxtm1in");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + 1),
fixedPoint,"fp$nxtxm1");
fprintf(debFile,"\n(**)");
fclose(debFile);
}


pathNewt(numberOfEquations,lags, leads,pathLength,
vecfunc, fdjac,params,shockVec,
fmats, fmatsj, fmatsi,
smats, smatsj, smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,intercept, linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
x,
check,lastDel,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
pathNewtMa50bdJob,
/*pathNewtMa50bdIq,*/
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru
);


if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
x,"x$nxtxtm1out");
fprintf(debFile,"\n(**)");
fclose(debFile);
}

bump(*maxNumberElements);
*maxNumberElements=maxElementsSpecified;

if(*check !=0){
printf("generateNextXTMinusOne:pathNewt zero shock failed. try homotopy\n");
for(i=0;i<*numberOfEquations*(*lags+*pathLength+ *leads);i++){
  x[i]=easyX[i];}
homotopyPathNewt(numberOfEquations,lags, leads,pathLength,
vecfunc, fdjac,params,shockVec,
fmats, fmatsj, fmatsi,
smats, smatsj, smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,intercept, linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
easyX,targetX,exogQ,
x,
check,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
pathNewtMa50bdJob,
/*pathNewtMa50bdIq,*/
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru
);
bump(*maxNumberElements);
if(*check !=0){
printf("generateNextXTminusone:homotopy failed. resetting x to computed values\n");
for(i=0;i<*numberOfEquations*(*lags+*pathLength+ *leads);i++){
  x[i]=safex[i];}
}

}
bump(*maxNumberElements);
*maxNumberElements=maxElementsSpecified;
if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
x,"genxtm1");
fprintf(debFile,"\n(**)");
fclose(debFile);
}
for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
zeroShockX[ii]=x[ii];}
if(numberVarsToShock==0){/* if 0 shock them all*/
for(ii=0;ii<*numberOfEquations;ii++){
shockVec[ii]=
 (shockScalar)*shockTable[ii+(*numberOfEquations*(*shockIndex-1))];}}
     else {for(ii=0;ii<numberVarsToShock;ii++){
         shockVec[*(&shockedVars+ii)]=
         (shockScalar)*shockTable[*(&shockedVars+ii)+
         (*numberOfEquations*(*shockIndex-1))];
printf("shocking %d with %e\n",
*(&shockedVars+ii),shockVec[*(&shockedVars+ii)]);
}}
compXEtm1(numberOfEquations,lags,leads,
vecfunc,fdjac,params,shockVec,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
zeroShockX,zeroShockX,exogQ,
x,check,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
compXMa50bdJob,
/*compXMa50bdIq,*/
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
);
if(type3Q){
for(ii=0;ii<*numberOfEquations;ii++){
diffFromZeroShock[ii]=x[ii+*numberOfEquations* *lags] - safex[ii+*numberOfEquations* *lags];
}
printf("nonzero differences from zero shock result\n");
cPrintMatrixNonZero(1,*numberOfEquations,diffFromZeroShock,1e-8);
}

for(ii=0;ii<numberVarsToMonitor;ii++){
printf("value for var number %d is %e\n",*((&monitoredVars)+ii),
x[*((&monitoredVars)+ii)+*numberOfEquations* *lags]);}

bump(*maxNumberElements);
*maxNumberElements=maxElementsEncountered;
free(shockVec);free(zeroShockX);free(diffFromZeroShock);free(lastDel);free(safex);
}



@}
@d generateNextXT signature
@{
void generateNextXT(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),void (* fdjac)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),double * params,
/*unsigned int * numberExog,*/
/*double * upsilonmat,unsigned int * upsilonmatj,unsigned int * upsilonmati,void (* exdfunc)(),*/
unsigned int * shockIndex,
double * shockTable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,
double * linearizationPoint,
//unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double easyX[],/*double targetX[],*/unsigned int * exogQ,
double x[],
unsigned int *check,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * pathNewtMa50bdJob,
//unsigned int * pathNewtMa50bdIq,
double * pathNewtMa50bdFact,
unsigned int * pathNewtMa50bdIrnf,
unsigned int * pathNewtMa50bdIptrl,
unsigned int * pathNewtMa50bdIptru/*,
unsigned int *compXMa50bdJob,
unsigned int * compXMa50bdIq,
double * compXMa50bdFact,
unsigned int * compXMa50bdIrnf,
unsigned int * compXMa50bdIptrl,
unsigned int * compXMa50bdIptru*/
)@}

@o generateNextXT.c -d
@{
#include <stdio.h>
#include  "useSparseAMA.h"
#include  "stackC.h"
using namespace stackC;
//#include "stochSims.h"
using namespace stochSims;

@<define assert bump@>




#include <stdio.h>

@<generateNextXT signature@>
{
FILE * debFile;
unsigned int maxElementsEncountered=0;
unsigned int maxElementsSpecified;
double * shockVec;double * lastDel;
unsigned int i;unsigned int ii;double *lclTargetX;double *lclEasyX;double * diffFromZeroShock;
double * safex;
shockVec= (double *) calloc(*numberOfEquations,sizeof(double));
lastDel= (double *) calloc(*numberOfEquations*(*lags+*leads+*pathLength),sizeof(double));
lclTargetX= (double *) calloc(*numberOfEquations*(*lags+*leads+*pathLength),sizeof(double));
lclEasyX= (double *) calloc(*numberOfEquations*(*lags+*leads+*pathLength),sizeof(double));
safex= (double *) calloc(*numberOfEquations*(*lags+*leads+*pathLength),sizeof(double));
for(i=0;i<*numberOfEquations*(*lags+*pathLength+ *leads);i++){
  safex[i]=x[i];}

printf("using T expectations, duh\n");
diffFromZeroShock= (double *) calloc(*numberOfEquations,sizeof(double));
maxElementsSpecified=*maxNumberElements;
switch(homotopyEasy){
case useBigEasy:
for(i=0;i<*numberOfEquations*(*lags+*leads+*pathLength);i++){
lclTargetX[i]=x[i];
lclEasyX[i]=easyX[i];
}
break;
case useBigX:
for(i=0;i<*numberOfEquations*(*lags+*leads+*pathLength);i++){
lclTargetX[i]=x[i];
lclEasyX[i]=x[i];
}}
if(numberVarsToShock==0){
for(ii=0;ii<*numberOfEquations;ii++){
shockVec[ii]=
 (shockScalar)*shockTable[ii+(*numberOfEquations*(*shockIndex-1))];}}
     else {for(ii=0;ii<numberVarsToShock;ii++){
         shockVec[*(&shockedVars+ii)]=
         (shockScalar)*shockTable[*(&shockedVars+ii)+
         (*numberOfEquations*(*shockIndex-1))];
printf("shocking %d with %e\n",
*(&shockedVars+ii),shockVec[*(&shockedVars+ii)]);
}}


if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
x,"x$nxtxin");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
fixedPoint,"fp$nxtx");
fprintf(debFile,"\n(**)");
fclose(debFile);
}


pathNewt(numberOfEquations,lags, leads,pathLength,
vecfunc, fdjac,params,shockVec,
fmats, fmatsj, fmatsi,
smats, smatsj, smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,intercept, linearizationPoint,//exogRows,exogCols,exogenizeQ,
x,
check,lastDel,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
pathNewtMa50bdJob,
//pathNewtMa50bdIq,
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru
);

if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
x,"x$nxtxout");
fprintf(debFile,"\n(**)");
fPrintMathDbl(debFile,*numberOfEquations,
shockVec,"shock");
fprintf(debFile,"\n(**)");
fclose(debFile);
}


bump(*maxNumberElements);
*maxNumberElements=maxElementsSpecified;

if(*check !=0){
printf("generateNextXT:pathNewt  shock failed. try homotopy\n");
for(i=0;i<*numberOfEquations*(*lags+*pathLength+ *leads);i++){
  x[i]=easyX[i];}
homotopyPathNewt(numberOfEquations,lags, leads,pathLength,
vecfunc, fdjac,params,shockVec,
fmats, fmatsj, fmatsi,
smats, smatsj, smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,intercept, linearizationPoint,//exogRows,exogCols,exogenizeQ,
lclEasyX,lclTargetX,exogQ,
x,
check,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
pathNewtMa50bdJob,
/*pathNewtMa50bdIq,*/
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru
);
bump(*maxNumberElements);
if(*check !=0){
printf("generateNextXt:homotopy failed. resetting x to computed values\n");
for(i=0;i<*numberOfEquations*(*lags+*pathLength+ *leads);i++){
  x[i]=safex[i];}
}

}
bump(*maxNumberElements);
*maxNumberElements=maxElementsSpecified;
if(debugQ&&(currentRequestedQ(intControlParameters,intOutputInfo))){
printf("writing  to codeGenDebFile\n");
debFile=fopen("codeGenDebFile","a");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
x,"genxtm1");
fprintf(debFile,"\n(**)");
fclose(debFile);
}
if(type3Q){
for(ii=0;ii<*numberOfEquations;ii++){
diffFromZeroShock[ii]=x[ii+*numberOfEquations* *lags] - safex[ii+*numberOfEquations* *lags];
}
printf("nonzero differences from zero shock result\n");
cPrintMatrixNonZero(1,*numberOfEquations,diffFromZeroShock,1e-8);
}

for(ii=0;ii<numberVarsToMonitor;ii++){
printf("value for var number %d is %e\n",*((&monitoredVars)+ii),
x[*((&monitoredVars)+ii)+*numberOfEquations* *lags]);}

bump(*maxNumberElements);
*maxNumberElements=maxElementsEncountered;
free(shockVec);free(lclTargetX);free(lclEasyX);free(diffFromZeroShock);free(lastDel);free(safex);
}



@}
@d generatePathX signature
@{
void generatePathX(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),void (* fdjac)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),double * params,
/*unsigned int * numberExog,
double * upsilonmat,unsigned int * upsilonmatj,unsigned int * upsilonmati,void (* exdfunc)(),*/
unsigned int * numberOfShocks,
unsigned int * shockIndices,
double * shockTable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,
double * linearizationPoint,
//unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double easyX[],double targetX[],unsigned int * exogQ,
double x[],
unsigned int *check,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * pathNewtMa50bdJob,
//unsigned int * pathNewtMa50bdIq,
double * pathNewtMa50bdFact,
unsigned int * pathNewtMa50bdIrnf,
unsigned int * pathNewtMa50bdIptrl,
unsigned int * pathNewtMa50bdIptru,
unsigned int * compXMa50bdJob,
//unsigned int * compXMa50bdIq,
double * compXMa50bdFact,
unsigned int * compXMa50bdIrnf,
unsigned int * compXMa50bdIptrl,
unsigned int * compXMa50bdIptru
)
@}


@o stochSims.h
@{
namespace stochSims {

@<generatePathX signature@>;
@<currentRequestedQ signature@>;
}
@}

@d currentRequestedQ signature
@{
unsigned int currentRequestedQ(unsigned int * intControlParameters,unsigned int * intOutputInfo)
@}

@o generatePathX.c -d
@{
#include <stdio.h>

@<define assert bump@>

#include "stackC.h"
using namespace stackC;

//#include "stochSims.h"
using namespace stochSims;
void failNextX(unsigned int * numberOfEquations,double * x)
{
unsigned int i;
for(i=0;i<*numberOfEquations;i++){x[i]=-999999.999999;}
}

@<currentRequestedQ signature@>
{
  unsigned int result=0;unsigned int ii;
for(ii=0;ii<numberOfDebugPairs;ii++){
result=result+((currentReplication==*((&debugPairs)+2*ii))&&(currentDate==*(&(debugPairs)+2*ii+1)));
}
return(result);
}
using namespace std;
@<generatePathX signature@>
{
double * lclFixedPoint;
unsigned int maxElementsEncountered=0;
unsigned int maxElementsSpecified;
unsigned int i;unsigned int j;unsigned int ii;


 lclFixedPoint=(double *)calloc(*numberOfEquations*(*lags+*leads+*pathLength),sizeof(double));
 for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
  lclFixedPoint[ii]=fixedPoint[ii];}
  maxElementsSpecified=*maxNumberElements;
 for(i=0;i<*numberOfEquations*(*lags+*pathLength+*numberOfShocks+*leads);i++){
 x[i]=targetX[i];}



for(i=0;i<*numberOfShocks;i++){
printf("for given draw, computing for date=%d\n",i);
currentDate=i;
if(check[i]!=0){
/* probably better to use failedQ and pass final path try back
so comment this out
failNextX(numberOfEquations,x+(*numberOfEquations*i));
*/
} else {


  switch (terminalConstraintSelection)
  {
  case useCrawlingDataPoint:
 for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
  lclFixedPoint[ii]=targetX[ii+(*numberOfEquations*(*pathLength-1+i))];}
  break;
  case useFixedDataPoint:
 for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
 lclFixedPoint[ii]=targetX[ii+(*numberOfEquations*(*pathLength-1))];}
  break;
  case useWaggingTail:
 for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
  lclFixedPoint[ii]=x[ii+(*numberOfEquations*(*pathLength-1+i))];}
  break;
  }


if(tMinusOneQ){
generateNextXTMinusOne(numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,
/*numberExog,
upsilonmat,upsilonmatj,upsilonmati,exdfunc,*/
shockIndices+i,
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
lclFixedPoint,intercept,linearizationPoint,//exogRows,exogCols,exogenizeQ,
easyX+(*numberOfEquations*i),targetX+(*numberOfEquations*i),exogQ,
x+(*numberOfEquations*i),check+i,intControlParameters,doubleControlParameters,
intOutputInfo,doubleOutputInfo,
pathNewtMa50bdJob,
//pathNewtMa50bdIq,
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru,
compXMa50bdJob,
//compXMa50bdIq,
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
);bump(*maxNumberElements);
} else {
generateNextXT(numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,
/*numberExog,
upsilonmat,upsilonmatj,upsilonmati,exdfunc,*/
shockIndices+i,
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
lclFixedPoint,intercept,linearizationPoint,//exogRows,exogCols,exogenizeQ,
easyX+(*numberOfEquations*i),/*targetX+(*numberOfEquations*i),*/exogQ,
x+(*numberOfEquations*i),check+i,intControlParameters,doubleControlParameters,
intOutputInfo,doubleOutputInfo,
pathNewtMa50bdJob,
//pathNewtMa50bdIq,
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru/*,
compXMa50bdJob,
compXMa50bdIq,
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
*/);bump(*maxNumberElements);
}
printf("QUICK PATCH to improve guess of next time period!!!!!!!!!!!!!!!\n");
for(j=0;j<*numberOfEquations;j++){
if(exogQ[j%*numberOfEquations]){
} else {
x[*numberOfEquations*(*lags+*leads+i)+j]= 
x[*numberOfEquations*(*lags+*leads+i-1)+j];
easyX[*numberOfEquations*(*lags+*leads+i)+j]= 
easyX[*numberOfEquations*(*lags+*leads+i-1)+j];
}
}

*maxNumberElements=maxElementsSpecified;
if(check[i]!=0)for(j=i;j<*numberOfShocks;j++){check[j]=1;}
}

}
*maxNumberElements=maxElementsEncountered;
free(lclFixedPoint);
}

@}
@o stochSims.h -d
@{
//#include <stdio.h>
namespace stochSims {

@<streamingGeneratePath signature@>;

}
@}
@d streamingGeneratePath signature
@{
void streamingGeneratePathX(FILE * streamShocksIn,FILE * streamEasyIn,FILE * streamTargetIn,
FILE * streamPathOut,
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),void (* fdjac)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),double * params,
/*unsigned int * numberExog,
double * upsilonmat,unsigned int * upsilonmatj,unsigned int * upsilonmati,void (* exdfunc)(),*/
unsigned int * numberOfShocks,
/*unsigned int * shockIndices,*/
double * shockTable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,
double * linearizationPoint,
//unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double easyX[],double targetX[],unsigned int * exogQ,
double x[],
unsigned int *check,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * pathNewtMa50bdJob,
//unsigned int * pathNewtMa50bdIq,
double * pathNewtMa50bdFact,
unsigned int * pathNewtMa50bdIrnf,
unsigned int * pathNewtMa50bdIptrl,
unsigned int * pathNewtMa50bdIptru,
unsigned int * compXMa50bdJob,
//unsigned int * compXMa50bdIq,
double * compXMa50bdFact,
unsigned int * compXMa50bdIrnf,
unsigned int * compXMa50bdIptrl,
unsigned int * compXMa50bdIptru)
@}

@o generatePathX.c -d
@{
#include <stdio.h>
//#include "stackC.h"
//using namespace stackC;
//#include "stochSims.h"
//using namespace stochSims;
@<streamingGeneratePath signature@>
{
double * lclFixedPoint;
unsigned int maxElementsEncountered=0;
unsigned int maxElementsSpecified;
unsigned int i;unsigned int j;unsigned int ii;
unsigned int aOne=1;

 lclFixedPoint=(double *)calloc(*numberOfEquations*(*lags+*leads+*pathLength),sizeof(double));
 for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
  lclFixedPoint[ii]=fixedPoint[ii];}
  maxElementsSpecified=*maxNumberElements;
 for(i=0;i<*numberOfEquations*(*lags+*pathLength+*leads);i++){
 x[i]=targetX[i];}



for(i=0;i<*numberOfShocks;i++){
getShocks(1,ICshockVecLength,*numberOfEquations,shockTable,streamShocksIn,0,1);
printf("for given draw, computing for date=%d\n",i);
currentDate=i;
if(check[0]!=0){
/* probably better to use failedQ and pass final path try back
so comment this out
failNextX(numberOfEquations,x+(*numberOfEquations*i));
*/
} else {


  switch (terminalConstraintSelection)
  {
  case useCrawlingDataPoint:
 for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
  lclFixedPoint[ii]=targetX[ii+(*numberOfEquations*(*pathLength-1))];}
  break;
  case useFixedDataPoint:
 for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
 lclFixedPoint[ii]=targetX[ii+(*numberOfEquations*(*pathLength-1))];}
  break;
  case useWaggingTail:
 for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
  lclFixedPoint[ii]=x[ii+(*numberOfEquations*(*pathLength-1))];}
  break;
  }


if(tMinusOneQ){
generateNextXTMinusOne(numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,
/*numberExog,
upsilonmat,upsilonmatj,upsilonmati,exdfunc,*/
&aOne,
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
lclFixedPoint,intercept,linearizationPoint,//exogRows,exogCols,exogenizeQ,
easyX,targetX,exogQ,
x,check,intControlParameters,doubleControlParameters,
intOutputInfo,doubleOutputInfo,
pathNewtMa50bdJob,
//pathNewtMa50bdIq,
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru,
compXMa50bdJob,
//compXMa50bdIq,
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
);bump(*maxNumberElements);
} else {
generateNextXT(numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,
/*numberExog,
upsilonmat,upsilonmatj,upsilonmati,exdfunc,*/
&aOne,
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
lclFixedPoint,intercept,linearizationPoint,
//exogRows,exogCols,exogenizeQ,
easyX,/*targetX,*/exogQ,
x,check,intControlParameters,doubleControlParameters,
intOutputInfo,doubleOutputInfo,
pathNewtMa50bdJob,
//pathNewtMa50bdIq,
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru/*(,
compXMa50bdJob,
//compXMa50bdIq,
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
*/);bump(*maxNumberElements);
}

putData(*numberOfEquations,x+(*lags)**numberOfEquations,streamPathOut);

printf("QUICK PATCH to improve guess of next time period!!!!!!!!!!!!!!!\n");
for(j=0;j<*numberOfEquations;j++){
if(exogQ[j%*numberOfEquations]){
} else {
x[*numberOfEquations*(*lags+*leads)+j]= 
x[*numberOfEquations*(*lags+*leads-1)+j];
easyX[*numberOfEquations*(*lags+*leads)+j]= 
easyX[*numberOfEquations*(*lags+*leads-1)+j];
}
}

*maxNumberElements=maxElementsSpecified;

printf("streaming shocks data and path\n");

}


for(ii=0;ii<(*lags+*leads+*pathLength-1)* *numberOfEquations;ii++){
easyX[ii]=easyX[ii+*numberOfEquations];}
getData(1,*numberOfEquations,*numberOfEquations,
easyX+(*lags+*leads+*pathLength-1)* *numberOfEquations,streamEasyIn,0,1);

for(ii=0;ii<(*lags+*leads+*pathLength-1)* *numberOfEquations;ii++){
targetX[ii]=targetX[ii+*numberOfEquations];}
getData(1,*numberOfEquations,*numberOfEquations,
targetX+(*lags+*leads+*pathLength-1)* *numberOfEquations,streamTargetIn,0,1);

for(ii=0;ii<(*lags+*leads+*pathLength-1)* *numberOfEquations;ii++){
x[ii]=x[ii+*numberOfEquations];}
for(ii=0;ii<*numberOfEquations;ii++){
x[ii+(*lags+*leads+*pathLength-1)**numberOfEquations]=
targetX[ii+(*lags+*leads+*pathLength-1)**numberOfEquations];}


}
*maxNumberElements=maxElementsEncountered;
free(lclFixedPoint);
}
@}
@d stochSim signature
@{
void stochSim(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (*vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
void (*fdjac)(double *, double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * params,
/*unsigned int * numberExog,double * upsilonmat,unsigned int * upsilonmatj,unsigned int * upsilonmati,
void (* exdfunc)(),*/
unsigned int * replications,
unsigned int * t0,unsigned int * tf,unsigned int * permVecs,
double * shockTable,/*unsigned int * shocksAvailable,*/
/*double * dataTable,unsigned int * dataAvailable,*/
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,double * linearizationPoint,
//unsigned int *exogRows,unsigned int * exogCols, unsigned int * exogenizeQ,
double easyX[],double targetX[],unsigned int * exogQ,
double x[],
unsigned int *failedQ,
controlInfo & theControlInfo, outputInfo & theOutputInfo,
unsigned int * pathNewtMa50bdJob,
//unsigned int * pathNewtMa50bdIq,
double * pathNewtMa50bdFact,
unsigned int * pathNewtMa50bdIrnf,
unsigned int * pathNewtMa50bdIptrl,
unsigned int * pathNewtMa50bdIptru,
unsigned int * compXMa50bdJob,
//unsigned int * compXMa50bdIq,
double * compXMa50bdFact,
unsigned int * compXMa50bdIrnf,
unsigned int * compXMa50bdIptrl,
unsigned int * compXMa50bdIptru
)
@}
@d freeStochSims signature
@{
void freeStochSim(unsigned int ** failedQ)
@}
@d allocStochSims signature
@{
void allocStochSim(unsigned int stochasticPathLength,unsigned int replications,unsigned int ** failedQ)
@}


@o stochSims.h -d
@{
namespace stochSims {


unsigned int * pathLength;unsigned int * replications;unsigned int * t0;unsigned int * stochasticPathLength;
//unsigned int * intControlParameters;double* doubleControlParameters;

char  * flnm;
}

@}

@o stochSims.h -d
@{
namespace stochSims {
@<processCommandLine signature@>;
}
@}

@d processCommandLine signature
@{
void processCommandLine(int argc, const char * argv[],const char ** namesArray,unsigned int modelNEQS,const char ** paramNamesArray,unsigned int numberOfParameters,double * parameters,
double * dataValues,unsigned int numberDataValues,unsigned int numShockValues,
unsigned int * pathLength,unsigned int * replications,unsigned int * t0,unsigned int * stochasticPathLength,
controlInfo & theControlInfo,char * flnm)
@}

@o stochSims.c -d

@{

#ifdef __APPLE__
#include<strings.h>
#endif
#ifdef __linux__
#include<string.h>
#endif
/* */
#include <cstdlib>
//#include "stochSims.h"

namespace stochSims {
void modData(unsigned int numberOfEquations,unsigned int numberDataValues,double * dataVals,
			 unsigned int vbl,unsigned int t0,unsigned int tf,double val1,double val2)
{
  unsigned int t;
  for(t=t0;t<=tf&&t<numberDataValues;t++){
	if(t0==tf) {
dataVals[t*numberOfEquations+vbl]=dataVals[t*numberOfEquations+vbl]+
  (val2+val1)/2;} else {
dataVals[t*numberOfEquations+vbl]=dataVals[t*numberOfEquations+vbl]+
  (t-t0)*val2/(tf-t0) + (tf-t)*val1/(tf-t0);}
  }
}
void modDataAbs(unsigned int numberOfEquations,unsigned int numberDataValues,double * dataVals,
			 unsigned int vbl,unsigned int t0,unsigned int tf,double val1,double val2)
{
  unsigned int t;
  for(t=t0;t<=tf&&t<numberDataValues;t++){
	if(t0==tf) {
dataVals[t*numberOfEquations+vbl]= (val2+val1)/2;} else {
	dataVals[t*numberOfEquations+vbl]=(t-t0)*val2/(tf-t0) + (tf-t)*val1/(tf-t0);}
  }
}/* */


@<processCommandLine signature@>
{
flnm=(char *)calloc(500,sizeof(char));

  float aFloat;unsigned int i;unsigned int anInt;
 unsigned int pl;unsigned int t1;unsigned int t2; double val1; double val2; unsigned int vbl;

useLnsrchQ=FALSE;
useStackQ=FALSE;
useIdentityQ=FALSE;
useFirstDiffQ=FALSE;
 alaminInput=1e-10;
 alfInput=1e-4;
 expandFactorInput=1.3;
 shrinkFactorInput=0.7;
 tolfInput=1e-10;
  tolxInput=1e-10;
 maxitsInput=100;
*replications=1;
*stochasticPathLength=1;
*pathLength=1;
*t0=0;
numberAlphas=1;
homotopyAlpha[0]=1.0;
numberBetas=1;
homotopyBeta[0]=1.0;
while(argc>1&&argv[1][0] == '-')
{
printf("processing command line args\n");
 switch(argv[1][1]){
   case 'p':
i=0;
while((strcmp(argv[2],paramNamesArray[i])) && (i <numberOfParameters))i++;
if(i==modelNEQS){
printf("i don't know the parameter %s: ignoring this parmeter value pair\n",argv[2]);} else {vbl = i;}
         sscanf(argv[3],"%f",&aFloat);val1=(double)aFloat;
         printf("got %d for param %s  and %f for value\n",vbl,paramNamesArray[i],
	val1);
		 parameters[i]=val1;
         argc--;argv++;
         argc--;argv++;
         break;
   case 'v':
i=0;
while((strcmp(argv[2],namesArray[i])) && (i <modelNEQS))i++;
if(i==modelNEQS){
printf("i don't know the variable %s: ignoring this variable value pair\n",argv[2]);} else {vbl = i;}
         t1=(unsigned int)atoi(argv[3]);
         t2=(unsigned int)atoi(argv[4]);
         sscanf(argv[5],"%f",&aFloat);val1=(double)aFloat;
         sscanf(argv[6],"%f",&aFloat);val2=(double)aFloat;
         printf("got %d for vbl %s and (%d,%d) (%f,%f)\n",vbl,namesArray[i],
	t1,t2,val1,val2);
		 modData(modelNEQS,numberDataValues,dataValues,vbl,t1,t2,val1,val2);
         argc--;argv++;
         argc--;argv++;
         argc--;argv++;
         argc--;argv++;
         argc--;argv++;
         break;
   case 'V':
i=0;
while((strcmp(argv[2],namesArray[i])) && (i <modelNEQS))i++;
if(i==modelNEQS){
printf("i don't know the variable %s: ignoring this variable value pair\n",argv[2]);} else {vbl = i;}
         t1=(unsigned int)atoi(argv[3]);
         t2=(unsigned int)atoi(argv[4]);
         sscanf(argv[5],"%f",&aFloat);val1=(double)aFloat;
         sscanf(argv[6],"%f",&aFloat);val2=(double)aFloat;
         printf("got %d for vbl %s and (%d,%d) (%f,%f)\n",vbl,namesArray[i],
	t1,t2,val1,val2);
		 modDataAbs(modelNEQS,numberDataValues,dataValues,vbl,t1,t2,val1,val2);
         argc--;argv++;
         argc--;argv++;
         argc--;argv++;
         argc--;argv++;
         argc--;argv++;
         break;
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
   case 'r':
         pl=atoi(argv[2]);
     printf("got %d for replications\n",pl);
     if(pl>REPLICATIONS||pl<1)
         {
       *replications=REPLICATIONS;
       printf("setting repetitions to maximum=%d\n",REPLICATIONS);
       } else { *replications = pl;}
         argc--;argv++;
         break;
   case 'a':
         pl=atoi(argv[2]);
     printf("got %d for t0\n",pl);
     if(pl>numShockValues)
         {
       *t0=numShockValues;
       printf("initial t0 to maximum=%d\n",numShockValues);
       } else { *t0 = pl;}
         argc--;argv++;
         break;
   case 's':
         *stochasticPathLength=atoi(argv[2]);
     printf("got %d for stochasticPathLength\n",*stochasticPathLength);
     if(*stochasticPathLength<1)
         {
       *stochasticPathLength=1;
       printf("setting tf to 1\n");
       }
         argc--;argv++;
         break;
   case 'N':
         (maxitsInput)=atoi(argv[2]);
     printf("got %d for maxits\n",(maxitsInput));
     if((maxitsInput)<1)
         {
       (maxitsInput)=1;
       printf("setting maxits to 1\n");
       }
         argc--;argv++;
         break;
   case 'F':
         sscanf(argv[2],"%f",&aFloat);
		 tolfInput=(double)aFloat;
     printf("got %e for tolfInput\n",tolfInput);
         argc--;argv++;
         break;
   case 'X':
         sscanf(argv[2],"%f",&aFloat);
		 tolxInput=(double)aFloat;
     printf("got %e for tolxInput\n",tolxInput);
         argc--;argv++;
         break;
   case 'H':
         sscanf(argv[2],"%u",&anInt);
		 numberAlphas=anInt;
		 if(anInt>10){numberAlphas=10; printf("only using first 10\n");}
         argc--;argv++;
		 for(i=0;i<anInt;i++){
         sscanf(argv[2],"%f",&aFloat);
		 if(i<10){homotopyAlpha[i]=(double)aFloat;}
     printf("got %e for homotopyAlpha[%d]\n",aFloat,i);
	 argc--;argv++;}
         break;
   case 'D':
         sscanf(argv[2],"%u",&anInt);
		 numberBetas=anInt;
		 if(anInt>10){numberBetas=10; printf("only using first 10\n");}
         argc--;argv++;
		 for(i=0;i<anInt;i++){
         sscanf(argv[2],"%f",&aFloat);
		 if(i<10){homotopyBeta[i]=(double)aFloat;}
     printf("got %e for homotopyBeta[%d]\n",aFloat,i);
	 argc--;argv++;}
         break;
   case 'K':
         sscanf(argv[2],"%f",&aFloat);
		 shrinkFactorInput=(double)aFloat;
     printf("got %e for shrinkFactorInput\n",shrinkFactorInput);
         argc--;argv++;
         break;
   case 'E':
         sscanf(argv[2],"%f",&aFloat);
		 expandFactorInput=(double)aFloat;
     printf("got %e for expandFactorInput\n",expandFactorInput);
         argc--;argv++;
         break;
   case 'L':
         useLnsrchQ=TRUE;
     printf("got flag for using lnsrch algorithm\n");
         break;
   case 'S':
         useStackQ=TRUE;
     printf("got flag for using stack algorithm\n");
         break;
   case 'I':
         useIdentityQ=TRUE;
         useFirstDiffQ=FALSE;
     printf("got flag for using identity matrix\n");
         break;
   case 'J':
         useFirstDiffQ=TRUE;
         useIdentityQ=FALSE;
     printf("got flag for using first difference matrix\n");
         break;
   case 'h':
     printf("\n-l <stack pathlength>\n"); 
     printf("-s <stochastic pathlength>\n");
     printf("-r <number of replications>\n");
     printf("-f <output filename>\n");
     printf("-L  use lnsearch\n");
     printf("-I  use identity matrix terminal condition\n");
     printf("-a <offset into datamatrix>\n");
     printf("-X <x convergence tolerance>\n");
     printf("-F <f(x) convergence tolerance>\n");
     printf("-N <maximum number of newton steps>\n");
     printf("-E <expansion factor>\n");
     printf("-K <shrinkage factor>\n");
     printf("-v <variableName> <dataPt 0 > <dataPtf> <incrementVal0> <incrementValf>\n");
     printf("-V <variableName> <dataPt 0 > <dataPtf> <val0> <valf>\n");
     printf("-p <parameterName> <valf>\n");
        break;
   case 'f':
         strcpy(flnm,argv[2]);
         printf("got %s for filename \n",flnm);
         argc--;argv++;
         break;
 default:
   printf("%s: unknown arg %s-- not processing any more args\n",
      argv[0],argv[1]);
 }
argc--;argv++;
}


printf("values for run:(pathLength=%d,replications=%d,t0=%d,stochasticPathLength=%d)\n",*pathLength,*replications,*t0,*stochasticPathLength);


}

/*
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

*/
}


@}


@o stochSims.h -d
@{
namespace stochSims {

@<freeStochSims signature@>;
@<allocStochSims signature@>;
}
#endif
@}
@o stochSims.c -d
@{
@<define assert bump@>
#include <cstdlib>
#include <cstdio>

#include <stdio.h>
#include "stackC.h"
using namespace stackC;


//#include "stochSims.h"

namespace stochSims {
@<allocStochSims signature@>
{
*failedQ=(unsigned int *)calloc(stochasticPathLength*replications,sizeof(int));
}
@<freeStochSims signature@>
{
free(*failedQ);
}


void streamingStochSim(
FILE * streamShocksIn,FILE * streamEasyIn,FILE * streamTargetIn,
FILE * streamPathOut,
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),void (* fdjac)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
double * params,
/*unsigned int * numberExog,double * upsilonmat,unsigned int * upsilonmatj,unsigned int * upsilonmati,
void (* exdfunc)(),*/
unsigned int * replications,
unsigned int * t0,unsigned int * tf,/*unsigned int * permVecs,*/
double * shockTable,/*unsigned int * shocksAvailable,*/
/*double * dataTable,unsigned int * dataAvailable,*/
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,double * linearizationPoint,
//unsigned int *exogRows,unsigned int * exogCols, unsigned int * exogenizeQ,
double easyX[],double targetX[],unsigned int * exogQ,
double x[],
unsigned int *failedQ,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * pathNewtMa50bdJob,
//unsigned int * pathNewtMa50bdIq,
double * pathNewtMa50bdFact,
unsigned int * pathNewtMa50bdIrnf,
unsigned int * pathNewtMa50bdIptrl,
unsigned int * pathNewtMa50bdIptru,
unsigned int * compXMa50bdJob,
//unsigned int * compXMa50bdIq,
double * compXMa50bdFact,
unsigned int * compXMa50bdIrnf,
unsigned int * compXMa50bdIptrl,
unsigned int * compXMa50bdIptru
)
{
//unsigned int check[1]={0};
unsigned int * stochasticPathLength;
unsigned int maxElementsEncountered=0;
unsigned int maxElementsSpecified;
unsigned int i/*,j*/;
FILE * debFile;
printf("streaming shocks data and path\n");
getData(*lags+*leads+*pathLength,*numberOfEquations,*numberOfEquations,
easyX,streamEasyIn,0,*lags+*leads+*pathLength);
getData(*lags+*leads+*pathLength,*numberOfEquations,*numberOfEquations,
targetX,streamTargetIn,0,*lags+*leads+*pathLength);
for(i=0;i<(*lags+*leads+*pathLength)* *numberOfEquations;i++){
x[i]=targetX[i];}
resetHomotopies;
if(debugQ){
debFile=fopen("codeGenDebFile","a");
fprintf(debFile,"begin stochSim\n");
fclose(debFile);
}


stochasticPathLength=(unsigned int *)calloc(1,sizeof(int));
*stochasticPathLength=*tf-*t0+1;
/*for(i=0;i<*stochasticPathLength;i++)failedQ[i]=0;*/
maxElementsSpecified=*maxNumberElements;
for(i=0;i<*replications;i++){
printf("computing for draw=%d\n",i);}
currentReplication=i;

streamingGeneratePathX(streamShocksIn,streamEasyIn,streamTargetIn,
streamPathOut,
numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,
/*numberExog,
upsilonmat,upsilonmatj,upsilonmati,exdfunc,*/
stochasticPathLength,/*permVecs+i**stochasticPathLength,*/
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,intercept,linearizationPoint,
//exogRows,exogCols,exogenizeQ,
easyX,targetX,exogQ,
x,failedQ,
intControlParameters,doubleControlParameters,
intOutputInfo,
doubleOutputInfo,
pathNewtMa50bdJob,
//pathNewtMa50bdIq,
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru,
compXMa50bdJob,
//compXMa50bdIq,
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
);
bump(*maxNumberElements);
*maxNumberElements=maxElementsSpecified;

*maxNumberElements=maxElementsEncountered;
if(debugQ){
debFile=fopen("codeGenDebFile","a");
fprintf(debFile,"end stochSim\n");
fclose(debFile);
}

free(stochasticPathLength);
}


@<stochSim signature@>
{
//unsigned int check[1]={0};
unsigned int * stochasticPathLength;
unsigned int maxElementsEncountered=0;
unsigned int maxElementsSpecified;
unsigned int i/*,j*/;
FILE * debFile;

resetHomotopies;
if(debugQ){
debFile=fopen("codeGenDebFile","a");
fprintf(debFile,"begin stochSim\n");
fclose(debFile);
}


stochasticPathLength=(unsigned int *)calloc(1,sizeof(int));
*stochasticPathLength=*tf-*t0+1;
for(i=0;i<*stochasticPathLength;i++)failedQ[i]=0;
maxElementsSpecified=*maxNumberElements;
for(i=0;i<*replications;i++){
printf("computing for draw=%d\n",i);
currentReplication=i;

generatePathX(numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,
/* numberExog,
upsilonmat,upsilonmatj,upsilonmati,exdfunc,*/
stochasticPathLength,permVecs+i**stochasticPathLength,
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,intercept,linearizationPoint,
//exogRows,exogCols,exogenizeQ,
easyX,targetX,exogQ,
x+(*numberOfEquations*(*lags+*stochasticPathLength)*i),failedQ+(i* *stochasticPathLength),
intControlParameters,doubleControlParameters,
intOutputInfo+(i* widthIntOutputInfo),
doubleOutputInfo+(i* widthDoubleOutputInfo),
pathNewtMa50bdJob,
//pathNewtMa50bdIq,
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru,
compXMa50bdJob,
//compXMa50bdIq,
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
);
bump(*maxNumberElements);
*maxNumberElements=maxElementsSpecified;
};
*maxNumberElements=maxElementsEncountered;
if(debugQ){
debFile=fopen("codeGenDebFile","a");
fprintf(debFile,"end stochSim\n");
fclose(debFile);
}

free(stochasticPathLength);
}
}
@}

@o distStochSims.c -d
@{
@<define assert bump@>
//#include "stochSims.h"
//#include "distStochSims.h"


void distStochSim(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),void (* fdjac)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
double * params,
unsigned int * numberExog,
double * upsilonmat,unsigned int * upsilonmatj,unsigned int * upsilonmati,void (* exdfunc)(),
unsigned int * replications,
unsigned int * t0,unsigned int * tf,unsigned int * permVecs,
double * shockTable,unsigned int * shocksAvailable,
double * dataTable,unsigned int * dataAvailable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,double * linearizationPoint,unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double easyX[],double targetX[],unsigned int * exogQ,
double x[],
unsigned int *failedQ,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * pathNewtMa50bdJob,
unsigned int * pathNewtMa50bdIq,
double * pathNewtMa50bdFact,
unsigned int * pathNewtMa50bdIrnf,
unsigned int * pathNewtMa50bdIptrl,
unsigned int * pathNewtMa50bdIptru,
unsigned int * compXMa50bdJob,
unsigned int * compXMa50bdIq,
double * compXMa50bdFact,
unsigned int * compXMa50bdIrnf,
unsigned int * compXMa50bdIptrl,
unsigned int * compXMa50bdIptru,
unsigned int draw,
FILE * outFile
)
{
unsigned int check[1]={0};
unsigned int * stochasticPathLength;
unsigned int j;
stochasticPathLength=(unsigned int *)calloc(1,sizeof(int));
*stochasticPathLength=*tf-*t0+1;

printf("computing for draw=%d\n",draw);


    /***** following is new, 2/26/01, ray *****/
    /* initialize path to zeroes, then fill in relevant values from dataTable[]. */
    for (j=0; j < *replications * (*numberOfEquations) *
	   (*lags + *leads + *pathLength + *stochasticPathLength); j++)
      x[j] = 0.0;
    /***** end of new code *****/



    /*********************************************
    fprintf(outFile, "x vector before running generatePathX in replication %d\n", draw);
    for (j=0; j < (*replications * (*numberOfEquations) *
		   (stochasticPathLength + (*lags))); j++)
      fprintf(outFile, "start value %d   %f\n", j, x[j]);
    ********************************************/

    /*********************************************
    printInputs(numberOfEquations,lags,leads,pathLength,
		vecfunc,fdjac,params,numberOfShocks,permVecs+i**numberOfShocks,
		shockTable,
		fmats,fmatsj,fmatsi,
		smats,smatsj,smatsi,
		maxNumberElements,qMat,qMatj,qMati,
		fixedPoint,
		x+(*numberOfEquations*(*lags+*numberOfShocks)*i),
		failedQ+(i* *numberOfShocks), outFile);
    ********************************************/


generatePathX(numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,
numberExog,
upsilonmat,upsilonmatj,upsilonmati,exdfunc,
stochasticPathLength,permVecs+draw**stochasticPathLength,
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,intercept,linearizationPoint,exogRows,exogCols,exogenizeQ,
easyX,targetX,exogQ,
x+(*numberOfEquations*(*lags+*stochasticPathLength)*draw),failedQ+(draw* *stochasticPathLength),
intControlParameters,doubleControlParameters,
intOutputInfo+(draw* widthIntOutputInfo),
doubleOutputInfo+(draw* widthDoubleOutputInfo),
pathNewtMa50bdJob,
pathNewtMa50bdIq,
pathNewtMa50bdFact,
pathNewtMa50bdIrnf,
pathNewtMa50bdIptrl,
pathNewtMa50bdIptru,
compXMa50bdJob,
compXMa50bdIq,
compXMa50bdFact,
compXMa50bdIrnf,
compXMa50bdIptrl,
compXMa50bdIptru
);



    /*************************************
 fprintf(outFile, "x vector after running generatePathX in replication %d\n", draw);
 for (j=0; j < (*replications * (*numberOfEquations) *
		(stochasticPathLength + (*lags))); j++)
		***********************/
   /***  for (j=0; j < *replications * (*numberOfEquations) *
	(*lags + *leads + *pathLength + stochasticPathLength); j++) ***/
    /******************
   fprintf(outFile, "%d   %f\n", j, x[j]);
   ********************/

free(stochasticPathLength);
}
@}


@o stochProto.m
@{
aimType2[model_List,horizon_,init_]:=fixedPathList[
nxtGuess[lags[model],
 func[model],
 drvFunc[model],
 qMat[model],
 fp[model],#]&,
 Join[init,Flatten[Table[fp[model][[Range[eqns[model]]]],{(horizon+leads[model])}]]],
 SameTest->(Max[Abs[#1-#2]]<10^(-15)&)]


altAimType2[model_,horizon_,init_]:=
With[{theArgs={{0., 0., 0., 0., 0., 0.0641495, -0.172672, -0.222054, -0.190656, 0.140763, 
  0., -0.0488314, 0.00496178, -0.102767, 0.043727}, 
 {0.0641495, -0.172672, -0.222054, -0.190656, 0.140763, 0., -0.0488314, 
  0.00496178, -0.102767, 0.043727, 0., -0.0158041, -0.00253689, -0.033196, 
  0.0142162}}},
fixedPathList[
nxtGuess[lags[model],
 func[model],
 drvFunc[model],
 qMat[model],
 fp[model],#]&,
 Join[init,Flatten[Table[theArgs[[2]][[5+Range[eqns[model]]]],{(horizon+leads[model])}]]],
 SameTest->(Max[Abs[#1-#2]]<10^(-15)&)]]


@}

The following routine computes the $x_t$ associated with a given set of expectations 
and a given shock. It uses the STACK algorithm code to compute the solution.



@o stochProto.m
@{

compXEtm1[model_List,expectations_List,shock_List]:=
With[{deFunc=func[model]},
With[{mlags=lags[model],mleads=leads[model],
  numEq=Length[func[model][[2]]],
  newFunc=
    Apply[Function,{deFunc[[1]],deFunc[[2]]-shock}]},
With[{init=expectations[[Range[numEq*mlags]]]},
 fixedPath[
  nxtGuess[lags[model],
  newFunc,
  drvFunc[model],
  BlockMatrix[{{ZeroMatrix[numEq*mleads,numEq*(mlags)],IdentityMatrix[numEq*mleads]}}],
  expectations,#]&,
  expectations,
  SameTest->(Max[Abs[#1-#2]]<10^(-15)&)][[numEq*mlags+Range[numEq]]]]]]


@}

This routine computes the $x_t$ consistent with time ``t'' expectations and a given
shock associated with the shockIndex.

@o stochProto.m
@{

generateNextXT[
{model_List,horizon_Integer,expType_Symbol,xPath_List},
	shockIndex_Integer]:=
With[{shock=shocks[model][[shockIndex]],
 numEq=Length[func[model][[2]]],
 mleads=leads[model],
 mlags=lags[model]},
With[{expectations=
(aimType2[model,horizon,
   Flatten[xPath][[
       Range[-numEq*mlags,-1]]]][[
             -1,Range[numEq*(mlags+mleads+1)]]])},
{model,horizon,expType,Append[xPath,aimType2Terror[shock,model,horizon,
   Flatten[xPath][[
       Range[-numEq*mlags,-1]]]][[
             -1,numEq*mlags+Range[numEq]]]]}]]



@}
This routine computes the update for the path using the stack routines for computing
time ``t'' period expectations.

@o stochProto.m
@{

nxtGuessTerror[shock_List,nlag_Integer,
theFunc_Function,theDrvFunc_Function,
termConstr_List,fp_List,guess_List]:=
With[{neq=Length[theDrvFunc[[2]]]},
With[{nlead=(Length[theDrvFunc[[2,1]]]/neq)-nlag-1},
Module[{nxlstC,nxlstD,lstC,lstD},
With[{appDim=(nlag+nlead+1)*neq,
theZeroMatsC=Table[ZeroMatrix[neq,neq],{nlag}],
theZeroMatsD=Table[ZeroMatrix[neq,1],{nlag}]
},
With[{theArgs=Partition[guess,appDim,neq]},
With[{theRes=Map[Apply[theFunc,#]& , theArgs]},
With[{newFunc=
    Apply[Function,{theFunc[[1]],theFunc[[2]]-shock}]},
With[{prime=nxtCDmats[{nlag,
Apply[theDrvFunc,theArgs[[1]]],
Apply[newFunc,theArgs[[1]]],
theZeroMatsC,theZeroMatsD}]},(*Print["prime=",prime,"applic=",Apply[newFunc,theArgs[[1]]],"theArgs[[1]]=",theArgs[[1]],"shock=",shock,"theArgs=",(*InputForm[Rationalize[#,1/10000000000000000]]& MapThingy*)theArgs,"theRes",theRes];*)
{nxlstC,nxlstD}=Fold[Function[{x,y},nxtCDmats[{nlag,
Apply[theDrvFunc,y],(*Print["applying theFunc=",Apply[theFunc,y],y];*)
Apply[theFunc,y],x[[1]],x[[2]]}]],
prime(*{theZeroMatsC,theZeroMatsD}*),Drop[theArgs,1]];
{lstC,lstD}=nxtCDmats[{nlag,
BlockMatrix[{{termConstr,ZeroMatrix[neq]}}],
termConstr . (Transpose[{guess[[Range[-neq*(nlag+nlead),-1]]]}]- 
Transpose[{fp[[Range[-neq*(nlag+nlead),-1]]]}]),nxlstC,nxlstD}];
bs=backSub[lstC,lstD];
guess-Flatten[bs]]]]]]]]]

@}
This routine applies nxtGuessTerror to 
compute the time t 
@o stochProto.m
@{

aimType2Terror[shock_,model_, horizon_, init_] := fixedPathList[
    nxtGuessTerror[shock,lags[model], func[model], drvFunc[model], qMat[model], 
      fp[model], #1] & , Join[init, Flatten[Table[fp[model][[Range[eqns[model]]]], 
       {horizon+1 + leads[model]}]]], SameTest -> 
     (Max[Abs[#1 - #2]] < 10^(-15) & )]


@}



This debugging routine extends the computation horizon to obtain a convergent
set of perfect foresight values.

@o stochProto.m
@{

aimType3[model_,initHoriz_,init_]:=
With[{numEq=Length[func[model][[2]]]},
fixedPathList[
{#[[1]]+1,aimType2[model,#[[1]]+1,init][[-1,numEq+Range[numEq]]]}&,{0,Table[0,{numEq}]},
SameTest->(Max[Abs[#1[[2]]-#2[[2]]]]<10^(-15)&)]];

altAimType3[model_,initHoriz_,init_]:=
With[{numEq=Length[func[model][[2]]]},
fixedPathList[
{#[[1]]+1,altAimType2[model,#[[1]]+1,init][[-1,numEq+Range[numEq]]]}&,{0,Table[0,{numEq}]},
SameTest->(Max[Abs[#1[[2]]-#2[[2]]]]<10^(-15)&)]];

@}

This debugging routine computes the $t-1$ error associated with a given set of expectations.
@o stochProto.m
@{

tMinusOneErrorFunc[model_,expectations_List]:=
With[{mlags=lags[model],mleads=leads[model],
  numEq=Length[func[model][[2]]]},
With[{xAtTZero=Table[Unique[],{numEq}]},
With[{argList=Join[expectations[[Range[numEq*mlags]]],xAtTZero,
  expectations[[numEq*(mlags+1)+Range[numEq*mleads]]]]},
Apply[Function,{xAtTZero,Apply[func[model],argList]}]]]]

@}

The example output uses the  following model definitions:

@o stochProtoTest.m
@{
$Path=Append[$Path,"/mq/home/m1gsa00/bestMath/"];
<<amsMatrices.m
(*testModel definition*)
(* pg 14 example juillard, laxton, mcadam, pioro *)
paperModel=aMod[{
ey[t],
pdot[t]  - (0.414 pdot[t+1] + (1 - 0.414)* pdot[t-1] + 
  0.196 (g^2/(g-y[t]) - g) + 0.276(g^2/(g - y[t-1]) - g)),
rr[t]- (rs[t] - 0.414 pdot[t+1] - (1-0.414)pdot[t-1]),
rs[t] - (3 pdot[t]+ y[t]),
y[t] - (0.304 y[t-1] - 0.98 rr[t] - 0.315 rr[t-1] - ey[t-1])
}/.g->0.049]

testModel=.;
lags[testModel]=1;
leads[testModel]=1;
eqns[testModel]=5;
(*
shocks[testModel]=
   RandomArray[NormalDistribution[0,0.1],{30,5}];
*)
<<rayMemoShocks0531.m
shocks[testModel]=shk;

xData[testModel]=Transpose[{Range[-0.3,0.4,0.01],Range[-0.6,0.1,0.01],
    Range[-0.3,0.4,0.01],Range[-0.6,0.1,0.01],Range[-0.3,0.4,0.01]}];
@}


The STACK algorithm routines require Mathematica functions characterizing
the model and it's derivatives.


@o stochProtoTest.m
@{


paperFunc = Apply[Function, 
With[{vbls=eqnVars[paperModel],ll=lagLead[paperModel]},
With[{argVal=Flatten[Table[Map[#[t+i]&,vbls],{i,-ll[[1]],ll[[2]]}]]},
With[{theSubs=Thread[argVal->Table[Unique[],{15}]]},
{argVal/.theSubs,paperModel[[1]]/.theSubs}]]]]

func[testModel]=paperFunc;


With[{allVars=eqnVars[paperModel],ll=lagLead[paperModel]},
drvsModel=((Map[Function[x,Map[D[x,#]&,
  Flatten[Table[Through[allVars[t+i]],{i,-ll[[1]],ll[[2]]}]]]],(paperModel[[1]])]))]


paperDrvFunc = Apply[Function, 
With[{vbls=eqnVars[paperModel],ll=lagLead[paperModel]},
With[{argVal=Flatten[Table[Map[#[t+i]&,vbls],{i,-ll[[1]],ll[[2]]}]]},
With[{theSubs=Thread[argVal->Table[Unique[],{15}]]},
{argVal/.theSubs,drvsModel/.theSubs}]]]]

drvFunc[testModel]=paperDrvFunc;

@}

The following code computes applies AIM to compute the asymptotic constraints for
the model.

@o stochProtoTest.m
@{



{pvbls,pssFuncs}=steadyState[paperModel]
pres=ssLinearizeAndMakeHMat[paperModel,pssFuncs[[1]]]
{{af,ar},{ab,rar}}=symbolicBiDirectionalAR[pres[]]
transMat=symbolicTransitionMatrix[ar];

{ubigvals,ubigEvs}=Eigensystem[Transpose[transMat]];
qmat=Join[af,ubigEvs[[{1}]]];
qMat[testModel]=qmat;
fp[testModel]=Table[0,{5*3}];

paperQ=BlockMatrix[{{qmat,ZeroMatrix[5]}}];


@}

Here are some example applications of the functions.
Also preparing splice for documentation.

@o stochProtoTest.m
@{
shk=shocks[testModel];
(*
matrixRowToAmsForm[shk[[1]],1,3]
*)
Splice["stochRes.mtex"]

(*
Timing[tryStoch=stochSim[2,10,30,testModel,1,tMinusOne]];
flt= Map[Flatten[#[[1]]]& , tryStoch];
generateNextXT[{testModel,1,try,Table[Range[5],{3}]},1];
generateNextXTMinusOne[{testModel,1,try,Table[Range[5],{3}]},1];

paperFunc @@ Join[{0,0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}]
type2=aimType2Terror[testModel,1,plugTM102]
huh=generateNextXT[{testModel,1,t,Table[{0,0,0,0,0},{1}]},1];
huh2=aimType2[testModel,1,huh[[4,2]]];
huh2a=altAimType2[testModel,1,huh[[4,2]]];
(*tail huh2 uses linear constraint*)
bmat=-Inverse[SubMatrix[qmat,{1,6},{5,5}]] .SubMatrix[qmat,{1,1},{5,5}];
Max[Abs[Transpose[{huh2[[-1,Range[5]+10]]}] - bmat . Transpose[{huh2[[-1,Range[5]+5]]}]]]
(*satisfies non linear constraint*)
Max[Abs[paperFunc @@ huh2[[-1]]]]

stochSim[2,2,5,testModel,5,tMinusOne];
huh0=generateNextXT[{testModel,0,t,Table[{0,0,0,0,0},{1}]},1];
huh02=aimType2[testModel,0,huh0[[4,2]]];
bmat . Transpose[{huh02[[-1,Range[5]+5]]}]
timeTEx02=generatePathX[{testModel,1,t,Table[{0,0,0,0,0},{1}]},{1,15,6}];
plugT00=timeTEx02[[1,2,Range[5]]];
pathOne=altAimType2[testModel,1,plugT00][[-1]];
pathTwo=aimType2[testModel,1,plugT00][[-1]];

(*different*)
pathOne00=altAimType2[testModel,1,0.5*plugT00][[-1]];
pathTwo00=aimType2[testModel,1,0.5*plugT00][[-1]];


(*same*)
pathOne01=altAimType2[testModel,1,0.1*plugT00][[-1]];
pathTwo01=aimType2[testModel,1,0.1*plugT00][[-1]];


pathOne02=altAimType2[testModel,1,0.25*plugT00][[-1]];
pathTwo02=aimType2[testModel,1,0.25*plugT00][[-1]];


*)

@}
\appendix


\include{stochRes}

\bibliographystyle{plain}
\bibliography{files}
\nocite{hollinger96}
\end{document}
