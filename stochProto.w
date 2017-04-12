\documentclass{article}
\newcommand{\myamp}{&}
\newcommand{\mywedge}{^}
%\usepackage{notebook}
\usepackage{amsmath}
\usepackage{latexsym}
\begin{document}
\author{Gary S. Anderson}
\title{Pseudo-code for Implementing Stochastic Simulations}
\maketitle

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
With[{laggedDataValues=xData[model][[t0Index+Range[-mlags,-1]]]},
With[{shockSeqList=
generateDraws[t0Index,tfIndex,replications,Length[shocks[model]]],
iterFunc=(generatePathX[{model,horizon,expType,laggedDataValues},#])&},
Map[iterFunc,shockSeqList]
]]]/;(t0Index>lags[model] && tfIndex>=t0Index && 
   horizon>0 && replications>0&&((expType == tMinusOne)||(expType == t)))

@}




@o stochProto.c -d
@{
#include "stochProto.h"
void free();
void pathNewt(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,double * shockVec,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,
double x[],
unsigned int *check);
long ignuin_(long *low,long *high);


void generateDraws(unsigned int t0Index,unsigned int tfIndex,unsigned int replications,unsigned int shocksAvailable,
unsigned int * iarray)
{
static  long K1=1;
unsigned int ntot,i;
long mxint;
ntot=(tfIndex-t0Index+1)*replications;
mxint=shocksAvailable;
    for(i=0; i<ntot; i++) {
        *(iarray+i) = (unsigned int )ignuin_(&K1,&mxint);
    }
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
Also, fp[model] provides the precomputed fixed point for the model. This version of the
routine initializes the path beyond the lagged values to the fp[model] values.

@o stochProto.c -d
@{

/*void * calloc(unsigned num,unsigned int amt);*/

void compXEtm1(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
void (* vecfunc)(),void (* fdjac)(),double * params,double * shockVec,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,
double x[],
unsigned int *check)
{
unsigned int i;
unsigned int aOne[1];
aOne[0]=1;
for(i=0;i<*numberOfEquations* *leads;i++)
qMat[i]=0.0;
for(i=0;i<*numberOfEquations* *leads;i++){
qMat[i]=1.0;
qMatj[i]=i+1+(*numberOfEquations* *lags);
qMati[i]=i+1;
}
qMati[i]=i+1;

for(i=0;i<*numberOfEquations* *leads;i++){
  fixedPoint[*numberOfEquations* (*lags)+i]=x[*numberOfEquations* (*lags+1)+i];}

pathNewt(numberOfEquations,lags, leads,aOne,
vecfunc, fdjac,params,shockVec,
fmats, fmatsj, fmatsi,
smats, smatsj, smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,
x,
check);
}
@}

@o stochProto.c -d
@{
void generateNextXTMinusOne(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,unsigned int * shockIndex,
double * shockTable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,
double x[],
unsigned int *check)
{
double * shockVec;double * tailVec;
unsigned int i;
shockVec= (double *) calloc(*numberOfEquations,sizeof(double));
tailVec= (double *) calloc(*numberOfEquations*(*lags+*leads+1),sizeof(double));
for(i=0;i<*numberOfEquations;i++)shockVec[i]=0;
pathNewt(numberOfEquations,lags, leads,pathLength,
vecfunc, fdjac,params,shockVec,
fmats, fmatsj, fmatsi,
smats, smatsj, smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,
x,
check);
for(i=0;i<*numberOfEquations*(*lags+*leads+1);i++)tailVec[i]=x[i];
compXEtm1(numberOfEquations,lags,leads,
vecfunc,fdjac,params,shockTable+(*numberOfEquations*(*shockIndex-1)),
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
tailVec,x,check);
free(tailVec);free(shockVec);
}
@}

@o stochProto.c -d
@{
void generatePathX(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,
unsigned int * numberOfShocks,
unsigned int * shockIndices,
double * shockTable,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,
double x[],
unsigned int *check)
{
unsigned int i;
for(i=0;i<*numberOfShocks;i++){
printf("for given draw, computing for date=%d\n",i);
generateNextXTMinusOne(numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,shockIndices+i,
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,
x+(*numberOfEquations*i),check);
}}
@}

@d stochSim argument list
@{unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
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
unsigned int *failedQ
@}


@o stochProto.c -d
@{
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
unsigned int *failedQ)
{
/*unsigned int check[1]={0};*/
unsigned int * numberOfShocks;
unsigned int i,j;
numberOfShocks=(unsigned int *)calloc(1,sizeof(unsigned int));
*numberOfShocks=*tf-*t0+1;
for(i=0;i<*replications;i++){
printf("computing for draw=%d\n",i);
for(j=0;j<*numberOfEquations*(*lags+*numberOfShocks+*leads);j++){
x[*numberOfEquations*(*lags+*numberOfShocks)*i+j]=
     dataTable[*t0* *numberOfEquations +j];}
/*initialize path to last simulation path*/
generatePathX(numberOfEquations,lags,leads,pathLength,
vecfunc,fdjac,params,numberOfShocks,permVecs+i**numberOfShocks,
shockTable,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,
x+(*numberOfEquations*(*lags+*numberOfShocks)*i),failedQ+i);
};
free(numberOfShocks);
}
@}


@o stochProto.m
@{
aimType2[model_List,horizon_,init_]:=FixedPointList[
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
FixedPointList[
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
 FixedPoint[
  nxtGuess[lags[model],
  newFunc,
  drvFunc[model],
  blockMatrix[{{zeroMatrix[numEq*mleads,numEq*(mlags)],IdentityMatrix[numEq*mleads]}}],
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
theZeroMatsC=Table[zeroMatrix[neq,neq],{nlag}],
theZeroMatsD=Table[zeroMatrix[neq,1],{nlag}]
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
blockMatrix[{{termConstr,zeroMatrix[neq]}}],
termConstr . (Transpose[{guess[[Range[-neq*(nlag+nlead),-1]]]}]- 
Transpose[{fp[[Range[-neq*(nlag+nlead),-1]]]}]),nxlstC,nxlstD}];
bs=backSub[lstC,lstD];
guess-Flatten[bs]]]]]]]]]

@}
This routine applies nxtGuessTerror to 
compute the time t 
@o stochProto.m
@{

aimType2Terror[shock_,model_, horizon_, init_] := FixedPointList[
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
FixedPointList[
{#[[1]]+1,aimType2[model,#[[1]]+1,init][[-1,numEq+Range[numEq]]]}&,{0,Table[0,{numEq}]},
SameTest->(Max[Abs[#1[[2]]-#2[[2]]]]<10^(-15)&)]];

altAimType3[model_,initHoriz_,init_]:=
With[{numEq=Length[func[model][[2]]]},
FixedPointList[
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

paperQ=blockMatrix[{{qmat,zeroMatrix[5]}}];


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

@o stochProto.h -d
@{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*
#include<time.h>
#include<sys/time.h>
*/
/* Performance Analysis simple code */
/*https://paolozaino.wordpress.com/2015/06/13/c-code-snippet-to-measure-function-execution-time-for-both-linux-and-mac-os-x/*/
#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#ifdef __APPLE__

  #define INIT_TIME struct timespec tsi, tsf; \
      double elaps_s; long elaps_ns; \
      clock_serv_t cclock; \
      mach_timespec_t mts;

#else

  #define INIT_TIME struct timespec tsi, tsf; \
      double elaps_s; long elaps_ns;

#endif

#ifdef __APPLE__

  #define START_GET_THE_TIME \
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock); \
    clock_get_time(cclock, &mts); \
    mach_port_deallocate(mach_task_self(), cclock); \
    tsi.tv_sec = mts.tv_sec; \
    tsi.tv_nsec = mts.tv_nsec;

  #define STOP_GET_THE_TIME \
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock); \
    clock_get_time(cclock, &mts); \
    mach_port_deallocate(mach_task_self(), cclock); \
    tsf.tv_sec = mts.tv_sec; \
    tsf.tv_nsec = mts.tv_nsec; \
    elaps_s = difftime(tsf.tv_sec, tsi.tv_sec); \
    elaps_ns = tsf.tv_nsec - tsi.tv_nsec;

#else

  #ifdef CLOCK_PROCESS_CPUTIME_ID
    /* cpu time in the current process */

    #define START_GET_THE_TIME clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tsi);

    #define STOP_GET_THE_TIME clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tsf); \
        elaps_s = difftime(tsf.tv_sec, tsi.tv_sec); \
        elaps_ns = tsf.tv_nsec - tsi.tv_nsec;

  #else

    /* this one should be appropriate to avoid errors on multiprocessors systems */

    #define START_GET_THE_TIME clock_gettime(CLOCK_MONOTONIC_RAW, &tsi);

    #define STOP_GET_THE_TIME clock_gettime(CLOCK_MONOTONIC_RAW, &tsf); \
        elaps_s = difftime(tsf.tv_sec, tsi.tv_sec); \
        elaps_ns = tsf.tv_nsec - tsi.tv_nsec;

  #endif

#endif




/*set maximal dimension constants*/
#define PATHLENGTH 1000
#define REPLICATIONS 1000
#define PRINTMAX 12
#define FALSE 0
#define TRUE 1


#define widthIntControlInfo 70
#define widthDoubleControlInfo 50
#define widthIntOutputInfo 20
#define widthDoubleOutputInfo 10
unsigned int intControlParameters[widthIntControlInfo];
double doubleControlParameters[widthDoubleControlInfo];
unsigned int intOutputInfo[REPLICATIONS*widthIntOutputInfo];
double doubleOutputInfo[REPLICATIONS*widthDoubleOutputInfo];


#define useStackQ intControlParameters[0]
#define useLnsrchQ intControlParameters[1]
#define numberOfDebugPairs intControlParameters[2]
#define useIdentityQ intControlParameters[3]
#define maxitsInput intControlParameters[4]
#define maxNumberStackElements intControlParameters[5]
#define maxNumberAimElements intControlParameters[6]
#define numberAlphas intControlParameters[7]
#define numberBetas intControlParameters[8]
#define useFirstDiffQ intControlParameters[9]
#define ma50PivotSearch intControlParameters[10]
#define useQQ intControlParameters[11]
#define terminalConstraintSelection intControlParameters[12]
#define useFixedPoint  0
#define useCrawlingDataPoint  1
#define useWaggingTail  2
#define useFixedDataPoint  3
#define useTailSolutionPlus  4
#define numberVarsToMonitor intControlParameters[13]
#define monitoredVars intControlParameters[14]
/*also reserve next 9 for monotored Vars*/
#define type3Q intControlParameters[24]
#define debugQ intControlParameters[25]
#define numberVarsToShock intControlParameters[26]
#define shockedVars intControlParameters[27]
/*also reserve next 9 for shocked Vars*/
#define tMinusOneQ intControlParameters[37]
#define debugPairs intControlParameters[38]
/*also reserve next 10 for debug pairs*/

#define homotopyXGuess intControlParameters[48]
#define homotopyEasy intControlParameters[49]
#define useBigX 0
#define useBigEasy 1
#define usePreviousHomotopyQ intControlParameters[50]
#define useShockFileQ intControlParameters[51]
#define shockFileOffset intControlParameters[52]
#define dataFileOffset intControlParameters[53]
#define streamingQ intControlParameters[54]
#define ICshockVecLength  intControlParameters[55]
#define ICnumberOfEquation intControlParameters[56]
#define ICnumberOfLags intControlParameters[57]
#define ICnumberOfLeads intControlParameters[58]
#define ICnumberOfParameters intControlParameters[59]
#define ICnumberOfDataValues intControlParameters[60]
#define ICnumberOfShocks intControlParameters[61]
#define ICnumberExog intControlParameters[62]
#define ignoreFailQ intControlParameters[63]

#define tolxInput doubleControlParameters[0]
#define tolfInput doubleControlParameters[1]
#define shrinkFactorInput doubleControlParameters[2]
#define expandFactorInput doubleControlParameters[3]
#define alaminInput doubleControlParameters[4]
#define alfInput doubleControlParameters[5]
#define ma50DropTol doubleControlParameters[6]

#define homotopyAlpha (doubleControlParameters+10)
/*also reserve next 9 for homotopyAlpha*/
#define homotopyBeta (doubleControlParameters+20)
/*also reserve next 9 for homotopyBeta*/
#define ma50Balance doubleControlParameters[30]
#define ma50DropEntry doubleControlParameters[31]
#define ma50DropCol doubleControlParameters[32]
#define shockScalar doubleControlParameters[33]



#define addOneToFailedQ (intOutputInfo[0])++
#define subOneFromFailedQ (intOutputInfo[0])--
#define resetFailedQ (intOutputInfo[0]=0)
#define addOneToNewtonSteps (intOutputInfo[1])++
#define resetNewtonSteps (intOutputInfo[1]=0)
#define addOneToFEvals (intOutputInfo[2])++
#define resetFEvals (intOutputInfo[2]=0)
#define addOneToFDrvEvals (intOutputInfo[3])++
#define resetFDrvEvals (intOutputInfo[3]=0)
#define addOneToShrinkSteps (intOutputInfo[4])++
#define resetShrinkSteps (intOutputInfo[4]=0)
#define addOneToExpandSteps (intOutputInfo[5])++
#define resetExpandSteps (intOutputInfo[5]=0)
#define addOneToLnsrchSteps (intOutputInfo[6])++
#define resetLnsrchSteps (intOutputInfo[6]=0)
#define addOneToHomotopies (intOutputInfo[7])++
#define resetHomotopies (intOutputInfo[7]=0)
#define addOneToHomotopyFailures (intOutputInfo[8])++
#define resetHomotopyFailures (intOutputInfo[8]=0)
#define currentReplication (intOutputInfo[9])
#define currentDate (intOutputInfo[10])

#define assignRealizedTolf doubleOutputInfo[0]
#define resetRealizedTolf (doubleOutputInfo[0]=0)
#define assignRealizedTolx doubleOutputInfo[1]
#define resetRealizedTolx (doubleOutputInfo[1]=0)


#define tMaOne (-1)
#define tMaTwo (-2)
#define tMaThree (-3)
#define tMaFour (-4)
#define tMaFive (-5)
#define tMaSix (-6)
#define tMaSeven (-7)
#define tMaEight (-8)
#define tMaNine (-9)
#define tMaTen (-10)
#define tMaEleven (-11)
#define tMaTwelve (-12)
#define tMaThirteen (-13)
#define tMaFourteen (-14)
#define tMaFifteen (-15)
#define tMaSixteen (-16)
#define tMaSeventeen (-17)
#define tMaEighteen (-18)
#define tMaNineteen (-19)
#define tMaTwenty (-20)
#define tMaTwentyOne (-21)
#define tMaTwentyTwo (-22)
#define tMaTwentyThree (-23)
#define tMaTwentyFour (-24)
#define tMaTwentyFive (-25)

#define tPaOne (1)
#define tPaTwo (2)
#define tPaThree (3)
#define tPaFour (4)
#define tPaFive (5)
#define tPaSix (6)
#define tPaSeven (7)
#define tPaEight (8)
#define tPaNine (9)
#define tPaTen (10)
#define tPaEleven (11)
#define tPaTwelve (12)
#define tPaThirteen (13)
#define tPaFourteen (14)
#define tPaFifteen (15)
#define tPaSixteen (16)
#define tPaSeventeen (17)
#define tPaEighteen (18)
#define tPaNineteen (19)
#define tPaTwenty (20)
#define tPaTwentyOne (21)
#define tPaTwentyTwo (22)
#define tPaTwentyThree (23)
#define tPaTwentyFour (24)
#define tPaTwentyFive (25)








@}

@o stochProto.h -d
@{

/*declare prototypes for functions*/
/*float  dtime_(float * userSystemTime);*/
double atof();
#include <stdio.h>
#include <stdlib.h>






/*unsigned int qRows=0;
 unsigned int auxInit=0;
 unsigned int aZero=0;*/

/*processCommandLine() determines defaults for these*/


/*timing routine variables*/


/*workspace*/
//double **fmats;unsigned int  **fmatsj;unsigned int  **fmatsi;
//double **smats;unsigned int  **smatsj;unsigned int  **smatsi;
/*success indicators for stochSims*/

/*csr q matrix*/
/*double * AMqMatrix;
unsigned int * AMqMatrixj;
unsigned int * AMqMatrixi;*/
/*csr b matrix*/
/*double * AMbMatrix;
unsigned int * AMbMatrixj;
unsigned int * AMbMatrixi;
 double * rootr;
 double * rooti;*/
/* double * brootr;
 double * brooti;
double*upsilonMatrix;unsigned int*upsilonMatrixj;unsigned int*upsilonMatrixi;
double*hMat;unsigned int*hMatj;unsigned int*hMati;
double*hzMat;unsigned int*hzMatj;unsigned int*hzMati;
double*cstar;unsigned int*cstarj;unsigned int*cstari;
double*phiInvMat;unsigned int*phiInvMatj;unsigned int*phiInvMati;
double*fmat;unsigned int*fmatj;unsigned int*fmati;
double*impact;unsigned int*impactj;unsigned int*impacti;*/
/*double*selectZmat;unsigned int*selectZmatj;unsigned int*selectZmati;*/
/*double*varthetaC;unsigned int*varthetaCj;unsigned int*varthetaCi;*/
/*double*varthetaZstar;unsigned int*varthetaZstarj;unsigned int*varthetaZstari;*/

 /*unsigned int * ma50bdJob;
 unsigned int * ma50bdIq;
 double * ma50bdFact;
 unsigned int * ma50bdIrnf;
 unsigned int * ma50bdIptrl;
 unsigned int * ma50bdIptru;*/
 /*unsigned int * cmpma50bdJob;
 unsigned int * cmpma50bdIq;
 double * cmpma50bdFact;
 unsigned int * cmpma50bdIrnf;
 unsigned int * cmpma50bdIptrl;
 unsigned int * cmpma50bdIptru;*/
/* unsigned int sysDim;*/




@}


@o stochProto.h -d
@{


//double FMAX(double a,double b);
//double FMIN(double a,double b);
double FABS(double a);
double doRightSmaller(double a,double b);
double doSign(double a);




@}


@o stochProto.h -d
@{

#ifdef __APPLE__
#include<strings.h>
#endif
#ifdef __linux__
#include<string.h>
#endif
/* */

void modData(unsigned int numberOfEquations,unsigned int numberDataValues,double * dataVals,
			 unsigned int vbl,unsigned int t0,unsigned int tf,double val1,double val2);

void modDataAbs(unsigned int numberOfEquations,unsigned int numberDataValues,double * dataVals,
			unsigned  int vbl,unsigned int t0,unsigned int tf,double val1,double val2);
/* */



void stochSim(@<stochSim argument list@>);
void generateDraws(unsigned int t0Index,unsigned int tfIndex,unsigned int replications,unsigned int shocksAvailable,
unsigned int * iarray);

void processCommandLine(unsigned int argc, char * argv[],char ** namesArray,unsigned int modelNEQS,char ** paramNamesArray,unsigned int numberOfParameters,double * parameters,
double * dataValues,unsigned int numberDataValues,unsigned int numShockValues,
unsigned int * pathLength,unsigned int * replications,unsigned int * t0,unsigned int * stochasticPathLength,
unsigned int * intControlParameters,double* doubleControlParameters,char * flnm);

@}

@o stochProto.c -d
@{
#ifdef __APPLE__
#include<strings.h>
#endif
#ifdef __linux__
#include<string.h>
#endif
/* */

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



void processCommandLine(unsigned int argc, char * argv[],char ** namesArray,unsigned int modelNEQS,char ** paramNamesArray,unsigned int numberOfParameters,double * parameters,
double * dataValues,unsigned int numberDataValues,unsigned int numShockValues,
unsigned int * pathLength,unsigned int * replications,unsigned int * t0,unsigned int * stochasticPathLength,
unsigned int * intControlParameters,double* doubleControlParameters,char * flnm)
{
  float aFloat;unsigned int i;unsigned int anInt;
 unsigned int pl;unsigned int t1;unsigned int t2; double val1; double val2; unsigned int vbl;
/*setup defaults*/
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
while((strcmp(argv[2],paramNamesArray[i])) && (i <modelNEQS))i++;
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
     if(pl>numShockValues||pl<0)
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
         sscanf(argv[2],"%d",&anInt);
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
         sscanf(argv[2],"%d",&anInt);
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

@}



\appendix


\include{stochRes}

\bibliographystyle{plain}
\bibliography{files}
\nocite{hollinger96}
\end{document}
