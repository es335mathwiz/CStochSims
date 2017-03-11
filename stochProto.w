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
void free();
void pathNewt(int * numberOfEquations,int * lags, int * leads,int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,double * shockVec,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,double * qMat,int * qMatj,int * qMati,
double * fixedPoint,
double x[],
int *check);
long ignuin_(long *low,long *high);


void generateDraws(int t0Index,int tfIndex,int replications,int shocksAvailable,
int * iarray)
{
static  long K1=1;
int ntot,i;
long mxint;
ntot=(tfIndex-t0Index+1)*replications;
mxint=shocksAvailable;
    for(i=0; i<ntot; i++) {
        *(iarray+i) = (int )ignuin_(&K1,&mxint);
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
#include <stdlib.h>
#include <stdio.h>

/*void * calloc(unsigned num,int amt);*/

void compXEtm1(int * numberOfEquations,int * lags, int * leads,
void (* vecfunc)(),void (* fdjac)(),double * params,double * shockVec,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,double * qMat,int * qMatj,int * qMati,
double * fixedPoint,
double x[],
int *check)
{
int i;
int aOne[1];
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
int * numberOfEquations,int * lags, int * leads,int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,int * shockIndex,
double * shockTable,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,double * qMat,int * qMatj,int * qMati,
double * fixedPoint,
double x[],
int *check)
{
double * shockVec;double * tailVec;
int i;
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
int * numberOfEquations,int * lags, int * leads,int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,
int * numberOfShocks,
int * shockIndices,
double * shockTable,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,double * qMat,int * qMatj,int * qMati,
double * fixedPoint,
double x[],
int *check)
{
int i;
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
@{int * numberOfEquations,int * lags, int * leads,int * pathLength,
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
int *failedQ
@}

@o stochProto.h -d
@{
void stochSim(@<stochSim argument list@>);
@}

@o stochProto.c -d
@{
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
int *failedQ)
{
/*int check[1]={0};*/
int * numberOfShocks;
int i,j;
numberOfShocks=(int *)calloc(1,sizeof(int));
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
