%	$Id: stackC.w,v 1.39 2000/03/28 21:05:12 m1gsa00 Exp m1gsa00 $	
\documentclass{article}
%\documentclass[html]{article}
@i miscLatexPkg.tex
\newcommand{\Treebox}[1]{
\Tr{\psframebox{#1}}}
\begin{document}

\title{``C'' Implementation of Stack:
}
\author{Gary Anderson}
\maketitle
 
\centerline{{\bf Abstract}}
\begin{quote}
  This paper describes the ``C'' implemention of the Stack Algorithm  outlined in\cite{juillard96}.
The document describes three components: the ``C'' subroutines for implementing
the stack algorithm, the makefile for creating the executables and libraries, 
and the ``C'' for a testing subroutines.
The code uses SPARSEKIT2 and HARWELL code for sparse matrices.
\end{quote}

\newpage
\tableofcontents
\newpage

This paper describes the ``C'' implemention of the Stack Algorithm  outlined in\cite{juillard96}.
The code uses SPARSEKIT2\cite{saad94}  and HARWELL\cite{nag95} code for sparse matrices.


\begin{description}
\item[{\bf nxtCDmats}]Given the current C and D matrices and
the curren H matrix, this routine computes the next C and D matrices.
\end{description}



The document describes three components: the ``C'' subroutines for implementing
the stack algorithm, the makefile for creating the executables and libraries, 
and the ``C'' for  testing the subroutines.


The function computes the next $C$ and $d$ matrices.
It assumes that the caller has allocated enough space to hold the pointers to the 
two matrices.




Figure \ref{tableau} presents a graphic characterization of the relevant
sets of linear constraints.
\begin{figure*}[htbp]
  \begin{center}
    \leavevmode
    
\fbox{
  \begin{pspicture}(14,5)
\rput[bl](0,0){\rnode{A}{\mbox{
\begin{pspicture}(6,5)
\psframe[fillstyle=solid,fillcolor=gray](0,5)(3,4)
\rput(0.5,4.5){$I$}
\psline[linewidth=0.6mm](1,4)(1,5)
\rput(2,4.5){$C_{-\tau}$}
\psframe[fillstyle=solid,fillcolor=gray](1,4)(4,3)
\rput(1.5,3.5){$I$}
\psline[linewidth=0.6mm](2,3)(2,4)
\rput(3,3.5){$C_{-\tau+1}$}
\psframe[fillstyle=solid,fillcolor=gray](2,2)(5,1)
\rput(2.5,1.5){$I$}
\psline[linewidth=0.6mm](3,1)(3,2)
\rput(4,1.5){$C_{-1}$}
\psframe[fillstyle=solid,fillcolor=gray](0,1)(6,0)
\rput(1.5,0.5){$H_{-}$}
\psline[linewidth=0.6mm](3,1)(3,0)
\rput(3.5,0.5){$H_{0}$}
\psline[linewidth=0.6mm](4,1)(4,0)
\rput(5,0.5){$H_{+}$}
%\psframe[fillstyle=solid,fillcolor=gray](9,0)(15,1)
\psline[linestyle=dashed]{<->}(1.5,3)(2.5,2)
%\psline[linestyle=dashed]{<->}(9,3.5)(11,1.5)
\psline[linewidth=0.6mm](0,1)(6,1)
%\psline[linewidth=0.6mm](5,0)(5,10)
%\psline[linewidth=0.6mm](10,0)(10,10)
\psset{arrows=->}
%\rput(7.5,8.0){\rnode{A}{$\tau+\theta$}}
%\rput(5,8.0){\rnode{B}{}}
%\rput(10,8.0){\rnode{C}{}}
%\ncline{A}{B}
%\ncline{A}{C}
%\rput(10.5,5.0){\rnode{D}{}}
%\rput(12.5,7.5){\rnode{E}{$H^{\sharp,0}$}}
%\nccurve[angleB=90,angleA=270]{E}{D}
\end{pspicture}}}}
\rput[bl](7,0){\rnode{B}{\mbox{
\begin{pspicture}(6,5)
\psframe[fillstyle=solid,fillcolor=gray](0,5)(3,4)
\rput(0.5,4.5){$I$}
\psline[linewidth=0.6mm](1,4)(1,5)
\rput(2,4.5){$C_{-\tau}$}
\psframe[fillstyle=solid,fillcolor=gray](1,4)(4,3)
\rput(1.5,3.5){$I$}
\psline[linewidth=0.6mm](2,3)(2,4)
\rput(3,3.5){$C_{-\tau+1}$}
\psframe[fillstyle=solid,fillcolor=gray](2,2)(5,1)
\rput(2.5,1.5){$I$}
\psline[linewidth=0.6mm](3,1)(3,2)
\rput(4,1.5){$C_{-1}$}
\psframe[fillstyle=solid,fillcolor=gray](3,1)(6,0)
\rput(3.5,0.5){$I$}
\psline[linewidth=0.6mm](4,0)(4,1)
\rput(5,0.5){$C_{0}$}
%\psframe[fillstyle=solid,fillcolor=gray](9,0)(15,1)
\psline[linestyle=dashed]{<->}(1.5,3)(2.5,2)
%\psline[linestyle=dashed]{<->}(9,3.5)(11,1.5)
\psline[linewidth=0.6mm](0,1)(6,1)
%\psline[linewidth=0.6mm](5,0)(5,10)
%\psline[linewidth=0.6mm](10,0)(10,10)
%\psset{arrows=->}
%\rput(7.5,8.0){\rnode{A}{$\tau+\theta$}}
%\rput(5,8.0){\rnode{B}{}}
%\rput(10,8.0){\rnode{C}{}}
%\ncline{A}{B}
%\ncline{A}{C}
%\rput(10.5,5.0){\rnode{D}{}}
%\rput(12.5,7.5){\rnode{E}{$H^{\sharp,0}$}}
%\nccurve[angleB=90,angleA=270]{E}{D}
\end{pspicture}}}}
\psset{arrows=->}
\ncline{A}{B}
  \end{pspicture}}
    \caption{Matrix Tableau Characterization of Algorithm: Initial Tableau}
    \label{tableau}
  \end{center}
\end{figure*}
In the figure the regions where the coefficients are potentially  non-zero are shaded gray. 


\section{An Example}







\subsection{Overview}
\label{sec:overview}



This example shows how to computes the path update for the model
presented in Julliard\cite{juillard96}. The example uses
the Anderson-Moore terminal constraints in place of the terminal identity
matrix used in the original Julliard paper.

This implementation of the  Stack algorithm requires two function.
\begin{description}
\item[{\bf fFunc}] Computes the value of $f(y,\beta)$
\item[{\bf dfFunc}] Computes the value of $\frac{\partial f(y,\beta)}{\partial y}$
\end{description}


\subsection{Include Files}
\label{sec:include}








\section{Stack Algorithm ``C'' Source Code Generation}
\label{sec:stackSource}

Assemble the components and output to the file {\bf stackC.c}.

@o  stackC.c -d
@{
@<define assert bump@>
@<define constants and specify include files@>
@<nxtCDmats definition@>
@<oneStepBack definition@>

@<compPathError definition@>
@<nxtGuess definition@>
@<nxtFPGuess definition@>
@<chkDrv definition@>
@}



\subsection{Defines and Includes}
\label{sec:defines}

@d define constants and specify include files
@{
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "useSparseAMA.h"
#include <sys/types.h>
#include <unistd.h>
@|stdio.h math.h
@}



\subsection{nxtCDmats Definition}
\label{sec:nxtCDmats}

The matrix represented by $(smatsA[0],smatsJA[0],smatsIA[0])$ contains the
$\frac{\partial f}{\partial x}$ for current time.\footnote{
Julliard denotes these by S.}
The matrices represented by $(oddSumCA,oddSumCJA,oddSumCIA)$ and
$(evenSumCA,evenSumCJA,evenSumCIA)$ represent the sum for the $C$
 computation at different
stages of the calculation.
\begin{gather*}
  \begin{bmatrix}
    H_{-\tau}&\dots&H_{\theta}
  \end{bmatrix}
\end{gather*}
The matrix represented by $(smatsA,smatsJA,smatsIA)$ contains the
$\frac{\partial f}{\partial x}$ for current time.\footnote{
Julliard denotes these by S.}
The matrices represented by $(oddSumDA,oddSumDJA,oddSumDIA)$ and
$(evenSumDA,evenSumDJA,evenSumDIA)$ represent the sum for the $d$ computation
at different
stages of the calculation.

Using odd and even to minimize need for allocating space for partial
sums. Beginning with total in $oddSum$.

\begin{programbox}
  |(oddSumCA,oddSumCJA,oddSumCIA)| =  \begin{bmatrix}H_{-\tau}&\dots&H_{\theta}  \end{bmatrix}
  |(oddSumDA,oddSumDJA,oddSumDIA)|=  \begin{bmatrix}f  \end{bmatrix}
\end{programbox}




With weighted sum of S matrices in hand(in the odd version of the variables)
Use the Harwell $MA50$ routines to factorize the matrices.

Since Harwell expects Compressed Sparse Columen(CSC) 
instead of Compressed Sparse Row(CSR) we'll have to
transpose when backsolving.

HARWELL documentation suggests setting SPARSEFACTOR to 3.
@d define constants and specify include files
@{
#define SPARSEFACTOR 3 @|  SPARSEFACTOR 
@}


@o stackC.h -d
@{
void nxtCDmats(@<nxtCDmats argument list@>);
@}
Function uses SPARSEKIT's CSR format.
Since HARWELL functions expects CSC,
MA50CD operates on the transpose.

@d nxtCDmats definition
@{
void nxtCDmats(@<nxtCDmats argument list@>){
#include "useSparseAMA.h"
#include "stackC.h"
@<nxtCDmats variable declarations@>
@<nxtCDmats scalar variable allocations@>
@<nxtCDmats array variable allocations@>
@<ma50xx array variable allocations@>
@<initialize sum with original fvec and smat@>
for(i=0;i<*lagss; i++){
@<loop over lagged C and d matrices@>
}
@<compute pivot sequence@>
@<factorize matrix@>
@<use factorization@>
@< nxtCDmats scalar variable deallocations@>
@< nxtCDmats array variable deallocations@>
@< ma50xx array variable deallocations@>
}
@}

@d nxtCDmats argument list
@{
unsigned int * numberOfEquations, unsigned int * lagss, unsigned int * leadss,
unsigned int * rowDim,
unsigned int * maxNumberHElements,
double * smatsA,unsigned int * smatsJA,unsigned int * smatsIA,
double * fvecA,unsigned int * fvecJA,unsigned int * fvecIA,
double ** cmatsA,unsigned int **cmatsJA,unsigned int **cmatsIA,
double **dmatsA,unsigned int **dmatsJA,unsigned int **dmatsIA
@| numberOfEquations lagss leadss 
smatsA smatsJA smatsIA 
fvecA fvecJA fvecIA 
cmatsA cmatsJA cmatsIA
dmatsA dmatsJA dmatsIA
@}



@d initialize sum with original fvec and smat
@{
for(i=0;i<smatsIA[*rowDim]-smatsIA[0];i++){
oddSumCA[i]=smatsA[i];
oddSumCJA[i]=smatsJA[i];}
for(i=0;i<*rowDim+1;i++){oddSumCIA[i]=smatsIA[i];}

for(i=0;i<*rowDim;i++)
{oddSumDA[i]=fvecA[i];oddSumDJA[i]=fvecJA[i];}
for(i=0;i<*rowDim+1;i++)
{oddSumDIA[i]=fvecIA[i];}
@}



@d loop over lagged C and d matrices
@{
timeOffset=i-*lagss;
@<obtain an rowDim by numberOfEquations sub matrix of s matrices@>
@<multiply c matrices by appropriate s matrix and subtract@>
@<multiply d matrices by appropriate s matrix and subtract@>
@<switch odd for even to avoid calloc@>
@}




% \psset{arrows=->}
% \pstree{\Tcircle{smats}}{\pstree{\Treebox{SUBMAT\_}}{\Tcircle{ao}}}

\vspace{1.0cm}

\begin{programbox}
  |(ao,jao,iao)| =  H_{-\tau+i}
\end{programbox}



@d obtain an rowDim by numberOfEquations sub matrix of s matrices
@{
*firstColumn= 1+(i * *numberOfEquations);
*lastColumn= *firstColumn+ *numberOfEquations-1;
extractSubmatrix(rowDim,aOne,aOne,rowDim,firstColumn,lastColumn,
oddSumCA,oddSumCJA,oddSumCIA,nr,nc,ao,jao,iao);
@}



\vspace{1.0cm}



\begin{programbox}
  |(b,jb,ib)| =  H_{-\tau+i} C_{-\tau+i}
  |(evenSumCA,evenSumCJA,evenSumCIA)| =  |(oddSumCA,oddSumCJA,oddSumCIA)| - |(b,jb,ib)|
\end{programbox}



@d multiply c matrices by appropriate s matrix and subtract
@{
printf("rowDim=%u,cColumns=%u,nzmax=%u,a=\n",*rowDim,*cColumns,*nzmax);
cPrintSparse(*rowDim,ao,jao,iao);
printf("rowDim=%u,cColumns=%u,nzmax=%u,a=\n",*rowDim,*cColumns,*nzmax);
printf("timeOffset=%d\n",timeOffset);
cPrintSparse(*rowDim,(cmatsA[timeOffset]),(cmatsJA[timeOffset]),(cmatsIA[timeOffset]));
sparseMult(rowDim,cColumns,nzmax,iw,aOne,ao,jao,iao,
(cmatsA[timeOffset]),(cmatsJA[timeOffset]),(cmatsIA[timeOffset]),
b,jb,ib,ierr);

pathNewtAssert(*ierr == 0);
bump((cmatsIA[timeOffset])[*rowDim]-(cmatsIA[timeOffset])[0]);
aSmallDouble=DBL_EPSILON;
dropSmallElements(rowDim,aOne,&aSmallDouble,nzmax,b,jb,ib,b,jb,ib,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);
/*actually want to subtract so mult elements by -1 also need to shift right*/
for(j=0;j<ib[*rowDim]-1;j++)
{b[j]=(-1)*b[j];jb[j]=jb[j]+(*numberOfEquations*(timeOffset+*lagss+1));};

printf("odd\n");
cPrintSparse(*rowDim,oddSumCA,oddSumCJA,oddSumCIA);
printf("b\n");
cPrintSparse(*rowDim,b,jb,ib);
sparseAdd(rowDim,cColumns,&nzmax,iw,aOne,oddSumCA,oddSumCJA,oddSumCIA,
b,jb,ib,evenSumCA,evenSumCJA,evenSumCIA,ierr);
printf("after sparseAdd\n");
pathNewtAssert(*ierr == 0);
bump(evenSumCIA[*rowDim]-evenSumCIA[0]);
@}



\vspace{1.0cm}


\begin{programbox}
  |(b,jb,ib)| =  H_{-\tau+i} C_{-\tau+i}
  |(evenSumDA,evenSumDJA,evenSumDIA)| =  |(oddSumDA,oddSumDJA,oddSumDIA)| - |(b,jb,ib)|
\end{programbox}



@d multiply d matrices by appropriate s matrix and subtract
@{

sparseMult(rowDim,aOne,nzmax,iw,aOne,ao,jao,iao,
(dmatsA[timeOffset]),(dmatsJA[timeOffset]),(dmatsIA[timeOffset]),
b,jb,ib,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);
aSmallDouble=DBL_EPSILON;
dropSmallElements(rowDim,aOne,&aSmallDouble,nzmax,b,jb,ib,b,jb,ib,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);

/*actually want to subtract so mult elements by -1*/
for(j=0;j<ib[*rowDim]-1;j++)b[j]=(-1)*b[j];

sparseAdd(rowDim,aOne,nzmax,iw,aOne,oddSumDA,oddSumDJA,oddSumDIA,
b,jb,ib,
evenSumDA,evenSumDJA,evenSumDIA,ierr);
pathNewtAssert(*ierr == 0);
bump(evenSumDIA[*rowDim]-evenSumDIA[0]);


@}

@d switch odd for even to avoid calloc
@{
tmp=oddSumCA;jtmp=oddSumCJA;itmp=oddSumCIA;
oddSumCA=evenSumCA;oddSumCJA=evenSumCJA;oddSumCIA=evenSumCIA;
evenSumCA=tmp;evenSumCJA=jtmp;evenSumCIA=itmp;
tmp=oddSumDA;jtmp=oddSumDJA;itmp=oddSumDIA;
oddSumDA=evenSumDA;oddSumDJA=evenSumDJA;oddSumDIA=evenSumDIA;
evenSumDA=tmp;evenSumDJA=jtmp;evenSumDIA=itmp;
@}



% \psset{arrows=<-}
% \pstree[treemode=U]{\Tcircle{oddsumD}}{\pstree{\Treebox{SUBMAT\_}}{\Tcircle{evensumD}\Tcircle{fvec}}}


\vspace{1.0cm}


Given a matrix $A$, MA50AD computes $P,Q$ such that
\begin{gather*}
  PAQ=LU
\end{gather*}

Given a matrix $A$, and $P,Q$ such that 
\begin{gather*}
  PAQ=LU
\end{gather*}
MA50BD computes the factorization.

@d compute pivot sequence
@{
/*still using CSR consequently doing everything to the 
transpose*/
/*copy submat of 
 oddSumC for subsequent use. note ma50ad modifies its A argument*/


for(i=0;i<*maxNumberHElements;i++)
{evenSumCA[i]=oddSumCA[i];evenSumCJA[i]=oddSumCJA[i];}
for(i=0;i<*rowDim +1;i++)
{evenSumCIA[i]=oddSumCIA[i];}
printf("at line 466\n");cPrintSparse(*rowDim,evenSumCA,evenSumCJA,evenSumCIA);
*firstColumn=(*numberOfEquations* *lagss)+1;
*lastColumn=*firstColumn + *rowDim-1;
extractSubmatrix(rowDim,aOne,aOne,rowDim,
firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,nr,nc,
oddSumCA,oddSumCJA,oddSumCIA);
*nonZeroNow=oddSumCIA[*rowDim]-oddSumCIA[0];
cPrintSparse(*rowDim,oddSumCA,oddSumCJA,oddSumCIA);
useMA50ID(cntl,icntl);
useMA50AD(rowDim,rowDim,nonZeroNow,nzmax,oddSumCA,oddSumCJA,jcn,oddSumCIA,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
/*wordybump(info[3]);*/
pathNewtAssert(info[0]>=0);
/* restore odd since ad is destructive*/
extractSubmatrix(rowDim,aOne,aOne,rowDim,
firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,nr,nc,
oddSumCA,oddSumCJA,jcn);

@}


@d factorize matrix
@{
useMA50BD(rowDim,rowDim,nonZeroNow,aOne,
oddSumCA,oddSumCJA,jcn,
cntl,icntl,ip,oddSumCIA,np,lfact,fact,irnf,iptrl,iptru,
w,iw,info,rinfo);
/*wordybump(info[3]);*/
pathNewtAssert(info[0]>=0);
@}

MA50CD applies the factoriation to solve
\begin{gather*}
  A^T x = b
\end{gather*}




@d use factorization

@{
    *trans = 1;

/*expand sum of c's use transpose since c colum major order */

itb[0]=1;cmatsExtent=0;
for(i=0;i<*cColumns;i++){


*lastColumn = *firstColumn=(1+i+*numberOfEquations *(1+*lagss));

extractSubmatrix(rowDim,aOne,
aOne,rowDim,firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,
nr,nc,
b,jb,ib);

csrToDns(rowDim,aOne,b,jb,ib,
nsSumC,rowDim,ierr);

useMA50CD(rowDim,rowDim,icntl,oddSumCIA,np,trans,
lfact,fact,irnf,iptrl,iptru,
nsSumC,x,
w,info);
bump(info[3]);
pathNewtAssert(info[0]>=0);

*nzmaxLeft= *nzmax-cmatsExtent-1;
dnsToCsr(aOne,rowDim,nzmaxLeft,x,
aOne,tb+(itb[i]-1),jtb+(itb[i]-1),itb+i,
ierr);
pathNewtAssert(*ierr == 0);

itb[i+1]=itb[i+1]+cmatsExtent;
itb[i]=itb[i]+cmatsExtent;
cmatsExtent=itb[i+1]-1;
}
bump(cmatsExtent);

aSmallDouble=DBL_EPSILON;
dropSmallElements(cColumns,aOne,&aSmallDouble,nzmax,tb,jtb,itb,tb,jtb,itb,ierr);
csrToCsc(cColumns,aOne,aOne,tb,jtb,itb,cmatsA[0],cmatsJA[0],cmatsIA[0]);

/*expand sum of d's*/
csrToDns(rowDim,aOne,oddSumDA,oddSumDJA,oddSumDIA,nsSumD,rowDim,ierr);
pathNewtAssert(*ierr == 0);
bump(*rowDim);
/*code should use info from previous call to set lfact
also can avoid calls to ma50ad once pattern settles down*/

useMA50CD(rowDim,rowDim,icntl,oddSumCIA,np,
trans,lfact,fact,irnf,iptrl,iptru,nsSumD,x,w,info);

dnsToCsr(rowDim,aOne,rowDim,x,
rowDim,dmatsA[0],dmatsJA[0],dmatsIA[0],ierr);
/*wordybump(info[3]);*/
pathNewtAssert(*ierr == 0);
@}

\subsection{nxtCDmats Variable Declarations}
\label{sec:nxtcddeclarations}


@d nxtCDmats variable declarations
@{
unsigned int maxElementsEncountered=0;
/*void * calloc(unsigned amt,unsigned int size);*/
double * evenSumCA;unsigned int * evenSumCJA;unsigned int * evenSumCIA;
double * evenSumDA;unsigned int * evenSumDJA;unsigned int * evenSumDIA;
double * oddSumCA;unsigned int * oddSumCJA;unsigned int * oddSumCIA;
double * oddSumDA;unsigned int * oddSumDJA;unsigned int * oddSumDIA;
double  *ao;unsigned int *jao,*iao;
double *b;unsigned int *jb,*ib;
double *tb;unsigned int *jtb,*itb;
double *tmp;unsigned int *jtmp,*itmp;
unsigned int *firstColumn,*lastColumn,*nr,*nc;
unsigned int *iw,*ierr,*nzmax,*nonZeroNow;unsigned int cmatsExtent;unsigned int *nzmaxLeft;
unsigned int i;
unsigned int j;
int timeOffset;
unsigned int *jcn;
double * cntl;
unsigned int * icntl;
unsigned int * ip ;
unsigned int * np;
unsigned int * jfirst;
unsigned int * lenr;
unsigned int * lastr;
unsigned int * nextr;
unsigned int * ifirst;
unsigned int * lenc;
unsigned int * lastc;
unsigned int * nextc;
unsigned int * info;
double * rinfo;
unsigned int *lfact;
double * fact;
unsigned int *irnf;
unsigned int * iptrl;
unsigned int * iptru;
unsigned int * aOne;
double aSmallDouble;
double * w;
double * x;
unsigned int * trans;
unsigned int * cColumns;
unsigned int * balColumns;
unsigned int *hColumns;
double * nsSumC;
double * nsSumD;@|
cPrintMatrixNonZero
evenSumCA   evenSumCJA   evenSumCIA 
  evenSumDA   evenSumDJA   evenSumDIA 
  oddSumCA   oddSumCJA   oddSumCIA 
  oddSumDA   oddSumDJA   oddSumDIA 
  firstColumn lastColumn nr nc jao iao 
 jb ib iw ierr jtmp itmp nzmax nonZeroNow 
 b tmp 
 ao 
i timeOffset 
 jcn 
  cntl 
  icntl 
  ip  
  np 
  jfirst 
  lenr 
  lastr 
  nextr 
  ifirst 
  lenc 
  lastc 
  nextc 
  info 
  rinfo 
 lfact 
  fact 
 irnf 
  iptrl 
  iptru 
  aOne 
  aTen 
  w 
  x 
  trans 
  cColumns 
  nsSumC hColumns
  nsSumD
@}


\subsection{nxtCDmats Storage Management}
\label{sec:nxtcdstorage}


@d nxtCDmats scalar variable allocations
@{

firstColumn = (unsigned int *)calloc(1,sizeof(unsigned int));
lastColumn = (unsigned int *)calloc(1,sizeof(unsigned int));
nr = (unsigned int *)calloc(1,sizeof(unsigned int));
nc = (unsigned int *)calloc(1,sizeof(unsigned int));
ierr = (unsigned int *)calloc(1,sizeof(unsigned int));
nzmax = (unsigned int *)calloc(1,sizeof(unsigned int));
nzmaxLeft = (unsigned int *)calloc(1,sizeof(unsigned int));
hColumns = (unsigned int *)calloc(1,sizeof(unsigned int));
nonZeroNow = (unsigned int *)calloc(1,sizeof(unsigned int));
cColumns = (unsigned int *)calloc(1,sizeof(unsigned int));
balColumns = (unsigned int *)calloc(1,sizeof(unsigned int));
*hColumns = *numberOfEquations*(*leadss+*lagss+1);
*cColumns = *numberOfEquations * *leadss;
*nzmax = *maxNumberHElements;
@|
nonZeroNow 
nzmax 
ierr 
nc 
nr 
lastColumn 
firstColumn 


@}

@d ma50xx array variable allocations
@{

/*for ma50ad*/
jcn = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
cntl= (double *)calloc(5,sizeof(double));
icntl= (unsigned int *)calloc(9,sizeof(unsigned int));
ip = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
np = (unsigned int *)calloc(1,sizeof(unsigned int));
jfirst = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
lenr = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
lastr = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
nextr = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
ifirst = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
lenc = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
lastc = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
nextc = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
info = (unsigned int *)calloc(7,sizeof(unsigned int));
rinfo = (double *)calloc(1,sizeof(double));
/* ma50bd*/
lfact =(unsigned int *)calloc(1,sizeof(unsigned int));
*lfact = ( *maxNumberHElements);/*pessimistic setting for filling*/
fact = (double *)calloc(*lfact,sizeof(double));
irnf = (unsigned int *)calloc(*lfact,sizeof(unsigned int));
iptrl = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
iptru = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
x = (double *)calloc(*rowDim * *numberOfEquations * (*leadss+1),sizeof(double));
trans = (unsigned int *) calloc(1,sizeof(unsigned int));
nsSumD = (double *)calloc(*rowDim,sizeof(double));
nsSumC = (double *)calloc(*rowDim /** *numberOfEquations * (*leadss+1+ *lagss)*/,sizeof(double));
aOne = (unsigned int *)calloc(1,sizeof(unsigned int));
*aOne = 1;
@}


@d nxtCDmats array variable allocations
@{

ib = (unsigned int *)calloc(*rowDim * *leadss + 1,sizeof(unsigned int));
jb = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
b = (double *)calloc(*maxNumberHElements,sizeof(double));
itb = (unsigned int *)calloc(*cColumns + 1,sizeof(unsigned int));
jtb = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
tb = (double *)calloc(*maxNumberHElements,sizeof(double));


evenSumCIA = (unsigned int *)calloc(( *rowDim * *leadss)+1,sizeof(unsigned int));
evenSumCJA = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
evenSumCA = (double *)calloc(*maxNumberHElements,sizeof(double));


/*larger than necessary now so that can use for transpose in csrcsc */
oddSumCIA = (unsigned int *)calloc(( *rowDim * *leadss)+ 1,sizeof(unsigned int));
oddSumCJA = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
oddSumCA = (double *)calloc(*maxNumberHElements,sizeof(double));


evenSumDIA = (unsigned int *)calloc(*rowDim+1,sizeof(unsigned int));/*MLK*/
evenSumDJA = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
evenSumDA = (double *)calloc(*rowDim,sizeof(double));


oddSumDIA = (unsigned int *)calloc(*rowDim+1,sizeof(unsigned int));/*MLK*/
oddSumDJA = (unsigned int *)calloc(*rowDim,sizeof(unsigned int));
oddSumDA = (double *)calloc(*rowDim,sizeof(double));


iao = (unsigned int *)calloc(*rowDim+1,sizeof(unsigned int));
jao = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
ao = (double *)calloc(*maxNumberHElements,sizeof(double));

/*work array needs elements equal to number of columns of matrix*/
iw = (unsigned int *)calloc(*rowDim * (*leadss + *lagss + 1),sizeof(unsigned int));
w = (double *)calloc(*rowDim * (*leadss + *lagss + 3),sizeof(double));/*MLK*/@|
oddSumCA 
oddSumCJA 
oddSumCIA 
evenSumDA 
evenSumDJA 
evenSumDIA 
evenSumCA 
evenSumCJA 
evenSumCIA 
b 
jb 
ib 

@}

@d nxtCDmats scalar variable deallocations
@{
free(firstColumn);
free(lastColumn);
free(nr);
free(nc);
free(ierr);
free(nzmax);
free(nzmaxLeft);
free(nonZeroNow);
free(hColumns);
free(cColumns);
free(balColumns);
@}


@d nxtCDmats array variable deallocations
@{
free(ib);
free(jb);
free(b);
free(itb);
free(jtb);
free(tb);
free(evenSumCIA);
free(evenSumCJA);
free(evenSumCA);
free(oddSumCIA);
free(oddSumCJA);
free(oddSumCA);
free(evenSumDIA);
free(evenSumDJA);
free(evenSumDA);
free(oddSumDIA);
free(oddSumDJA);
free(oddSumDA);
free(iao);
free(jao);
free(ao);
free(iw);
free(w);
@}

@d ma50xx array variable deallocations
@{

/*for ma50ad*/
free(jcn);
free(cntl);
free(icntl);
free(ip);
free(np);
free(jfirst);
free(lenr);
free(lastr);
free(nextr);
free(ifirst);
free(lenc);
free(lastc);
free(nextc);
free(info);
free(rinfo);
/* ma50bd*/
free(lfact);
free(fact);
free(irnf);
free(iptrl);
free(iptru);
free(x);
free(trans);
free(nsSumD);
free(nsSumC);
free(aOne);
@}


\subsection{backSub Definition}
\label{sec:backSub}
@d oneStepBack argument list
@{
unsigned int * rowDim,
double ** yvecA,unsigned int ** yvecJA,unsigned int ** yvecIA,
double ** cmatsA,unsigned int **cmatsJA,unsigned int **cmatsIA,
double **dmatsA,unsigned int **dmatsJA,unsigned int **dmatsIA
 @| numberOfEquations lagss leadss 
smatsA smatsJA smatsIA 
yvecA yvecJA yvecIA 
cmatsA cmatsJA cmatsIA
dmatsA dmatsJA dmatsIA
@}



\subsection{oneStepBack Definition}
\label{sec:oneStepBack}

@o stackC.h -d
@{
void oneStepBack(@<oneStepBack argument list@>);
@}



@d oneStepBack definition
@{
#include "stackC.h"
void oneStepBack(@<oneStepBack argument list@>)
{
@<oneStepBack variable definitions@>
@<oneStepBack variable allocations@>
/*cmat non zero then multiply else product is zero*/
if(cmatsIA[0][*rowDim]-cmatsIA[0][0]) {
  sparseMult(rowDim,aOne,rowDim,iw,aOne,cmatsA[0],cmatsJA[0],cmatsIA[0],
  yvecA[0+1 ],yvecJA[0+1 ],yvecIA[0+1 ],
  rcy,rcyj,rcyi,ierr);
pathNewtAssert(*ierr == 0);

aSmallDouble=DBL_EPSILON;

dropSmallElements(rowDim,aOne,&aSmallDouble,nzmax,rcy,rcyj,rcyi,rcy,rcyj,rcyi,ierr);
pathNewtAssert(*ierr == 0);

  for(i=0;i<rcyi[*rowDim]-rcyi[0];i++)rcy[i]=(-1)*rcy[i];
  sparseAdd(rowDim,aOne,nzmax,iw,aOne,
  dmatsA[0],dmatsJA[0],dmatsIA[0],
  rcy,rcyj,rcyi,
  yvecA[(0) ],yvecJA[(0) ],yvecIA[(0)],ierr);
pathNewtAssert(*ierr == 0);
} else {
  for(i=0;i<*rowDim;i++)
  {yvecA[0][i]=dmatsA[0][i];}
  for(i=0;i<*rowDim;i++)
  {yvecJA[0][i]=dmatsJA[0][i];}
  for(i=0;i<=*rowDim;i++)
  {yvecIA[0][i]=dmatsIA[0][i];}
  }



@<oneStepBack variable deallocations@>
}
@}
@d oneStepBack variable definitions
@{
double aSmallDouble;
unsigned int * aOne;
double *rcy;
unsigned int *rcyj ,*rcyi;
unsigned int *ierr;
unsigned int i;
unsigned int *nzmax;
unsigned int * iw;
@}
@d oneStepBack variable allocations
@{
nzmax=(unsigned int *)calloc(1,sizeof(unsigned int));
aOne=(unsigned int *)calloc(1,sizeof(unsigned int));
ierr=(unsigned int *)calloc(1,sizeof(unsigned int));
*aOne=1;
*nzmax=*rowDim;
rcy= (double *)calloc(*rowDim,sizeof(double));
rcyj=(unsigned int *) calloc(*rowDim,sizeof(unsigned int));
rcyi=(unsigned int *)calloc(*rowDim+1,sizeof(unsigned int));
iw = (unsigned int *)calloc(1,sizeof(unsigned int));
@}
@d oneStepBack variable deallocations
@{
free(ierr);
free(aOne);
free(rcy);
free(rcyj);
free(rcyi);
free(iw);
free(nzmax);
@}

\subsection{terminalConditions}
\label{sec:term}
@d computeR argument list
@{
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
unsigned int * maxNumberElements,
void (* func)(),void (* dfunc)(),double * params,
double * expansionPoint,
double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * impact, unsigned int * impactj, unsigned int * impacti,
double * rvec
@}
@d computeR
@{
void computeR(
@<computeR argument list@>
)
{
impactPart1=calloc(numberOfEquations*leads,sizeof(double));
impactPart2=calloc(numberOfEquations*leads,sizeof(double));
/*xxxxxxxxx add code for deviations*/
rowDim=numberOfEquations*leads;

sparseMatTimesVec(&rowDim,termConstr,termConstrj,termConstri,expansionPoint,rvec);
sparseMatTimesVec(&rowDim,
impactr,impactrj,impactri,expansionPoint,impactPart1+(numberOfEuations*lags));
sparseMatTimesVec(&rowDim,
impactr,impactrj,impactri,impactPart1,impactPart2);
for(i=0;i<numberOfEquations*leads;i++){rvec[i]=rvec[i]-impactPart2[i];
}
free(impactPart1);
free(impactPart2);
}
@}




\subsection{nxtGuess Definition}
\label{sec:nxtGuess}


@o stackC.h -d
@{
void nxtGuess(@<nxtGuess argument list@>);
void newNxtGuess(@<nxtGuess argument list@>);
@}


@d nxtGuess argument list
@{
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads, unsigned int * capT,
double **fmats,unsigned int  **fmatsj,unsigned int  **fmatsi,
double **smats,unsigned int  **smatsj,unsigned int  **smatsi,
unsigned int *maxNumberHElements,
double * termConstr,unsigned int * termConstrj,unsigned int * termConstri,double * fp,
double * initialX,
double * shockVec,
double * updateDirection @| termConstr fp initialX shockVec
theFunc theDrvFunc capT 
@}

@o stackC.h -d
@{
void chkDrv(@<chkDrv argument list@>);
@}

@d chkDrv argument list
@{
unsigned int n, double * fdrv,unsigned int * fdrvj,unsigned int * fdrvi,
double * fvec,double * delxvec
@}

@o stackC.h -d
@{
void constructFdrv(@<constructFdrv argument list@>);
@}

@d constructFdrv argument list
@{
unsigned int numberOfEquations,unsigned int lags, unsigned int leads,unsigned int pathLength,
double * xvec,double * params,void (* vFunc)(),void (* vFuncDrv)(),
double * termConstr,unsigned int * termConstrj,unsigned int * termConstri,
double * fixedPoint,double * intercept,double * linearizationPoint,unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double * shockVec,
double * fvec,
double * fdrv,unsigned int * fdrvj,unsigned int * fdrvi,unsigned int ihomotopy,
unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo
@}


@d chkDrv definition
@{
#include <math.h>
#define NEGLIGIBLEDOUBLE 1.0e-9
void chkDrv(@<chkDrv argument list@>)
{
unsigned int i;
/*unsigned int aOne=1;*/
double * fvals;
fvals = (double * ) calloc(n,sizeof(double));
#ifdef DEBUG 
printf("chkDrv:beginning\n");
#endif


sparseMatTimesVec(&n,fdrv,fdrvj,fdrvi,delxvec,fvals);
for(i=0;i<n;i++){
if(fabs(fvals[i]-fvec[i])>NEGLIGIBLEDOUBLE){
#ifdef DEBUG 
printf("chkDrv:discrepancy for %d,(%e,%e)\n",i,fvals[i],fvec[i]);
#endif
}
}
free(fvals);
#ifdef DEBUG 
printf("chkDrv:done\n");
#endif

}
#define addOneToFEvals (intOutputInfo[2])++
#define homotopyAlpha (doubleControlParameters+10)
#define addOneToFDrvEvals (intOutputInfo[3])++
#define maxNumberStackElements intControlParameters[5]
void constructFdrv(@<constructFdrv argument list@>)
{
double * deviations;
unsigned int * ignore;/*double  dignore[1]={1.0};*/
unsigned int rowDim;unsigned int * fvecj;unsigned int * fveci;
unsigned int i;unsigned int j;unsigned int soFar;double * zeroShockVec;
ignore = (unsigned int *)calloc(numberOfEquations+1,sizeof(unsigned int));
fvecj = (unsigned int *)calloc(numberOfEquations+1,sizeof(unsigned int));
fveci = (unsigned int *)calloc(numberOfEquations+1,sizeof(unsigned int));
deviations = (double * ) calloc(numberOfEquations*(leads+lags),sizeof(double));
zeroShockVec = (double * ) calloc(numberOfEquations,sizeof(double));
/*identity matrix for lagged values*/
for(i=0;i<numberOfEquations*lags;i++){
fvec[i]=0;
fdrv[i]=1;
fdrvj[i]=i+1;
fdrvi[i]=i+1;}
fdrvi[i]=i+1;
soFar=numberOfEquations*lags;

/*fill in derivative matrices*/
{/*shock only applies for first period*/
addOneToFEvals;
vFunc(xvec,params,shockVec,
fvec+numberOfEquations*lags,fvecj,fveci,homotopyAlpha+ihomotopy,linearizationPoint,exogRows,exogCols,exogenizeQ);
addOneToFDrvEvals;
vFuncDrv(xvec,params,shockVec,
fdrv+soFar,fdrvj+soFar,fdrvi+numberOfEquations*lags,homotopyAlpha+ihomotopy,linearizationPoint,exogRows,exogCols,exogenizeQ);
for(j=0;j<numberOfEquations+1;j++){
fdrvi[numberOfEquations*lags+j]=
fdrvi[numberOfEquations*lags+j]+soFar;
}
for(j=0;j<fdrvi[(lags+1)*numberOfEquations]-
fdrvi[(lags)*numberOfEquations];j++){
fdrvj[j+fdrvi[(lags)*numberOfEquations]-1]=
fdrvj[j+fdrvi[(lags)*numberOfEquations]-1];
}
soFar=fdrvi[numberOfEquations*(lags+1)]-1;
}
for(i=1;i<pathLength;i++){
addOneToFEvals;
vFunc(xvec+i*numberOfEquations,params,zeroShockVec,
fvec+numberOfEquations*lags+i*numberOfEquations,fvecj,fveci,homotopyAlpha+ihomotopy,linearizationPoint,exogRows,exogCols,exogenizeQ);
addOneToFDrvEvals;
vFuncDrv(xvec+i*numberOfEquations,params,zeroShockVec,
fdrv+soFar,fdrvj+soFar,fdrvi+numberOfEquations*lags+(i*numberOfEquations),homotopyAlpha+ihomotopy,linearizationPoint,exogRows,exogCols,exogenizeQ);
for(j=0;j<numberOfEquations+1;j++){
fdrvi[(i*numberOfEquations)+numberOfEquations*lags+j]=
fdrvi[(i*numberOfEquations)+numberOfEquations*lags+j]+soFar;
}
for(j=0;j<fdrvi[(i+lags+1)*numberOfEquations]-
fdrvi[(i+lags)*numberOfEquations];j++){
fdrvj[j+fdrvi[(i+lags)*numberOfEquations]-1]=
fdrvj[j+fdrvi[(i+lags)*numberOfEquations]-1]+i*numberOfEquations;
}
soFar=fdrvi[numberOfEquations*(lags+i+1)]-1;
}
pathNewtAssert(soFar<maxNumberStackElements*pathLength);
/*fill in terminal constraint*/
for(i=0;i<termConstri[numberOfEquations*leads]-termConstri[0];i++){
fdrv[soFar+i]=termConstr[i];
fdrvj[soFar+i]=termConstrj[i]+(numberOfEquations*pathLength);
}
for(i=0;i<numberOfEquations*leads+1;i++){
fdrvi[numberOfEquations*(lags+pathLength)+i]=termConstri[i]+soFar;
}

/*xxxxxxxxx add code for deviations*/
for(i=0;i<numberOfEquations* (lags+ leads);i++){
deviations[i]=xvec[numberOfEquations* pathLength+i]-fixedPoint[i+numberOfEquations];}
rowDim=numberOfEquations*leads;
sparseMatTimesVec(&rowDim,
termConstr,termConstrj,termConstri,deviations,fvec+(numberOfEquations*(lags+pathLength)));
for(i=0;i<numberOfEquations*leads;i++){fvec[numberOfEquations*(lags+pathLength)+i]=fvec[numberOfEquations*(lags+pathLength)+i]-intercept[i];}
free(ignore);
free(fvecj);
free(fveci);
free(deviations);
free(zeroShockVec);
}
@}

@d newNxtGuess argument list
@{
unsigned int * sysDim,
unsigned int * maxNumberHElements,
double * fvec,
double * fdrv,unsigned int * fdrvj,unsigned int * fdrvi,
double * xdel,
unsigned int * ma50bdJob,
unsigned int * ma50bdIq,
double * ma50bdFact,
unsigned int * ma50bdIrnf,
unsigned int * ma50bdIptrl,
unsigned int * ma50bdIptru,
unsigned int * intControlParameters,double * doubleControlParameters
@}

@d newNxtGuess definition
@{
#define sysDimSwitchLevel 30000
#define ma50DropThreshold (1e-8)

void newNxtGuess(@<newNxtGuess argument list@>)
{


unsigned int i;
unsigned int maxElementsEncountered=0;
double * copychkfdrv;unsigned int * copychkfdrvj;unsigned int * copychkfdrvi;
unsigned int * jcn;
double * cntl;
unsigned int * icntl;
unsigned int * ip ;
unsigned int * np;
unsigned int * jfirst;
unsigned int * lenr;
unsigned int * lastr;
unsigned int * nextr;
unsigned int * ifirst;
unsigned int * lenc;
unsigned int * lastc;
unsigned int * nextc;
unsigned int * info;
double * rinfo;
unsigned int *lfact;
double * fact;
unsigned int *irnf;
unsigned int * iptrl;
unsigned int * iptru;
unsigned int *iw;
double * w;
unsigned int nzmax;
unsigned int  * aOne;
unsigned int nonZeroNow;
unsigned int  trans;
@}
@d newNxtGuess definition
@{

copychkfdrv = (double * )calloc(*maxNumberHElements,sizeof(double));
copychkfdrvj = (unsigned int * )calloc(*maxNumberHElements,sizeof(unsigned int));
copychkfdrvi=(unsigned int * ) calloc(*sysDim+1,sizeof(unsigned int));
jcn = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
cntl= (double *)calloc(5,sizeof(double));
icntl= (unsigned int *)calloc(9,sizeof(unsigned int));
ip = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
np = (unsigned int *)calloc(1,sizeof(unsigned int));
jfirst = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
lenr = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
lastr = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
nextr = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
w = (double *)calloc(*sysDim,sizeof(double));
iw = (unsigned int *)calloc(3**sysDim,sizeof(unsigned int));
ifirst = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
lenc = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
lastc = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
nextc = (unsigned int *)calloc(*sysDim,sizeof(unsigned int));
info = (unsigned int *)calloc(7,sizeof(unsigned int));
rinfo = (double *)calloc(1,sizeof(double));
lfact = (unsigned int *)calloc(1,sizeof(unsigned int));
*lfact = ( *maxNumberHElements);/*pessimistic setting for filling*/
aOne = (unsigned int *)calloc(1,sizeof(unsigned int));
*aOne=1;
@}
@d newNxtGuess definition
@{
copmat_(sysDim,fdrv,fdrvj,fdrvi,
copychkfdrv,copychkfdrvj,copychkfdrvi,aOne,aOne);

useMA50ID(cntl,icntl);
cntl[1]=ma50Balance;
cntl[2]=ma50DropEntry;
cntl[3]=ma50DropCol;
icntl[3]=ma50PivotSearch;
/*if(*sysDim>sysDimSwitchLevel){cntl[2]=ma50DropThreshold;}*/
nzmax=*maxNumberHElements;
nonZeroNow=copychkfdrvi[*sysDim]-copychkfdrvi[0];
useMA50AD(sysDim,sysDim,&nonZeroNow,
&nzmax,copychkfdrv,copychkfdrvj,jcn,copychkfdrvi,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
wordybump(info[3]);
pathNewtAssert(info[0]>=0);

#ifdef DEBUG 
printf("\n ma50ad info\n");
for(i=0;i<7;i++)printf(" %d ",info[i]);
printf("\n ma50ad info\n");
#endif

if(*ma50bdJob!=2){
useMA50BD(sysDim,sysDim,&nonZeroNow,ma50bdJob,
fdrv,fdrvj,fdrvi,
cntl,icntl,ip,copychkfdrvi,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump(info[3]);
#ifdef DEBUG 
printf("\n ma50bd info\n");
for(i=0;i<7;i++)printf(" %d ",info[i]);
printf("\n ma50bd info\n");
#endif
pathNewtAssert(info[0]>=0);
if(*ma50bdJob=1)*ma50bdJob=1;
/* if it was 1 promote 
                             to 3(ie conservative alternative)*/
/*unless we're dropping terms with cntl[2]>0*/
if(cntl[2]>0){*ma50bdJob=1;}
} else {
useMA50BD(sysDim,sysDim,&nonZeroNow,ma50bdJob,
fdrv,fdrvj,fdrvi,
cntl,icntl,ip,copychkfdrvi,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump(info[3]);
pathNewtAssert(info[0]>=0);
if(info[0]<-7){/* small pivot values,
reset to 3*/
#ifdef DEBUG 
printf("small pivots, resetting to 3\n");
#endif
*ma50bdJob=1;
/*unless we're dropping terms with cntl[2]>0*/
if(cntl[2]>0){*ma50bdJob=1;}
useMA50BD(sysDim,sysDim,&nonZeroNow,ma50bdJob,
fdrv,fdrvj,fdrvi,
cntl,icntl,ip,copychkfdrvi,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump(info[3]);
pathNewtAssert(info[0]>=0);
}
}
    trans = 1;
useMA50CD(sysDim,sysDim,icntl,copychkfdrvi,np,&trans,
lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
fvec,xdel,
w,info);
*maxNumberHElements=maxElementsEncountered;
pathNewtAssert(info[0]>=0);
@}
@d newNxtGuess definition
@{
free(aOne);
free(jcn );
free(cntl);
free(icntl);
free(ip );
free(np );
free(jfirst );
free(lenr );
free(lastr );
free(nextr );
free(w);
free(iw);
free(ifirst );
free(lenc );
free(lastc );
free(nextc );
free(info );
free(rinfo );
free(lfact);
free(copychkfdrv);free(copychkfdrvj);free(copychkfdrvi);
}
@}


@d nxtGuess definition
@{
void nxtGuess(@<nxtGuess argument list@>)
{
@<nxtGuess variable declarations@>
@<nxtGuess variable initializations@>
@<nxtGuess storage allocations@>
@<nxtGuess initialize lagged C and d matrices@>

for(tNow=0;tNow<*capT;tNow++){
if(tNow< *leads){*ma50bdJob=1;} else {*ma50bdJob=1;}
@<nxtGuess obtain sparse representation and compute next C and d@>
}
if(*leads>0){
@<nxtGuess use terminal constraint@>
}
@<nxtGuess backsolve@>
@<nxtGuess storage deallocations@>

}
@}


@d nxtGuess variable declarations
@{
unsigned int maxElementsEncountered=0;
unsigned int * ma50bdJob;
unsigned int * ma50bdIq;
double * ma50bdFact;
unsigned int * ma50bdIrnf;
unsigned int * ma50bdIptrl;
unsigned int * ma50bdIptru;

unsigned int tNow;
double * deviations;
double **cmats;unsigned int  **cmatsj;unsigned int  **cmatsi;
double **dmats;unsigned int **dmatsj;unsigned int **dmatsi;
double **ymats;unsigned int  **ymatsj;unsigned int  **ymatsi;
double *gmats;unsigned int  *gmatsj;unsigned int  *gmatsi;
unsigned int i,j;
unsigned int *hColumns;
unsigned int *qColumns;
unsigned int *rowDim;
double *fullfvec;
double *fulldfvec;
double *fullXvec;
unsigned int * aOne;
unsigned int * ierr;
@}


@d nxtGuess variable initializations
@{
aOne=(unsigned int *)calloc(1,sizeof(unsigned int));
*aOne=1;
ierr=(unsigned int *)calloc(1,sizeof(unsigned int));
deviations=(double *)calloc(*numberOfEquations*(*lags+*leads),sizeof(double));
hColumns=(unsigned int *)calloc(1,sizeof(unsigned int));
qColumns=(unsigned int *)calloc(1,sizeof(unsigned int));
rowDim=(unsigned int *)calloc(1,sizeof(unsigned int));
*hColumns=*numberOfEquations*(*lags+1+*leads);
*qColumns=*numberOfEquations*(*lags+*leads);
*rowDim=*numberOfEquations* (*leads? *leads:1);

fullfvec=(double *)calloc(*numberOfEquations**leads,sizeof(double));
fulldfvec=(double *)calloc(*numberOfEquations*(*lags+*leads+ 1),sizeof(double));
fullXvec = (double *)calloc(*numberOfEquations *  (*lags+*leads+*capT),sizeof(double));
for(i=0;i<(*lags+*leads+*capT) * *numberOfEquations ;i++){fullXvec[i]=1.0;};
double ** ptrToPtrToDouble = NULL;
unsigned int ** ptrToPtrToInt = NULL;

cmats =(double **)calloc(*capT+(*lags+*leads)+1,sizeof(ptrToPtrToDouble));
cmatsj =(unsigned int **)calloc(*capT+(*lags+*leads)+1,sizeof(ptrToPtrToInt));
cmatsi =(unsigned int **)calloc(*capT+(*lags+*leads)+1,sizeof(ptrToPtrToInt));
dmats =(double **)calloc(*capT+(*lags+*leads)+1,sizeof(ptrToPtrToDouble));
dmatsj =(unsigned int **)calloc(*capT+(*lags+*leads)+1,sizeof(ptrToPtrToInt));
dmatsi =(unsigned int **)calloc(*capT+(*lags+*leads)+1,sizeof(ptrToPtrToInt));
gmats =(double *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(double));
gmatsj =(unsigned int *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(unsigned int));
gmatsi =(unsigned int *)calloc(*numberOfEquations*(*leads?*leads:1)+1,sizeof(unsigned int));
ymats =(double **)calloc(*capT+(*leads+*lags)+1,sizeof(double *));
ymatsj =(unsigned int **)calloc(*capT+(*leads+*lags)+1,sizeof(unsigned int *));
ymatsi =(unsigned int **)calloc(*capT+(*leads+*lags)+1,sizeof(unsigned int *));
for(i=0;i<*capT+(*lags+*leads)+1;i++){
ymats[i] =(double *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(double));
ymatsj[i] =(unsigned int *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(unsigned int));
ymatsi[i] =(unsigned int *)calloc(*numberOfEquations*(*leads?*leads:1)+1,sizeof(unsigned int));
cmats[i] =(double *)calloc(*maxNumberHElements,sizeof(double));
cmatsj[i] =(unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
cmatsi[i] =(unsigned int *)calloc(*numberOfEquations*(*leads?*leads:1)+1,sizeof(unsigned int));
dmats[i] =(double *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(double));
dmatsj[i] =(unsigned int *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(unsigned int));
dmatsi[i] =(unsigned int *)calloc(*numberOfEquations*(*leads?*leads:1)+1,sizeof(unsigned int));
}

@}
@d nxtGuess storage allocations
@{
ma50bdIptru = (unsigned int *)calloc(*numberOfEquations* (*leads?*leads:1),sizeof(unsigned int));
ma50bdIptrl = (unsigned int *)calloc(*numberOfEquations* (*leads?*leads:1),sizeof(unsigned int));
ma50bdIrnf = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
ma50bdFact = (double *)calloc(*maxNumberHElements,sizeof(double));
ma50bdIq = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
ma50bdJob = (unsigned int *)calloc(1,sizeof(unsigned int));
@}
@d nxtGuess initialize lagged C and d matrices
@{
for(i=0;i<*lags;i++){
for(j=0;j<*rowDim;j++){
cmatsi[i][j]=1;
dmatsi[i][j]=1;
cmatsj[i][j]=1;
dmatsj[i][j]=1;
cmats[i][j]=0;
dmats[i][j]=0;
}
cmatsi[i][*rowDim]=1;
dmatsi[i][*rowDim]=1;
}
printf("at init nxtCDmats call\n");
cPrintSparse(*rowDim,(cmats[*lags]),(cmatsj[*lags]),cmatsi[*lags]);

@}

@d nxtGuess obtain sparse representation and compute next C and d
@{

printf("before nxtCDmats call\n");
cPrintSparse(*rowDim,(cmats[*lags]),(cmatsj[*lags]),cmatsi[*lags]);


nxtCDmats(numberOfEquations,lags,leads,
numberOfEquations,maxNumberHElements,
   smats[tNow],smatsj[tNow],smatsi[tNow],
   fmats[tNow],fmatsj[tNow],fmatsi[tNow],
   cmats+*lags+tNow,cmatsj+*lags+tNow,cmatsi+*lags+tNow,
   dmats+*lags+tNow,dmatsj+*lags+tNow,dmatsi+*lags+tNow);

@}
@d nxtGuess use terminal constraint
@{

nxtCDmats(numberOfEquations,lags,leads,
rowDim,maxNumberHElements,
smats[*capT],smatsj[*capT],smatsi[*capT],
fmats[*capT],fmatsj[*capT],fmatsi[*capT],
   cmats+*lags+*capT,cmatsj+*lags+*capT,cmatsi+*lags+*capT,
   dmats+*lags+*capT,dmatsj+*lags+*capT,dmatsi+*lags+*capT);

cmatsi[*lags+*capT+1][0]=cmatsi[*lags+*capT+1][1]=1;
dmatsi[*lags+*capT+1][0]=dmatsi[*lags+*capT+1][1]=1;


@}
@d nxtGuess backsolve
@{
if(*leads>0){
oneStepBack(rowDim,
   ymats+*lags+*capT,ymatsj+*lags+*capT,ymatsi+*lags+*capT,
   cmats+*lags+*capT,cmatsj+*lags+*capT,cmatsi+*lags+*capT,
   dmats+*lags+*capT,dmatsj+*lags+*capT,dmatsi+*lags+*capT);
}
for(i=*capT-1;i>-1;i--){
oneStepBack(numberOfEquations,
   ymats+*lags+i,ymatsj+*lags+i,ymatsi+*lags+i,
   cmats+*lags+i,cmatsj+*lags+i,cmatsi+*lags+i,
   dmats+*lags+i,dmatsj+*lags+i,dmatsi+*lags+i);
}

for(i=0;i<*capT;i++){
csrToDns(numberOfEquations,aOne,ymats[i+*lags],ymatsj[i+*lags],ymatsi[i+*lags],
updateDirection+((*lags + i) * *numberOfEquations),numberOfEquations,ierr);
pathNewtAssert(*ierr == 0);
bump(ymatsi[i+*lags][*numberOfEquations]-ymatsi[i+*lags][0]);
}
if(*leads>0){
csrToDns(rowDim,aOne,ymats[*capT+*lags],ymatsj[*capT+*lags],ymatsi[*capT+*lags],
updateDirection+((*lags + *capT) * *numberOfEquations),rowDim,ierr);
pathNewtAssert(*ierr == 0);
bump(ymatsi[*capT+*lags][*numberOfEquations]-ymatsi[*capT+*lags][0]);
}

@}
@d nxtGuess storage deallocations
@{
free(ma50bdIptru);
free(ma50bdIptrl);
free(ma50bdIrnf);
free(ma50bdFact);
free(ma50bdIq);
free(ma50bdJob);
free(deviations);
free(hColumns);
free(qColumns);
free(rowDim);
free(fullXvec);
free(fullfvec);
free(fulldfvec);
free(ierr);
free(aOne);

for(i=0;i<*capT+(*lags+*leads)+1;i++){
free(cmats[i]);
free(cmatsj[i]);
free(cmatsi[i]);
free(dmats[i]);
free(dmatsj[i]);
free(dmatsi[i]);
free(ymats[i]);
free(ymatsj[i]);
free(ymatsi[i]);
}
free(cmats);
free(cmatsj);
free(cmatsi);
free(dmats);
free(dmatsj);
free(dmatsi);
free(gmats);
free(gmatsj);
free(gmatsi);
free(ymats);
free(ymatsj);
free(ymatsi);

@}

\subsection{nxtFPGuess Definition}
\label{sec:nxtFPGuess}


@o stackC.h -d
@{
void nxtFPGuess(@<nxtFPGuess argument list@>);
@}



@d nxtFPGuess argument list
@{
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
double * fmats,unsigned int * fmatsj,unsigned int * fmatsi,
double * smats,unsigned int * smatsj,unsigned int * smatsi,
unsigned int *maxNumberHElements,
double * initialX,double * updateDirection @| termConstr fp initialX 
theFunc theDrvFunc capT 
@}


@d nxtFPGuess definition
@{
void nxtFPGuess(@<nxtFPGuess argument list@>)
{
@<nxtFPGuess variable declarations@>
@<nxtFPGuess variable initializations@>
@<nxtFPGuess storage allocations@>
@<nxtFPGuess initialize C and d matrices@>
@<nxtFPGuess obtain sparse representation and compute next C and d@>
@<nxtFPGuess use terminal negative identities@>
@<nxtFPGuess backsolve@>
@<nxtFPGuess storage deallocations@>

}
@}


@d nxtFPGuess variable declarations
@{
unsigned int * ma50bdJob;
unsigned int * ma50bdIq;
double * ma50bdFact;
unsigned int * ma50bdIrnf;
unsigned int * ma50bdIptrl;
unsigned int * ma50bdIptru;

unsigned int j;
double **cmats;unsigned int  **cmatsj;unsigned int  **cmatsi;
double **dmats;unsigned int **dmatsj;unsigned int **dmatsi;
double **ymats;unsigned int  **ymatsj;unsigned int  **ymatsi;
double *gmats;unsigned int  *gmatsj;unsigned int  *gmatsi;
unsigned int i;
unsigned int *hColumns;
unsigned int *qColumns;
unsigned int *rowDim;
double *fullfvec;
double *fulldfvec;
/*double *fullXvec;*/
unsigned int * aOne;
unsigned int * ierr;
@}


@d nxtFPGuess variable initializations
@{
aOne=(unsigned int *)calloc(1,sizeof(unsigned int));
*aOne=1;
ierr=(unsigned int *)calloc(1,sizeof(unsigned int));
hColumns=(unsigned int *)calloc(1,sizeof(unsigned int));
qColumns=(unsigned int *)calloc(1,sizeof(unsigned int));
rowDim=(unsigned int *)calloc(1,sizeof(unsigned int));
*hColumns=*numberOfEquations*(*lags+1+*leads);
*qColumns=*numberOfEquations*(*lags+*leads);
*rowDim=*numberOfEquations* *leads;

fullfvec=(double *)calloc(*numberOfEquations**leads,sizeof(double));
fulldfvec=(double *)calloc(*numberOfEquations*(*lags+*leads+ 1),sizeof(double));

double ** ptrToPtrToDouble = NULL;
unsigned int ** ptrToPtrToInt = NULL;

cmats =(double **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToDouble));
cmatsj =(unsigned int **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToInt));
cmatsi =(unsigned int **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToInt));
dmats =(double **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToDouble));
dmatsj =(unsigned int **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToInt));
dmatsi =(unsigned int **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToInt));
gmats =(double *)calloc(*numberOfEquations,sizeof(double));
gmatsj =(unsigned int *)calloc(*numberOfEquations,sizeof(unsigned int));
gmatsi =(unsigned int *)calloc(*numberOfEquations,sizeof(unsigned int));
ymats =(double **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToDouble));
ymatsj =(unsigned int **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToInt));
ymatsi =(unsigned int **)calloc((*lags+*leads)+2,sizeof(ptrToPtrToInt));
for(i=0;i<(*lags+*leads)+2;i++){
ymats[i] =(double *)calloc(*numberOfEquations* *leads,sizeof(double));
ymatsj[i] =(unsigned int *)calloc(*numberOfEquations* *leads,sizeof(unsigned int));
ymatsi[i] =(unsigned int *)calloc(*numberOfEquations* *leads+1,sizeof(unsigned int));
cmats[i] =(double *)calloc(*maxNumberHElements,sizeof(double));
cmatsj[i] =(unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
cmatsi[i] =(unsigned int *)calloc((*numberOfEquations  * *leads) +1,sizeof(unsigned int));
dmats[i] =(double *)calloc(*numberOfEquations* *leads,sizeof(double));
dmatsj[i] =(unsigned int *)calloc(*numberOfEquations* *leads,sizeof(unsigned int));
dmatsi[i] =(unsigned int *)calloc(*numberOfEquations* *leads+1,sizeof(unsigned int));
}

@}
@d nxtFPGuess storage allocations
@{
ma50bdIptru = (unsigned int *)calloc(*numberOfEquations* *leads,sizeof(unsigned int));
ma50bdIptrl = (unsigned int *)calloc(*numberOfEquations* *leads,sizeof(unsigned int));
ma50bdIrnf = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
ma50bdFact = (double *)calloc(*maxNumberHElements,sizeof(double));
ma50bdIq = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
ma50bdJob = (unsigned int *)calloc(1,sizeof(unsigned int));
@}
@d nxtFPGuess initialize C and d matrices
@{
/*construct negative identity matrices*/
for(i=0;i<*lags;i++){
for(j=0;j<*numberOfEquations;j++){
cmats[i][j]=-1;
cmatsi[i][j]=j+1;
cmatsj[i][j]=j+1;
dmatsi[i][j]=1;
}
cmatsi[i][*numberOfEquations]=*numberOfEquations+1;
dmatsi[i][*numberOfEquations]=1;
dmatsj[i][0]=0;
}


@}
@d nxtFPGuess obtain sparse representation and compute next C and d
@{
/*
dnsToCsr(numberOfEquations,aOne,numberOfEquations,fullfvec,
numberOfEquations,
fmats,fmatsj,fmatsi,
ierr);
pathNewtAssert(*ierr == 0);
bump(fmatsi[*numberOfEquations]-fmatsi[0]);
dnsToCsr(numberOfEquations,hColumns,maxNumberHElements,
fulldfvec,
numberOfEquations,
smats,smatsj,smatsi,
ierr);
pathNewtAssert(*ierr == 0);
bump(smatsi[*numberOfEquations]-smatsi[0]);
*/
*ma50bdJob=1;

printf("before nxtCDmats call\n");
cPrintSparse(*rowDim,(cmats[*lags]),(cmatsj[*lags]),cmatsi[*lags]);


nxtCDmats(numberOfEquations,lags,leads,
numberOfEquations,maxNumberHElements,
   smats,smatsj,smatsi,
   fmats,fmatsj,fmatsi,
   cmats+*lags,cmatsj+*lags,cmatsi+*lags,
   dmats+*lags,dmatsj+*lags,dmatsi+*lags);

@}
@d nxtFPGuess use terminal negative identities
@{

/*construct negative identity matrices*/
for(j=0;j<*numberOfEquations* *leads;j++){
smats[2*j]=1;
smatsj[2*j]=(*numberOfEquations * (*lags-1))+j+1;


smats[2*j+1]=-1;
smatsj[2*j+1]=j+1+(*numberOfEquations * (*lags ));

smatsi[j]=2*j+1;
}
smatsi[*numberOfEquations* *leads]=2* (*numberOfEquations * *leads)+1;

for(j=0;j<=*numberOfEquations* *leads;j++)fmatsi[j]=1;
fmatsj[0]=fmatsj[1]=1;fmats[0]=0;

*ma50bdJob=1;
nxtCDmats(numberOfEquations,lags,leads,
rowDim,maxNumberHElements,
smats,smatsj,smatsi,
fmats,fmatsj,fmatsi,
   cmats+*lags+1,cmatsj+*lags+1,cmatsi+*lags+1,
   dmats+*lags+1,dmatsj+*lags+1,dmatsi+*lags+1);


for(j=0;j<=*numberOfEquations;j++){
cmatsi[*lags+*leads+1][j]=1;
dmatsi[*lags+*leads+1][j]=1;
}
@}
@d nxtFPGuess backsolve
@{
/*
oneStepBack(rowDim,
   ymats+*lags+*leads,ymatsj+*lags+*leads,ymatsi+*lags+*leads,
   cmats+*lags+*leads,cmatsj+*lags+*leads,cmatsi+*lags+*leads,
   dmats+*lags+*leads,dmatsj+*lags+*leads,dmatsi+*lags+*leads);
*/
oneStepBack(rowDim,
   ymats+*lags+1,ymatsj+*lags+1,ymatsi+*lags+1,
   cmats+*lags+1,cmatsj+*lags+1,cmatsi+*lags+1,
   dmats+*lags+1,dmatsj+*lags+1,dmatsi+i+1);
for(i=*lags;i>-1;i--){
oneStepBack(numberOfEquations,
   ymats+i,ymatsj+i,ymatsi+i,
   cmats+i,cmatsj+i,cmatsi+i,
   dmats+i,dmatsj+i,dmatsi+i);
}

for(i=0;i<*lags+1;i++){
csrToDns(numberOfEquations,aOne,ymats[i],ymatsj[i],ymatsi[i],
updateDirection+(( i) * *numberOfEquations),numberOfEquations,ierr);
pathNewtAssert(*ierr == 0);


}
csrToDns(rowDim,aOne,ymats[*lags+1],ymatsj[*lags+1],ymatsi[*lags+1],
updateDirection+((*lags+1 ) * *numberOfEquations),rowDim,ierr);
pathNewtAssert(*ierr == 0);




@}
@d nxtFPGuess storage deallocations
@{
free(ma50bdIptru);
free(ma50bdIptrl);
free(ma50bdIrnf);
free(ma50bdFact);
free(ma50bdIq);
free(ma50bdJob);
free(hColumns);
free(qColumns);
free(rowDim);

free(fullfvec);
free(fulldfvec);
free(ierr);
free(aOne);
for(i=0;i<(*lags+*leads)+2;i++){
free(cmats[i]);
free(cmatsj[i]);
free(cmatsi[i]);
free(dmats[i]);
free(dmatsj[i]);
free(dmatsi[i]);
free(ymats[i]);
free(ymatsj[i]);
free(ymatsi[i]);
}
free(cmats);
free(cmatsj);
free(cmatsi);
free(dmats);
free(dmatsj);
free(dmatsi);
free(gmats);
free(gmatsj);
free(gmatsi);
free(ymats);
free(ymatsj);
free(ymatsi);

@}



\subsection{compPathError Definition}
\label{sec:compPathError}


@o stackC.h -d
@{
void compPathError(@<compPathError argument list@>);
@}



@d compPathError argument list
@{
unsigned int * numberOfEquations,unsigned int * lags,unsigned int * leads,
void (* theFunction)(),
double * termConstr,unsigned int * termConstrj,unsigned int * termConstri,double * fp,
double * initialX,
double * shockVec,
unsigned int * capT,
unsigned int * maxNumberHElements,
double ** fmats,unsigned int ** fmatsj,unsigned int **fmatsi,
double ** smats,unsigned int ** smatsj,unsigned int **smatsi
@}


@d compPathError definition
@{
void compPathError(@<compPathError argument list@>)
{
unsigned int maxElementsEncountered=0;
unsigned int * rowDim;
unsigned int * qColumns;
unsigned int tNow;
unsigned int * aOne;
unsigned int * aZero;
unsigned int * ierr;
unsigned int i;
double * deviations;
double * zeroShockVec;
double * fullfvec;
ierr = (unsigned int *) calloc(1,sizeof(unsigned int));
aOne= (unsigned int *) calloc(1,sizeof(unsigned int));
*aOne=1;
aZero= (unsigned int *) calloc(1,sizeof(unsigned int));
*aZero=1;
rowDim= (unsigned int *) calloc(1,sizeof(unsigned int));
qColumns= (unsigned int *) calloc(1,sizeof(unsigned int));
*rowDim=*numberOfEquations* *leads;
*qColumns=*numberOfEquations* (*leads+*lags);
deviations = (double *) calloc(*numberOfEquations* (*lags+*leads+1),sizeof(double));
zeroShockVec = (double *) calloc(*numberOfEquations,sizeof(double));
fullfvec = (double *) calloc(*numberOfEquations,sizeof(double));

{/*use shock for first period only*/
 theFunction(initialX,initialX,
shockVec,
  fmats[0],fmatsj[0],fmatsi[0]);
}
for(tNow=1;tNow<*capT;tNow++) {
 theFunction(initialX+((tNow/*+(*lags-1)*/) * *numberOfEquations),initialX,
zeroShockVec,
  fmats[tNow],fmatsj[tNow],fmatsi[tNow]);
}

copyMatrix(rowDim,aOne,termConstr,termConstrj,termConstri,aOne,
smats[*capT],smatsj[*capT],smatsi[*capT]);
/*
copmat_(rowDim,termConstr,termConstrj,termConstri,
smats[*capT],smatsj[*capT],smatsi[*capT],aOne,aOne);
*/
/*xxxxxxxxx add code for deviations using gmat*/
for(i=0;i<*numberOfEquations* (*lags+ *leads);i++){
deviations[i]=initialX[*numberOfEquations* *capT+i]-fp[i+*numberOfEquations];}
sparseMatTimesVec(rowDim,smats[*capT],smatsj[*capT],smatsi[*capT],deviations,fullfvec);
dnsToCsr(rowDim,aOne,rowDim,fullfvec,
aOne,
fmats[*capT],fmatsj[*capT],fmatsi[*capT],ierr);
pathNewtAssert(*ierr == 0);
bump(fmatsi[*capT][*rowDim]-fmatsi[*capT][0]);
free(zeroShockVec);
}


@}
\subsection{B matrix computation}
\label{sec:bmat}

@d applySparseReducedForm argument list
@{
unsigned int rowDim,unsigned int colDim,
double * initialX, double * fp,
double * bmat,unsigned int * bmatj,unsigned int * bmati,
double * resultX
@}
@d applySparseReducedForm declarations
@{
double * deviations;
unsigned int i;
@}
@d applySparseReducedForm allocations
@{
deviations = (double *) calloc(colDim,sizeof(double));
@}
@d applySparseReducedForm frees
@{
free(deviations);
@}
@d applySparseReducedForm
@{
void applySparseReducedForm(@<applySparseReducedForm argument list@>)
{
@<applySparseReducedForm declarations@>
@<applySparseReducedForm allocations@>
for(i=0;i<colDim;i++){deviations[i]=initialX[i]-fp[i%rowDim];}
sparseMatTimesVec(&rowDim,bmat,bmatj,bmati,deviations,resultX);
for(i=0;i<rowDim;i++){resultX[i]=resultX[i]+fp[i%rowDim];}
@<applySparseReducedForm frees@>
}
@}

@d obtainSparseReducedForm argument list
@{
 unsigned int * maxNumberHElements, 
 unsigned int qrows,  unsigned int qcols, double * qmat,  unsigned int * qmatj, unsigned int * qmati,
 double * bmat,  unsigned int * bmatj,  unsigned int * bmati
@}

@d obtainSparseReducedForm allocations
@{

jcn = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
cntl= (double *)calloc(5,sizeof(double));
icntl= (unsigned int *)calloc(9,sizeof(unsigned int));
ip = (unsigned int *)calloc(qrows,sizeof(unsigned int));
np = (unsigned int *)calloc(1,sizeof(unsigned int));
jfirst = (unsigned int *)calloc(qrows,sizeof(unsigned int));
lenr = (unsigned int *)calloc(qrows,sizeof(unsigned int));
lastr = (unsigned int *)calloc(qrows,sizeof(unsigned int));
nextr = (unsigned int *)calloc(qrows,sizeof(unsigned int));
w = (double *)calloc(qrows,sizeof(double));
iw = (unsigned int *)calloc(3*qrows,sizeof(unsigned int));
ifirst = (unsigned int *)calloc(qrows,sizeof(unsigned int));
lenc = (unsigned int *)calloc(qrows,sizeof(unsigned int));
lastc = (unsigned int *)calloc(qrows,sizeof(unsigned int));
nextc = (unsigned int *)calloc(qrows,sizeof(unsigned int));
info = (unsigned int *)calloc(7,sizeof(unsigned int));
rinfo = (double *)calloc(1,sizeof(double));
/* ma50bd*/

qrmat = (double *) calloc(*maxNumberHElements,sizeof(double));
qrmatj = (unsigned int *) calloc(*maxNumberHElements,sizeof(unsigned int));
qrmati = (unsigned int *) calloc(qrows+1,sizeof(unsigned int);
tb = (double *) calloc(*maxNumberHElements,sizeof(double));
jtb = (unsigned int *) calloc(*maxNumberHElements,sizeof(unsigned int));
itb = (unsigned int *) calloc(qcols-qrows+1,sizeof(unsigned int));
b = (double *) calloc(*maxNumberHElements,sizeof(double));
jb = (unsigned int *) calloc(*maxNumberHElements,sizeof(unsigned int *));
ib = (unsigned int *) calloc(qcols-qrows+1,sizeof(unsigned int *));


lfact =(unsigned int *)calloc(1,sizeof(unsigned int));
*lfact = ( *maxNumberHElements);/*pessimistic setting for filling*/
fact = (double *)calloc(*lfact,sizeof(double));
irnf = (unsigned int *)calloc(*lfact,sizeof(unsigned int));
iptrl = (unsigned int *)calloc(qrows,sizeof(unsigned int));
iptru = (unsigned int *)calloc(qrows,sizeof(unsigned int));
x = (double *)calloc(  qcols,sizeof(double));
nsSumC = (double *)calloc(qrows ,sizeof(double));
@}
@d obtainSparseReducedForm frees
@{
free(w);
free(iw);
free(b);
free(jb);
free(ib);
free(tb);
free(jtb);
free(itb);
free(jcn );
free(cntl);
free(icntl);
free(ip );
free(np );
free(jfirst );
free(lenr );
free(lastr );
free(nextr );
free(ifirst );
free(lenc );
free(lastc );
free(nextc );
free(info );
free(rinfo );
free(/* ma50bd*/

qrmat );
free(qrmatj );
free(qrmati );
free(lfact );
free(fact );
free(irnf );
free(iptrl );
free(iptru );
free(x );
free(nsSumC );
@}


@d obtainSparseReducedForm declarations
@{
void * calloc(unsigned amt,unsigned int size);
double * nsSumC;unsigned int ierr;double * x;
unsigned int nzmaxLeft;double aSmallDouble;double nsSumD;
unsigned int cmatsExtent;unsigned int i;unsigned int cColumns;
double *b;unsigned int *jb,*ib;
double *tb;unsigned int *jtb,*itb;
unsigned int  trans;
 double * qrmat; unsigned int * qrmatj; unsigned int * qrmati;
unsigned int *iw;double * w;
unsigned int  aOne; unsigned int  firstColumn;unsigned int  lastColumn;
unsigned int  nr;unsigned int  nc;unsigned int nonZeroNow;unsigned int nzmax;
unsigned int * jcn;
double * cntl;
unsigned int * icntl;
unsigned int * ip ;
unsigned int * np;
unsigned int * jfirst;
unsigned int * lenr;
unsigned int * lastr;
unsigned int * nextr;
unsigned int * ifirst;
unsigned int * lenc;
unsigned int * lastc;
unsigned int * nextc;
unsigned int * info;
double * rinfo;
unsigned int *lfact;
double * fact;
unsigned int *irnf;
unsigned int * iptrl;
unsigned int * iptru;
@}

@o reducedForm.c -d
@{
#include <float.h>
@<obtainSparseReducedForm@>
@<applySparseReducedForm@>
@}


@d obtainSparseReducedForm
@{
void obtainSparseReducedForm(
@<obtainSparseReducedForm argument list@>
)
{
@<obtainSparseReducedForm declarations@>



@<obtainSparseReducedForm allocations@>


/*solve relation Qr xr = Ql xl and change sign later note xl are just
elements of identity matrix so that  solving Qr xr = Ql will give us
Bmatrix but with wrong sign*/
@<bmat compute pivot sequence@>
@<bmat factorize matrix@>
@<bmat use factorization@>
@<obtainSparseReducedForm frees@>
}

@}



@d bmat compute pivot sequence
@{
/*still using CSR consequently doing everything to the 
transpose*/
/*note ma50ad modifies its A argument*/


firstColumn=(qcols-qrows+1);
lastColumn=qcols;
aOne=1;
extractSubmatrix(&qrows,&aOne,&aOne,&qrows,
&firstColumn,&lastColumn,
qmat,qmatj,qmati,&nr,&nc,
qrmat,qrmatj,qrmati);

nonZeroNow=qrmati[qrows]-qrmati[0];

useMA50ID(cntl,icntl);
nzmax=*maxNumberHElements;
useMA50AD(&qrows,&qrows,&nonZeroNow,
&nzmax,qrmat,qrmatj,jcn,qrmati,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
/* restore odd since ad is destructive*/
extractSubmatrix(&qrows,&aOne,&aOne,&qrows,
&firstColumn,&lastColumn,
qmat,qmatj,qmati,&nr,&nc,
qrmat,qrmatj,jcn);

@}


@d bmat factorize matrix
@{
useMA50BD(&qrows,&qrows,&nonZeroNow,&aOne,
qrmat,qrmatj,jcn,
cntl,icntl,ip,qrmati,np,lfact,fact,irnf,iptrl,iptru,
w,iw,info,rinfo);
@}

MA50CD applies the factoriation to solve
\begin{gather*}
  A^T x = b
\end{gather*}




@d bmat use factorization

@{
    trans = 1;

/*expand sum of c's use transpose since c colum major order */

itb[0]=1;cmatsExtent=0;
cColumns=qcols-qrows;
for(i=0;i<cColumns;i++){


lastColumn = firstColumn=(1+i);

extractSubmatrix(&qrows,&aOne,
&aOne,&qrows,&firstColumn,&lastColumn,
qmat,qmatj,qmati,
&nr,&nc,
b,jb,ib);

csrToDns(&qrows,&aOne,b,jb,ib,
nsSumC,&qrows,&ierr);

useMA50CD(&qrows,&qrows,icntl,qrmati,np,&trans,
lfact,fact,irnf,iptrl,iptru,
nsSumC,x,
w,info);

nzmaxLeft= nzmax-cmatsExtent-1;
dnsToCsr(&aOne,&qrows,&nzmaxLeft,x,
&aOne,tb+(itb[i]-1),jtb+(itb[i]-1),itb+i,
&ierr);
itb[i+1]=itb[i+1]+cmatsExtent;
itb[i]=itb[i]+cmatsExtent;
cmatsExtent=itb[i+1]-1;
}
aSmallDouble=DBL_EPSILON;
dropSmallElements(&cColumns,&aOne,&aSmallDouble,nzmax,tb,jtb,itb,tb,jtb,itb,&ierr);
csrToCsc(&cColumns,&aOne,&aOne,tb,jtb,itb,bmat,bmatj,bmati);
/*change sign*/
for(i=0;i<bmati[qrows]-bmati[0];i++)bmat[i]=(-1)*bmat[i];
@}

\section{Makefile}
\label{sec:makefile}


@o makeNotStack -t
@{
@<makefile definitions@>
@<makefile generic dependencies@>
@<makefile specific dependencies@>
@}

@d makefile definitions
@{
#don't forget unsetenv TMPDIR
SPARSELIB = -L$(HOME)/dataHome/sparse/SPARSKIT2 -lskit
GEF = $(HOME)/dataHome/gef/GEF-4.0.5

CFLAGS =  -g -I$(GEF)/dist/usr/include/
NUWEBFLAGS = -n
LINTFLAGS = -b -c  -h 
LINKFLAGS =  -g -L/opt/nag/hsml/double/source -lharwell \
        -L/msu/res2/m1gsa00/dataHome/sparse/SPARSKIT2 -lskit\
        /msu/res2/m1gsa00/lapack/LAPACK/blas_os5.a \
        -lc -ldl -lm 

WEBSOURCE =  stackC.w
CSOURCE   = $(WEBSOURCE:.w=.c)
TEXSOURCE   = $(WEBSOURCE:.w=.tex)
HTMLSOURCE   = $(WEBSOURCE:.hw=.tex)
atIFiles    = 
AIMLIB  = -L /msu/res2/m1gsa00/aim/summer98/aimCCode -lfaim -lsrrit lapack_os5.a blas_os5.a

OTHERSOURCE =  miscLatexPkg.tex

SOURCE = $(CSOURCE)  $(OTHERSOURCE)

OBJECT = $(SOURCE:.c=.o)
LINTFILE = $(CSOURCE:.c=.ln)

@}
@d makefile generic dependencies
@{

.SUFFIXES:  .tex .dvi .w .print .c .hw .html .ps

.hw.html:
    nuweb $*.hw
    latex2html -split 0 -no_reuse -no_navigation $*.tex

.w.hw:
    cp $*.w $*.hw

.w.tex:
    nuweb $(NUWEBFLAGS) $*
.dvi.ps:
    dvips -o $*.ps $*

.tex.dvi:
    latex $*
#   bibtex $*
    nuweb  $(NUWEBFLAGS) $*
    latex $*
    latex $*


.w.dvi:
    nuweb  $(NUWEBFLAGS) $*
    latex $*
#   bibtex $*
    nuweb  $(NUWEBFLAGS) $*
    latex $*
    latex $*

.w.c:
    nuweb -t $*

.w.o:
    make $*.c
    make $*.o

.c.o:
    gcc $(CFLAGS) -c $*.c

.dvi.print:
    $(PLPR) $*

# Dependency rules

@}
@d makefile specific dependencies
@{

stackC.w:   $(atIFiles) 
    touch stackC.w


lint:   $(LINTFILE)

printall:   

testStack.c:    stackC.w
	nuweb $(NUWEBFLAGS) -t stackC.w

testStack:   testStack.o    stackC.o testStackModel.o
	f77 -o testStack -g  testStack.o stackC.o testStackModel.o \
	$(SPARSELIB) $(LINKFLAGS) $(AIMLIB)


@}

\section{Numerical Recipes Modifications of FPnewt}

\subsection{array allocation program}
\label{sec:allocarray}

@o stackC.c -d
@{
#include "stochProto.h"
void allocMa50(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int pathLength,unsigned int maxElements,
unsigned int **ma50bdIptru,
unsigned int **ma50bdIptrl,
unsigned int **ma50bdIrnf,
double **ma50bdFact,
unsigned int **ma50bdIq,
unsigned int **ma50bdJob)
{
unsigned int sysDim;
sysDim= numberOfEquations*(lags+pathLength+leads);
*ma50bdIptru = (unsigned int *)calloc(sysDim,sizeof(unsigned int));
*ma50bdIptrl = (unsigned int *)calloc(sysDim,sizeof(unsigned int));
*ma50bdIrnf = (unsigned int *)calloc(maxElements,sizeof(unsigned int));
*ma50bdFact = (double *)calloc(maxElements,sizeof(double));
*ma50bdIq = (unsigned int *)calloc(sysDim,sizeof(unsigned int));
*ma50bdJob = (unsigned int *)calloc(1,sizeof(unsigned int));
}

void allocFPNewt(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int pathLength,unsigned int maxElements,
double ** genericFP,
double ** genericIntercept,
double***fmats,unsigned int***fmatsj,unsigned int***fmatsi,
double***smats,unsigned int***smatsj,unsigned int***smatsi)
{unsigned int i;
if(pathLength<1)pathLength=1;
*genericFP=(double *) calloc(numberOfEquations*(lags+leads+1+pathLength),
sizeof(double));
*genericIntercept=(double *) calloc(numberOfEquations*(leads),
sizeof(double));
*fmats =(double **)calloc((pathLength)+lags+1,sizeof(double *));
*fmatsj =(unsigned int **)calloc((pathLength)+lags+1,sizeof(unsigned int *));
*fmatsi =(unsigned int **)calloc((pathLength)+lags+1,sizeof(unsigned int *));
*smats =(double **)calloc((pathLength)+lags+1,sizeof(double *));
*smatsj =(unsigned int **)calloc((pathLength)+lags+1,sizeof(unsigned int *));
*smatsi =(unsigned int **)calloc((pathLength)+lags+1,sizeof(unsigned int *));
for(i=0;i<(pathLength)+lags+1;i++){
(*fmats)[i] =(double *)calloc(numberOfEquations*(lags+leads),sizeof(double));
(*fmatsj)[i] =(unsigned int *)calloc(numberOfEquations*(lags+leads),sizeof(unsigned int));
(*fmatsi)[i] =(unsigned int *)calloc(
     numberOfEquations*(lags+leads)+1,sizeof(unsigned int));
(*smats)[i] =(double *)calloc(maxElements,sizeof(double));
(*smatsj)[i] =(unsigned int *)calloc(maxElements,sizeof(unsigned int));
(*smatsi)[i] =(unsigned int *)calloc(
     numberOfEquations*(lags+leads)+1,sizeof(unsigned int));
}
}
void freeMa50(
unsigned int **ma50bdIptru,
unsigned int **ma50bdIptrl,
unsigned int **ma50bdIrnf,
double **ma50bdFact,
unsigned int **ma50bdIq,
unsigned int **ma50bdJob)
{
free(*ma50bdIptru);
free(*ma50bdIptrl);
free(*ma50bdIrnf);
free(*ma50bdFact);
free(*ma50bdIq);
free(*ma50bdJob);
}
void freeFPNewt(unsigned int lags, unsigned int pathLength,
double ** genericFP,
double ** genericIntercept,
double***fmats,unsigned int***fmatsj,unsigned int***fmatsi,
double***smats,unsigned int***smatsj,unsigned int***smatsi)
{unsigned int i;
free(*genericFP);
free(*genericIntercept);
for(i=0;i<(pathLength)+lags+1;i++){
free((*fmats)[i]);
free((*fmatsj)[i]);
free((*fmatsi)[i]);
free((*smats)[i]);
free((*smatsj)[i]);
free((*smatsi)[i]);
}
free(*fmats);
free(*fmatsj);
free(*fmatsi);
free(*smats);
free(*smatsj);
free(*smatsi);
}
void allocAltComputeAsymptoticQ(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int maxElements,double**AMqMatrix,unsigned int**AMqMatrixj,unsigned int**AMqMatrixi,
double** rootr,double**rooti)
{
*AMqMatrix=(double *)
   calloc(maxElements,sizeof(double));
*AMqMatrixj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*AMqMatrixi=(unsigned int *)
   calloc((numberOfEquations*(leads+lags)+1),
        sizeof(unsigned int));
*rootr=(double *) calloc((numberOfEquations)*((lags)+(leads)),sizeof(double));
*rooti=(double *) calloc((numberOfEquations)*((lags)+(leads)),sizeof(double));

}


void freeAltComputeAsymptoticQ(
double**AMqMatrix,unsigned int**AMqMatrixj,unsigned int**AMqMatrixi,
double**rootr,double**rooti)
{
free(*AMqMatrix);
free(*AMqMatrixj);
free(*AMqMatrixi);
free(*rootr);
free(*rooti);
}



void allocPhiF(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int numberExogenous,
unsigned int maxElements,
double**psiMatrix,unsigned int**psiMatrixj,unsigned int**psiMatrixi,
double**upsilonMatrix,unsigned int**upsilonMatrixj,unsigned int**upsilonMatrixi,
double**phiMatrix,unsigned int**phiMatrixj,unsigned int**phiMatrixi,
double**fMatrix,unsigned int**fMatrixj,unsigned int**fMatrixi,
double**vartheta,unsigned int**varthetaj,unsigned int**varthetai,
double**impact,unsigned int**impactj,unsigned int**impacti
)
{
*psiMatrix=(double *)
   calloc(maxElements,sizeof(double));
*psiMatrixj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*psiMatrixi=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));

*upsilonMatrix=(double *)
   calloc(maxElements,sizeof(double));
*upsilonMatrixj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*upsilonMatrixi=(unsigned int *)
   calloc((numberExogenous+1),
        sizeof(unsigned int));

*phiMatrix=(double *)
   calloc(maxElements,sizeof(double));
*phiMatrixj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*phiMatrixi=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));

*fMatrix=(double *)
   calloc(maxElements,sizeof(double));
*fMatrixj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*fMatrixi=(unsigned int *)
   calloc((numberOfEquations*(leads+lags)+1),
        sizeof(unsigned int));
*vartheta=(double *)
   calloc(maxElements,sizeof(double));
*varthetaj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*varthetai=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));
*impact=(double *)
   calloc(maxElements,sizeof(double));
*impactj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*impacti=(unsigned int *)
   calloc(((1+leads)*numberOfEquations+1),
        sizeof(unsigned int));
}

void freePhiF(
double**psiMatrix,unsigned int**psiMatrixj,unsigned int**psiMatrixi,
double**upsilonMatrix,unsigned int**upsilonMatrixj,unsigned int**upsilonMatrixi,
double**phiMatrix,unsigned int**phiMatrixj,unsigned int**phiMatrixi,
double**fMatrix,unsigned int**fMatrixj,unsigned int**fMatrixi,
double**vartheta,unsigned int**varthetaj,unsigned int**varthetai,
double**impact,unsigned int**impactj,unsigned int**impacti
)
{
free(*psiMatrix);
free(*psiMatrixj);
free(*psiMatrixi);
free(*upsilonMatrix);
free(*upsilonMatrixj);
free(*upsilonMatrixi);
free(*phiMatrix);
free(*phiMatrixj);
free(*phiMatrixi);
free(*fMatrix);
free(*fMatrixj);
free(*fMatrixi);
free(*vartheta);
free(*varthetaj);
free(*varthetai);
free(*impact);
free(*impactj);
free(*impacti);
}

 
void allocLinearTerminator(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int numberExogenous,
unsigned int maxElements,
double**upsilonMatrix,unsigned int**upsilonMatrixj,unsigned int**upsilonMatrixi,
double**hMat,unsigned int**hMatj,unsigned int**hMati,
double**hzMat,unsigned int**hzMatj,unsigned int**hzMati,
double**cstar,unsigned int**cstarj,unsigned int**cstari,
double**AMqMatrix,unsigned int**AMqMatrixj,unsigned int**AMqMatrixi,
double** rootr,double**rooti,
double**bMat,unsigned int**bMatj,unsigned int**bMati,
double**phiInvMat,unsigned int**phiInvMatj,unsigned int**phiInvMati,
double**fmat,unsigned int**fmatj,unsigned int**fmati,
double**varthetaZstar,unsigned int**varthetaZstarj,unsigned int**varthetaZstari,
double**impact,unsigned int**impactj,unsigned int**impacti,
double**varthetaC,unsigned int**varthetaCj,unsigned int**varthetaCi,
double**selectZmat,unsigned int**selectZmatj,unsigned int**selectZmati
)
{
*upsilonMatrix=(double *)
   calloc(maxElements,sizeof(double));
*upsilonMatrixj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*upsilonMatrixi=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));
*hMat=(double *)
   calloc(maxElements,sizeof(double));
*hMatj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*hMati=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));
*hzMat=(double *)
   calloc(maxElements,sizeof(double));
*hzMatj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*hzMati=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));
*bMat=(double *)
   calloc(maxElements,sizeof(double));
*bMatj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*bMati=(unsigned int *)
   calloc((numberOfEquations*leads+1),
        sizeof(unsigned int));
*phiInvMat=(double *)
   calloc(maxElements,sizeof(double));
*phiInvMatj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*phiInvMati=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));
*fmat=(double *)
   calloc(maxElements,sizeof(double));
*fmatj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*fmati=(unsigned int *)
   calloc((numberOfEquations*leads+1),
        sizeof(unsigned int));
*impact=(double *)
   calloc(maxElements,sizeof(double));
*impactj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*impacti=(unsigned int *)
   calloc((numberOfEquations*(leads+1)+1),
        sizeof(unsigned int));
*varthetaC=(double *)
   calloc(maxElements,sizeof(double));
*varthetaCj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*varthetaCi=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));
*selectZmat=(double *)
   calloc(maxElements,sizeof(double));
*selectZmatj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*selectZmati=(unsigned int *)
   calloc((numberExogenous+1),
        sizeof(unsigned int));
*varthetaZstar=(double *)
   calloc(maxElements,sizeof(double));
*varthetaZstarj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*varthetaZstari=(unsigned int *)
   calloc((numberOfEquations+1),
        sizeof(unsigned int));
*AMqMatrix=(double *)
   calloc(maxElements,sizeof(double));
*AMqMatrixj=(unsigned int *)
   calloc(maxElements,sizeof(unsigned int));
*AMqMatrixi=(unsigned int *)
   calloc((numberOfEquations*(leads+lags)+1),
        sizeof(unsigned int));
*cstar=(double *) calloc((numberOfEquations),sizeof(double));
*cstarj=(unsigned int *) calloc((numberOfEquations),sizeof(unsigned int));
*cstari=(unsigned int *) calloc((numberOfEquations)+1,sizeof(unsigned int));
*rootr=(double *) calloc((numberOfEquations)*((lags)+(leads)),sizeof(double));
*rooti=(double *) calloc((numberOfEquations)*((lags)+(leads)),sizeof(double));
}

void freeLinearTerminator(
double**upsilonMatrix,unsigned int**upsilonMatrixj,unsigned int**upsilonMatrixi,
double**hMat,unsigned int**hMatj,unsigned int**hMati,
double**hzMat,unsigned int**hzMatj,unsigned int**hzMati,
double**cstar,unsigned int**cstarj,unsigned int**cstari,
double**AMqMatrix,unsigned int**AMqMatrixj,unsigned int**AMqMatrixi,
double** rootr,double**rooti,
double**bMat,unsigned int**bMatj,unsigned int**bMati,
double**phiInvMat,unsigned int**phiInvMatj,unsigned int**phiInvMati,
double**fmat,unsigned int**fmatj,unsigned int**fmati,
double**varthetaZstar,unsigned int**varthetaZstarj,unsigned int**varthetaZstari,
double**impact,unsigned int**impactj,unsigned int**impacti,
double**varthetaC,unsigned int**varthetaCj,unsigned int**varthetaCi,
double**selectZmat,unsigned int**selectZmatj,unsigned int**selectZmati
)
{
free(*upsilonMatrix);
free(*upsilonMatrixj);
free(*upsilonMatrixi);
free(*hMat);
free(*hMatj);
free(*hMati);
free(*hzMat);
free(*hzMatj);
free(*hzMati);
free(*bMat);
free(*bMatj);
free(*bMati);
free(*phiInvMat);
free(*phiInvMatj);
free(*phiInvMati);
free(*fmat);
free(*fmatj);
free(*fmati);
free(*impact);
free(*impactj);
free(*impacti);
free(*varthetaC);
free(*varthetaCj);
free(*varthetaCi);
free(*selectZmat);
free(*selectZmatj);
free(*selectZmati);
free(*varthetaZstar);
free(*varthetaZstarj);
free(*varthetaZstari);
free(*AMqMatrix);
free(*AMqMatrixj);
free(*AMqMatrixi);
free(*cstar);
free(*cstarj);
free(*cstari);
free(*rootr);
free(*rooti);
}


void allocPathNewt(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int pathLength,unsigned int replications,unsigned int stochasticPathLength,
double**genericPath,
double**genericZeroPath,
double**genericEasyPath,
double**genericTargetPath
)
{
*genericPath=(double *)calloc(
    replications*
    numberOfEquations*(lags+leads+pathLength+stochasticPathLength),
    sizeof(double));
*genericZeroPath=(double *)calloc(
    replications*
    numberOfEquations*(lags+leads+pathLength+stochasticPathLength),
    sizeof(double));
*genericEasyPath=(double *)calloc(
    replications*
    numberOfEquations*(lags+leads+pathLength+stochasticPathLength+1),
    sizeof(double));
*genericTargetPath=(double *)calloc(
    replications*
    numberOfEquations*(lags+leads+pathLength+stochasticPathLength),
    sizeof(double));
}
void freePathNewt(double ** genericPath)
{
free(*genericPath);
}
void allocShockVec(unsigned int numberOfEquations,double**shockVec)
{
*shockVec=(double *)calloc(
    numberOfEquations,sizeof(double));
}
void freeShockVec(double ** shockVec)
{
free(*shockVec);
}
void allocShocksData(unsigned int numberOfEquations,unsigned int numberOfShocks,unsigned int numberOfData,
double**shockVec,double ** dataVec,double ** zeroShockVec)
{
unsigned int i;
*shockVec=(double *)calloc(
    numberOfShocks*numberOfEquations,sizeof(double));
*dataVec=(double *)calloc(
    numberOfData*numberOfEquations,sizeof(double));
*zeroShockVec=(double *)calloc(
    numberOfEquations,sizeof(double));
for(i=0;i<numberOfEquations;i++){(*zeroShockVec)[i]=0.0;}
}
void freeShocksData(double ** shockVec,double ** dataVec,
double ** zeroShockVec)
{
free(*shockVec);
free(*dataVec);
free(*zeroShockVec);
}


@}

\subsection{FPnewt.c}
\label{sec:FPnewt.c}

Numerical Recipes has a globally convergent
routine for find the root of a system of equations\cite{press92}
@d define assert bump
@{

#define wordybump(potentialMaxValue) \
   if(potentialMaxValue>maxElementsEncountered) \
   maxElementsEncountered=potentialMaxValue;\
printf("bump stuff(%d,%d) at line %d",potentialMaxValue,maxElementsEncountered,\
__LINE__);


#define bump(potentialMaxValue) \
   if(potentialMaxValue>maxElementsEncountered) \
   maxElementsEncountered=potentialMaxValue;

#include <signal.h>
#define pathNewtAssert(expression)  \
  if(!(expression))\
		   __pathNewtAssert (expression, __FILE__, __LINE__)

#define __pathNewtAssert(expression, file, lineno)  \
  {printf("pathNewtAssert: processid=%ld\n",(long)getpid());\
   printf ("%s:%u: failed assertion\n", file, lineno);\
   	kill(getpid(),SIGUSR2);}


@}

@o stackC.h -d
@{
void FPnewt(@<FPNewt argument list@>);
@}


@d FPNewt argument list
@{
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
void (* func)(),void (* dfunc)(),double * params,
double x[],
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
unsigned int *check
@}


@o myNewt.c -d 
@{
@<FPnewt defines@>
void FPnewt(@<FPNewt argument list@>)
{
@<FPnewt declarations@>
@<evaluate func and find max element@>
for (its=1;its<=MAXITS;its++) {
@<get newton update@>
@<evaluate func update norm@>
@<check for convergence@>
	}
	printf("MAXITS exceeded in FPnewt");FREERETURN
}
@}
@d FPnewt defines
@{
#include <stdio.h>
#include <math.h>
#include "useSparseAMA.h"
#include "stackC.h"
#define NRANSI
/*#include "./nrutil.h"*/

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

#define MAXITS 200
#define TOLF 1.0e-12
#define TOLMIN 1.0e-6
#define TOLX 1.0e-10
#define STPMX 100.0
#define GAMMA 1.0

unsigned int nn;
double *fvec;
#define FREERETURN {free(fvec);free(xold);free(shockVec);\
	free(p);free(g);free(aOne);free(ierr);\
	free(indx);return;}
#define PFREERETURN {free(fvec);free(xold);free(xoldls);free(xdel);\
free(deviations);free(fullfvec);\
	free(p);free(g);free(aOne);free(ierr);\
	free(indx);free(rowDim);free(qColumns);return;}
@}

@d FPnewt declarations
@{
    unsigned int n;/*double * xdel;*/
/*	double fmin(double x[]);*/
	void lnsrch(unsigned int n,unsigned int np,unsigned int reps,
    double xold[], double   * fold, double g[], double p[],double * params,
		 double * shockVec,double * f, double stpmax, unsigned int *check, 
         void (*func)(double*,double*,double*,double*,unsigned int*,unsigned int*),double *x);
	void lubksb(double **a, unsigned int n, unsigned int *indx, double b[]);
	void ludcmp(double **a, unsigned int n, unsigned int *indx, double *d);
	unsigned int i,its/*,j*/,*indx,*aOne/*,*ndns*/,*ierr;
	double /*d,*/den,f,fold,stpmax,sum,temp,test,*g,*p,*xold,normSum;
    double * shockVec;
    n=*numberOfEquations*(*lags+*leads+1);
	aOne=(unsigned int *)calloc(1,sizeof(unsigned int));
    *aOne=1;
	ierr=(unsigned int *)calloc(1,sizeof(unsigned int));
	indx=(unsigned int *)calloc(n+1,sizeof(unsigned int));
	g=(double *)calloc(n+1,sizeof(double));
	p=(double *)calloc(n+1,sizeof(double));
	xold=(double *)calloc(n+1,sizeof(double));
	fvec=(double *)calloc(n+1,sizeof(double));
	shockVec=(double *)calloc(n+1,sizeof(double));
    for(i=0;i<n+1;i++)shockVec[i]=0.0;
	nn=n;
@}


@d evaluate func and find max element
@{


/* modification begin */
      func(x,params,shockVec,fmats[0],fmatsj[0],fmatsi[0]);
        for (normSum=0.0,i=0;i<=fmatsi[0][*numberOfEquations]-fmatsi[0][0];i++) 
            normSum += SQR(fmats[0][i]);
        f= 0.5 * normSum * (*lags+*leads+1);
        csrToDns(numberOfEquations,
        aOne,fmats[0],fmatsj[0],fmatsi[0],fvec+(*numberOfEquations * *lags),numberOfEquations,ierr);
/* modification end */
	test=0.0;
	for (i=0;i<*numberOfEquations;i++)
		if (fabs(fvec[(*numberOfEquations**lags)+i]) > test) 
        test=fabs(fvec[(*numberOfEquations**lags)+i]);
dfunc(x,params,smats[0],smatsj[0],smatsi[0]);
	if (test<0.01*TOLF) FREERETURN
	for (sum=0.0,i=0;i<*numberOfEquations;i++) sum += SQR(x[i]);
	stpmax=(*lags+*leads+1)*STPMX*FMAX(sqrt(sum),(double)n);
@}

@d get newton update
@{
dfunc(x,params,smats[0],smatsj[0],smatsi[0]);printf("delete this in stackC.w");
		for (i=0;i<n;i++) xold[i]=x[i];
		fold=f;
		/*modification begin*/
nxtFPGuess(numberOfEquations,lags,leads,
fmats[0],fmatsj[0],fmatsi[0],
smats[0],smatsj[0],smatsi[0],
maxNumberElements,x,p);


/*
printf("here is check=%d, and f=%f and x before\n",*check,f);
cPrintMatrixNonZero(n,1,x);
*/


		lnsrch(n,*numberOfEquations,1,
        xold,&fold,g,p,params,shockVec,&f,stpmax,check,func,x);

/*
printf("here is check=%d, and f=%f and x after\n",*check,f);
cPrintMatrixNonZero(n,1,x);
*/
/*
for(i=0;i<*numberOfEquations*(*lags+*leads+1);i++)x[i]=xold[i];
for(i=0;i<n;i++)x[i]=x[i]-p[i];
*/


@}

@d evaluate func update norm
@{

      func(x,params,shockVec,fmats[0],fmatsj[0],fmatsi[0]);
        for (normSum=0.0,i=0;i<=fmatsi[0][*numberOfEquations]+fmatsi[0][0];i++) 
            normSum += SQR(fmats[0][i]);
        f= 0.5 * normSum * (*lags+*leads+1);
        csrToDns(numberOfEquations,
        aOne,fmats[0],fmatsj[0],fmatsi[0],fvec+(*numberOfEquations**lags),
        numberOfEquations,ierr);

        sparseMatTimesVec(numberOfEquations,
        smats[0],smatsj[0],smatsi[0],fvec,
        g+(*numberOfEquations * *lags));
        for(i=0;i>*numberOfEquations;i++){
          g[(*numberOfEquations* (*lags-1))+i]=fvec[(*numberOfEquations* *lags)+i];
          g[(*numberOfEquations* (*lags+1))+i]= -fvec[(*numberOfEquations* *lags)+i];
          }


@}

@d check for convergence
@{
		test=0.0;
		for (i=0;i<*numberOfEquations;i++)
			if (fabs(fvec[(*numberOfEquations**lags)+i]) > test) test=fabs(fvec[(*numberOfEquations**lags)+i]);
		if (test < TOLF) {
			*check=0;
			FREERETURN
		}
		if (*check) {
			test=0.0;
			den=FMAX(f,0.5*n);
			for (i=0;i<n;i++) {
				temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TOLMIN ? 0 : 1);
			FREERETURN
		}
		test=0.0;
		for (i=0;i<n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) FREERETURN
@}


\subsection{pathNewt}
\label{sec:pathNewt}




@o myNewt.c -d 
@{
void pathNewt(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,double * shockVec,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
 unsigned int * maxNumberElements,double * qMat, unsigned int * qMatj, unsigned int * qMati,
double * fixedPoint,
double x[],
 unsigned int *check)
{
@<pathNewt declarations@>
@<pathNewt initializations@>
@<q terminal constraint computation@>
@<pathNewt check for convergence @>
	for (its=1;its<=MAXITS;its++) {
@<pathNewt update path@>
@<q terminal constraint computation@>
@<pathNewt check for convergence@>
		for (i=0;i<n;i++) xold[i]=x[i];
	}
 printf("MAXITS exceeded in pathNewt");PFREERETURN
}
@}


@d q terminal constraint computation
@{
/* modification begin */

@<compute fmats for path@>

@<apply q terminal constraint@>
/*modifications begin */
for(tNow=0;tNow<*pathLength;tNow++) {
      fdjac(x+(tNow* *numberOfEquations),params,
      smats[tNow],smatsj[tNow],smatsi[tNow]);
}
/*modifications end */


@}





@d compute fmats for path
@{
normSum=0.0;
for(tNow=0;tNow<*pathLength;tNow++) {
      vecfunc(x+(tNow* *numberOfEquations),params,shockVec,
      fmats[tNow],fmatsj[tNow],fmatsi[tNow]);
        for (i=0;
        i<=fmatsi[tNow][*numberOfEquations]-fmatsi[tNow][0];i++)
            normSum += SQR(fmats[0][i]);
}
@}


@d apply q terminal constraint
@{

copyMatrix(rowDim,aOne,qMat,qMatj,qMati,aOne,
smats[*pathLength],smatsj[*pathLength],smatsi[*pathLength]);

/*
copmat_(rowDim,qMat,qMatj,qMati,
smats[*pathLength],smatsj[*pathLength],smatsi[*pathLength],
aOne,aOne);
*/

/*xxxxxxxxx add code for deviations using gmat*/
for(i=0;i<*numberOfEquations* (*lags+ *leads);i++){
deviations[i]=x[*numberOfEquations* *pathLength+i]-fixedPoint[i];}
sparseMatTimesVec(rowDim,smats[*pathLength],smatsj[*pathLength],smatsi[*pathLength],deviations,fullfvec);
dnsToCsr(rowDim,aOne,rowDim,fullfvec,
aOne,
fmats[*pathLength],fmatsj[*pathLength],fmatsi[*pathLength],ierr);
for(i=0;i<*rowDim;i++){
  fvec[(*pathLength * *numberOfEquations)+i]=fullfvec[i];
            normSum += SQR(fullfvec[i]);}

        f= 0.5 * normSum;
for(tNow=0;tNow<*pathLength;tNow++) {
        csrToDns(numberOfEquations,
        aOne,fmats[tNow],fmatsj[tNow],fmatsi[tNow],
        fvec+(tNow * *numberOfEquations),numberOfEquations,ierr);
        }
/* modification end */
@}



@d Anderson-Moore algorithm  variable deallocations
@{
free(asymptoticLinearization);
free(cond);
free(epsi);
free(inform);
free(iq);
free(itsbad);
free(nbig);
free(nexa);
free(nnum);
free(nroot);
free(rooti);
free(rootr);
free(uprbnd);
free(qColumns);
free(hColumns);
@}


@d Anderson-Moore algorithm  variable allocations
@{
qColumns=(unsigned int *)calloc(1,sizeof(unsigned int));
*qColumns=(*numberOfEquations)*((*lags)+(*leads));
hColumns=(unsigned int *)calloc(1,sizeof(unsigned int));
*hColumns=(*numberOfEquations)*((*lags)+(*leads)+1);
asymptoticLinearization=(double *)
   calloc( (*numberOfEquations)*((*numberOfEquations)*((*leads)+(*lags)+1)),
   sizeof(double));

cond = (double *)calloc(1,sizeof(double));
epsi = (double *)calloc(1,sizeof(double));
inform = (unsigned int *)calloc(1,sizeof(unsigned int));
iq = (unsigned int *)calloc(1,sizeof(unsigned int));
itsbad = (unsigned int *)calloc(1,sizeof(unsigned int));
nbig = (unsigned int *)calloc(1,sizeof(unsigned int));
nexa = (unsigned int *)calloc(1,sizeof(unsigned int));
nnum = (unsigned int *)calloc(1,sizeof(unsigned int));
nroot = (unsigned int *)calloc(1,sizeof(unsigned int));
rootr=(double *) calloc((*numberOfEquations)*((*lags)+(*leads)),sizeof(double));
rooti=(double *) calloc((*numberOfEquations)*((*lags)+(*leads)),sizeof(double));
uprbnd = (double *)calloc(1,sizeof(double));
@}


@d apply the Anderson-Moore algorithm  to get AMqMatrix
@{
#define NEGLIGIBLEDOUBLE 1.0e-9
#define AIMROOTBOUND 1.0
  *cond = NEGLIGIBLEDOUBLE ;
  *epsi = NEGLIGIBLEDOUBLE ;
  *uprbnd = AIMROOTBOUND + NEGLIGIBLEDOUBLE; 
  *iq = 0;
/*
faim_(asymptoticLinearization,numberOfEquations,lags,leads,epsi, cond, uprbnd,
AMqMatrix,iq,rootr, rooti,nroot,nexa, nnum, nbig, itsbad, inform);
*/
/*
for(i=0;i<*qColumns;i++){printf("roots,(%f,%f)\n",rootr[i],rooti[i]);}
*/
printf("q matrix has %d rows.\n",*iq);
printf("q matrix has %d rows associated with large roots.\n",*nbig);
printf("q matrix has %d rows associated with auxiliary initial conditions.\n",
*nnum+ *nexa);
if(*inform == 0)
{printf("Caller has terminated normally with inform =%d.\n",*inform);}
else {  printf("Caller has terminated with inform =%d.\n",*inform);}
@|NEGLIGIBLEDOUBLE AIMROOTBOUND
@}

@d alt apply the Anderson-Moore algorithm  to get AMqMatrix
@{
auxInit=qRows=0;
maxHElementsForSparseAMA=maxHElements;
/*void * aPointerToVoid;adding since all the sparseAMA.h files have this arg*/
printf("from stackC.w line 2194\n");
cPrintSparse(*numberOfEquations,smats[0],smatsj[0],smatsi[0]);
sparseAMA(&maxHElementsForSparseAMA,DISCRETE_TIME,*numberOfEquations,
*numberOfEquations*(*lags+1+*leads),*leads,
smats[0],smatsj[0],smatsi[0],
newH,newHj,newHi,
&auxInit,&qRows,
qMat,qMatj,qMati,
&essential,rootr,rooti,&returnCode);
@|NEGLIGIBLEDOUBLE AIMROOTBOUND
@}

\subsection{computeAsymptoticQMatrix}
\label{sec:computeAsymptoticQMatrix}


@d create asymptotic linearization for Anderson-Moore algorithm
@{
dfunc(canadaFP,params,
smats[0],smatsj[0],smatsi[0]);






/**hColumns=(*lags+*leads+1)* *numberOfEquations;*/

@}

@o stackC.h -d
@{
void computeAsymptoticQMatrix(@<computeAsymptoticQMatrix argument list@>);
@}


@d computeAsymptoticQMatrix argument list
@{
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
void (* func)(),void (* dfunc)(),double * params,
double canadaFP[],unsigned int * pthLngth,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
double * AMqMatrix,
unsigned int * ierr
@}


@o myNewt.c -d
@{
void computeAsymptoticQMatrix(@<computeAsymptoticQMatrix argument list@>)
{
@<computeAsymptoticQMatrix variable declarations@>
@<computeAsymptoticQMatrix variable allocations@>
  *ierr=0;
@<Anderson-Moore algorithm  variable allocations@>
@<create asymptotic linearization for Anderson-Moore algorithm@>
@<apply the Anderson-Moore algorithm  to get AMqMatrix@>
@<Anderson-Moore algorithm  variable deallocations@>
@<computeAsymptoticQMatrix variable deallocations@>
}
@}
@d computeAsymptoticQMatrix variable deallocations
@{
@}
@d computeAsymptoticQMatrix variable allocations
@{
@}
@d computeAsymptoticQMatrix variable declarations
@{
double * asymptoticLinearization;
double * cond;
double * epsi;
unsigned int * inform;
unsigned int * iq;
unsigned int * itsbad;
unsigned int * nbig;
unsigned int * nexa;
unsigned int * nnum;
unsigned int * nroot;
/*double * zeroVector;*/
unsigned int * qColumns;
double * rootr;
double * rooti;
double * uprbnd;
/*double * wts;*/
/*double * err;*/
unsigned int * hColumns;
/*unsigned int i;*/
@}
@d altComputeAsymptoticQMatrix argument list
@{
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
void (* func)(),void (* dfunc)(),double * params,
double canadaFP[],unsigned int * pthLngth,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
double * qMat,unsigned int * qMatj,unsigned int * qMati,
unsigned int * ierr
@}

@o stackC.h -d
@{
void altComputeAsymptoticQMatrix(@<altComputeAsymptoticQMatrix argument list@>);
@}

@o myNewt.c -d
@{
void altComputeAsymptoticQMatrix(@<altComputeAsymptoticQMatrix argument list@>)
{
@<altComputeAsymptoticQMatrix variable declarations@>
@<altComputeAsymptoticQMatrix variable allocations@>
  *ierr=0;

@<create asymptotic linearization for Anderson-Moore algorithm@>
@<alt apply the Anderson-Moore algorithm  to get AMqMatrix@>

@<altComputeAsymptoticQMatrix variable deallocations@>
}
@}
@d altComputeAsymptoticQMatrix variable deallocations
@{
free(newH);
free(newHj);
free(newHi);
free(rootr);
free(rooti);
@}
@d altComputeAsymptoticQMatrix variable allocations
@{
newH=(double *)calloc(maxHElements,sizeof(double));
newHj=(unsigned int *)calloc(maxHElements,sizeof(unsigned int));
newHi=(unsigned int *)calloc(maxHElements,sizeof(unsigned int));
rootr=(double *) calloc((*numberOfEquations)*((*lags)+(*leads)),
sizeof(double));
rooti=(double *) calloc((*numberOfEquations)*((*lags)+(*leads)),
sizeof(double));

@}
@d altComputeAsymptoticQMatrix variable declarations
@{
unsigned int auxInit;
unsigned int qRows;
unsigned int maxHElements=50000;
unsigned int maxHElementsForSparseAMA;
unsigned int essential;
unsigned int returnCode=0;
/*double * zeroVector;*/
/*int * qColumns;*/
double * newH;
unsigned int * newHj;
unsigned int * newHi;
double * rootr;
double * rooti;
/*double * uprbnd;*/
/*double * wts;*/
/*double * err;*/
/*int * hColumns;*/
/*int i;*/
@}


@d pathNewt undefines
@{

#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef PFREERETURN
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software '>9m_L31.. */

@}
@d pathNewt check for convergence
@{
		test=0.0;
		for (i=0;i<n;i++)
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
		if (test < TOLF) {
			*check=0;
			PFREERETURN
		}
		if (*check) {
			test=0.0;
			den=FMAX(f,0.5*n);
			for (i=0;i<n;i++) {
				temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TOLMIN ? 0 : 1);
			PFREERETURN
		}
		test=0.0;
		for (i=0;i<n;i++) {
			temp=(fabs(x[i]-xold[i]))/(FMIN(fabs(x[i]),fabs(xold[i]))+GAMMA);
			if (temp > test) test=temp;
		}
		if (test < TOLX) PFREERETURN

@}
@d pathNewt check for convergence 
@{
/*
	test=0.0;
	for (i=0;i<*numberOfEquations;i++)
		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	if (test<0.01*TOLF) PFREERETURN
*/
	for (sum=0.0,i=0;i<*numberOfEquations;i++) sum += SQR(x[i]);
	stpmax=(*lags+*leads+1)*STPMX*FMAX(sqrt(sum),(double)n);

@}


@d pathNewt update path
@{
nxtGuess(numberOfEquations,lags,leads,pathLength,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,fixedPoint,x,shockVec,xdel);
		for (i=0;i<n;i++) xoldls[i]=x[i];
		fold=f;


/*printf("here is path check=%d, and f=%f and x before\n",*check,f);
cPrintMatrixNonZero(n,1,x);
*/

if(*pathLength)  {
		lnsrch(n,*numberOfEquations,*pathLength,
        xoldls,&fold,g,xdel,params,shockVec,&f,stpmax,check,vecfunc,x);
        } else {
for(i=*numberOfEquations* *lags;i<n;i++)x[i]=x[i]-xdel[i];
}

/*
printf("here is path check=%d, and f=%f and x after\n",*check,f);
cPrintMatrixNonZero(n,1,x);
*/

/*
for(i=*numberOfEquations* *lags;i<n;i++)x[i]=x[i]-xdel[i];
*/
@}
@d pathNewt initializations
@{
  *check=0;
    n=*numberOfEquations*(*lags+*leads+*pathLength);
	rowDim=(unsigned int *)calloc(1,sizeof(unsigned int));
	*rowDim=*numberOfEquations**leads;
    aOne=(unsigned int *)calloc(1,sizeof(unsigned int));
    *aOne=1;
    qColumns=(unsigned int *)calloc(1,sizeof(unsigned int));
    *qColumns=*numberOfEquations*(*leads+*lags);
    deviations=(double *)calloc(*qColumns,sizeof(double));
    fullfvec=(double *)calloc(*rowDim,sizeof(double));
	ierr=(unsigned int *)calloc(1,sizeof(unsigned int));
	indx=(unsigned int *)calloc(n+1,sizeof(unsigned int));
	g=(double *)calloc(n+1,sizeof(double));
	p=(double *)calloc(n+1,sizeof(double));
	xdel=(double *)calloc(n+1,sizeof(double));
	xold=(double *)calloc(n+1,sizeof(double));
	xoldls=(double *)calloc(n+1,sizeof(double));
	fvec=(double *)calloc(n+1,sizeof(double));


	nn=n;

@}
@d pathNewt declarations
@{    unsigned int  n; 
    unsigned int tNow;unsigned int * rowDim;unsigned int * qColumns;
    double * deviations;double * fullfvec;
/*	double fmin(double x[]);*/
	void lnsrch(unsigned int  n, unsigned int np,unsigned int reps,double xold[], double  *fold, double g[], double p[],double * params,double * shockVec,
		 double * f, double stpmax, unsigned int *check, void (*func)(),double * x);


	unsigned int i,its/*,j*/,*indx,*aOne/*,*ndns*/,*ierr;
	double /*d,*/den,f,fold,stpmax,sum,temp,test,*g,*p,*xold,*xoldls,*xdel,normSum;
@}


\subsection{lnsrch}
\label{sec:lnsrch}

@o myNewt.c -d
@{
#include <math.h>
#define NRANSI
/*#include "./nrutil.h"*/
#define ALF 1.0e-4
#define TOLX 1.0e-10

extern void dgemm_(char * noTransp1,char * noTransp2,
unsigned int * neq,unsigned int * aOne,unsigned int *BMatrixColumns1,
double *aFloatOne,
double *__asymptoticBmatrices,unsigned int *BMatrixRows1,
double *deviationsFromPeriodicPath,unsigned int *BMatrixColumns2,
double *aZero,
double *pathNow,
unsigned int *BMatrixRows2);


extern double  dnrm2_();
#include <string.h>

void lnsrch(unsigned int  n,unsigned int np,unsigned int reps,
double xold[], double * fold, double g[], double p[], 
		 double * params,double * shockVec,double * f,double stpmax, unsigned int *check,
			void (*func)(double*,double*,double*,double*,unsigned int*,unsigned int*),double * x)
{
@< lnsrch declarations@>
@< lnsrch preloop@>
@< lnsrch loop@>
@< lnsrch postloop@>
free(aOne);free(aZero);free(aTwo);free(fvec);free(fvecj);free(fveci);
}
#undef ALF
#undef TOLX
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software '>9m_L31.. */

@}

@d lnsrch declarations
@{
	unsigned int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;
	unsigned int * aOne;unsigned int * aTwo,*ierr,*aZero;
	double * fnow,*fvec/*,*fnorm*/, *xorig,*aDoubleZero;
    unsigned int * fvecj;
    double *dir,*aDoubleOne;
    char * transp, *noTransp;
	unsigned int * fveci;
    transp = (char *)calloc(2,sizeof(char));
    strcpy(transp,"T");
    noTransp = (char *)calloc(2,sizeof(char));
    strcpy(noTransp,"N");
    dir=(double *)calloc(1,sizeof(double));
    aDoubleOne=(double *)calloc(1,sizeof(double));
    *aDoubleOne=1.0;
    fnow=(double *)calloc(n,sizeof(double));
    aDoubleZero=(double *)calloc(1,sizeof(double));
    *aDoubleZero=0.0;
	ierr=(unsigned int *)calloc(1,sizeof(unsigned int));
	fveci=(unsigned int *)calloc(n+1,sizeof(unsigned int));
	aOne=(unsigned int *)calloc(1,sizeof(unsigned int));
	*aOne=1;
	aZero=(unsigned int *)calloc(1,sizeof(unsigned int));
	*aOne=0;
	aTwo=(unsigned int *)calloc(1,sizeof(unsigned int));
	*aTwo=2;
	fvec=(double *)calloc(n,sizeof(double));
	xorig=(double *)calloc(n,sizeof(double));
	fvecj=(unsigned int *)calloc(n,sizeof(double));

@}
@d lnsrch preloop
@{


/*
    fold = dnrm2_(&np,foldp,aOne);
*/

*fold=0;
		for(i=0;i<reps;i++){
((*func)(xold+(i*np),params,shockVec,fvec,fvecj,fveci));
        useCNRMS(&np,aZero,fvec,fvecj,fveci,xorig);
        *fold += xorig[0];
        }
      *fold *= 0.5;

/*
csrToDns(&n,aOne,fvec,fvecj,fveci,xorig,&n,ierr);

    dgemm_(transp,noTransp,
           aOne,aOne,&n,aDoubleOne,
           xorig,&n,
           p,&n,
           aDoubleZero,dir,aOne);
if(*dir<0)
*/

for (i=0;i<n;i++) p[i]= (-p[i]);

	*check=0;
	for (sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=0;i<n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=0;i<n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	if(test!=0){alamin=TOLX/test;} else {alamin=1.0e16;}
	alam=1.0;
@}
@d lnsrch free storage
@{

			free(transp);
            free(ierr);
			free(noTransp);
			free(dir);
			free(aDoubleZero);
			free(aDoubleOne);
			free(fnow);
			free(xorig);
			free(fvec);
			free(fvecj);
			free(fveci);
			free(aOne);
			free(aZero);
			free(aTwo);
@}

@d lnsrch loop
@{
	for (;;) {
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
		for(i=0;i<1;i++)((*func)(x+(i*np),params,shockVec,fvec,fvecj,fveci));
        useCNRMS(&np,aTwo,fvec,fvecj,fveci,xorig);*f = xorig[0] *  0.5;
/*printf("fnormvals %f\n",*f);*/
		if (alam < alamin) {
		  /*			for (i=0;i<n;i++) x[i]=xold[i];*/
			*check=1;
@<lnsrch free storage@>
			return;
		} else if (*f <= *fold+ALF*alam*slope) {
@<lnsrch free storage@>
		  return;}
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-*fold-slope));
			else {
				rhs1 = *f-*fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0) tmplam=0.5*alam;
					else tmplam=(-b+sqrt(disc))/(3.0*a);
			}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=*fold;
		alam=FMAX(tmplam,0.1*alam);
	}

@}
@d lnsrch postloop
@{
@<lnsrch free storage@>
@}



\section{Development Notes}
\label{sec:devnotes}


\subsection{Known Bugs}
\label{sec:knownbugs}


\subsection{Things to Do}
\label{sec:thingstodo}

\begin{itemize}
\item code does not actually allow starting with rows in Q non zero even though it is an argument
\item shocks
\item hands off stochastic sims
\item nl aim really non linear except for stability
\end{itemize}

\appendix
\section{Appendix}
\label{sec:appendix}
\begin{description}
\item[{\bf  FPNewt}] FPNewt(NEQS,NLAGS,NLEADS,
theSparseFunction, theSparseFunctionDerivative, paramVector, 
fp, 
fmats, fmatsj, fmatsi, 
smats, smatsj, smatsi, 
maxNumberElements, 
ierr)
\item[{\bf  computeAsymptoticQMatrix}] computeAsymptoticQMatrix(NEQS, NLAGS, NLEADS, 
theSparseFunction, theSparseFunctionDerivative, paramVector, 
fp, 
fmats, fmatsj, fmatsi, 
smats, smatsj, smatsi, 
maxNumberElements, 
linearConstrainMatrix, 
ierr
);
)
\item[{\bf  pathNewt}] pathNewt(NEQS, NLAGS, NLEADS, pathLength, 
theSparseFunction, theSparseFunctionDerivative, paramVector, 
fmats, fmatsj, fmatsi, 
smats, smatsj, smatsi, 
maxNumberElements, 
linearConstrainMatrix, fp
path, 
ierr)
\end{description}



\subsection{Purify Output}
\label{sec:purify}

\begin{verbatim}
 Copyright (C) 1992-1996 Pure Software Inc. All rights reserved. 
  * For contact information type: "purify -help"
  * Command-line: pureStack 
  * Options settings: -purify -windows=no -log-file=stackViewPurify \
    -purify-home=/opt/pure/bin//../purify 
PureLA: 3 simple licenses, 4 users.  Please remedy.
  * Purify licensed to Federal Reserve Board
  * Purify checking enabled.

****  Purify instrumented pureStack (pid 23033)  ****
Current file descriptors in use: 5
FIU: file descriptor 0: <stdin>
FIU: file descriptor 1: <stdout>
FIU: file descriptor 2: <stderr>
FIU: file descriptor 26: <reserved for Purify internal use>
FIU: file descriptor 27: <reserved for Purify internal use>

****  Purify instrumented pureStack (pid 23033)  ****
Purify: Searching for all memory leaks...

Memory leaked: 0 bytes (0%); potentially leaked: 0 bytes (0%)

Purify Heap Analysis (combining suppressed and unsuppressed chunks)
                   Chunks      Bytes
              Leaked          0          0
  Potentially Leaked          0          0
              In-Use          2     100144
  ----------------------------------------
     Total Allocated          2     100144

****  Purify instrumented pureStack (pid 23033)  ****
  * Program exited with status code 1.
  * 0 access errors, 0 total occurrences.
  * 0 bytes leaked.
  * 0 bytes potentially leaked.
  * Basic memory usage (including Purify overhead):
     1135703 code
      111809 data/bss
     5644295 heap (peak use)
        3336 stack
  * Shared library memory usage (including Purify overhead):
      805145 libc.so.1_pure_p1_c0_032_551.so.1 (shared code)
       34908 libc.so.1_pure_p1_c0_032_551.so.1 (private data)
        1204 libdl.so.1_pure_p1_c0_032_551.so.1 (shared code)
           0 libdl.so.1_pure_p1_c0_032_551.so.1 (private data)
      364984 libF77.so.3_pure_p1_c0_032_551.so.3 (shared code)
       69574 libF77.so.3_pure_p1_c0_032_551.so.3 (private data)
        1889 libinternal_stubs.so.1 (shared code)
         208 libinternal_stubs.so.1 (private data)


\end{verbatim}

\subsection{Diagnostic Printing}
\label{sec:diag}

Use emacs to enter:
\begin{verbatim}
gdb testStack -x gdbCommands
\end{verbatim}

@o gdbCommands
@{
break nxtCDmats
commands
p "cmats[0] has a array"
p *(cmats[0])@@12
p *(cmatsj[0])@@12
p *(cmatsi[0])@@12
p "dmats[0] has a array"
p *(dmats[0])@@12
p *(dmatsj[0])@@12
p *(dmatsi[0])@@12
p "cmats[1] has a array"
p *(cmats[1])@@12
p *(cmatsj[1])@@12
p *(cmatsi[1])@@12
p "dmats[1] has a array"
p *(dmats[1])@@12
p *(dmatsj[1])@@12
p *(dmatsi[1])@@12
p "cmats[2] has a array"
p *(cmats[2])@@12
p *(cmatsj[2])@@12
p *(cmatsi[2])@@12
p "dmats[2] has a array"
p *(dmats[2])@@12
p *(dmatsj[2])@@12
p *(dmatsi[2])@@12
p "cmats[3] has a array"
p *(cmats[3])@@12
p *(cmatsj[3])@@12
p *(cmatsi[3])@@12
p "dmats[3] has a array"
p *(dmats[3])@@12
p *(dmatsj[3])@@12
p *(dmatsi[3])@@12
end
run
@}




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




%\bibliographystyle{authordate2}
\bibliographystyle{plainnat}
\bibliography{files,anderson}
\end{document}

