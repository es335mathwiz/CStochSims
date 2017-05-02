%$Id: stackC.w,v 1.44 2002/06/10 19:21:25 m1gsa00 Exp m1gsa00 $
\documentclass[html]{article} \include{miscLatexPkg}
\newcommand{\Treebox}[1]{ \Tr{\psframebox{#1}}} \begin{document}

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
@<define constants and specify include files@>
@<define assert bump@>
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
#include "stochSims.h"
//#include "stochProto.h"
#include "useSparseAMA.h"

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
#include <stdio.h>
#include <stdlib.h>
//void free(void * ptr);
//void * calloc(size_t amt,size_t size);
void nxtCDmats(@<nxtCDmats argument list@>);
@}
Function uses SPARSEKIT's CSR format.
Since HARWELL functions expects CSC,
MA50CD operates on the transpose.

@d nxtCDmats definition
@{
#include <unistd.h>
//pid_t getpid(void);

void nxtCDmats(@<nxtCDmats argument list@>){
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
double **dmatsA,unsigned int **dmatsJA,unsigned int **dmatsIA,
unsigned int * ma50bdJob,
unsigned int * ma50bdIq,
double * ma50bdFact,
unsigned int * ma50bdIrnf,
unsigned int * ma50bdIptrl,
unsigned int * ma50bdIptru
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


% \psset{arrows=<-}
% \pstree[treemode=U]{\Tcircle{evensumC}}{\pstree{\Treebox{APLB\_}}{\pstree[treemode=U]{\Tcircle{b}}{\pstree{\Treebox{AMUB\_}}{\Tcircle{ao}\Tcircle{cmats}}}\Tcircle{oddsumC}}}

\vspace{1.0cm}



\begin{programbox}
  |(b,jb,ib)| =  H_{-\tau+i} C_{-\tau+i}
  |(evenSumCA,evenSumCJA,evenSumCIA)| =  |(oddSumCA,oddSumCJA,oddSumCIA)| - |(b,jb,ib)|
\end{programbox}



@d multiply c matrices by appropriate s matrix and subtract
@{
sparseMult(rowDim,cColumns,aOne,ao,jao,iao,
(cmatsA[timeOffset]),(cmatsJA[timeOffset]),(cmatsIA[timeOffset]),
b,jb,ib,nzmax,iw,ierr);
pathNewtAssert(*ierr == 0);
bump((cmatsIA[timeOffset])[*rowDim]-(cmatsIA[timeOffset])[0]);
aSmallDouble=DBL_EPSILON;
dropSmallElements(rowDim,aOne,&aSmallDouble,nzmax,b,jb,ib,b,jb,ib,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);
/*actually want to subtract so mult elements by -1 also need to shift right*/
for(j=0;j<ib[*rowDim]-1;j++)
{b[j]=(-1)*b[j];jb[j]=jb[j]+(*numberOfEquations*(timeOffset+*lagss+1));};

sparseAdd(rowDim,cColumns,
nzmax,iw,aOne,oddSumCA,oddSumCJA,oddSumCIA,
b,jb,ib,evenSumCA,evenSumCJA,evenSumCIA,ierr);
pathNewtAssert(*ierr == 0);
bump(evenSumCIA[*rowDim]-evenSumCIA[0]);
@}


% \psset{arrows=<-}
% \pstree[treemode=U]{\Tcircle{evensumD}}{\pstree{\Treebox{APLB\_}}{\pstree[treemode=U]{\Tcircle{b}}{\pstree{\Treebox{AMUB\_}}{\Tcircle{ao}\Tcircle{dmats}}}\Tcircle{oddsumD}}}

\vspace{1.0cm}


\begin{programbox}
  |(b,jb,ib)| =  H_{-\tau+i} C_{-\tau+i}
  |(evenSumDA,evenSumDJA,evenSumDIA)| =  |(oddSumDA,oddSumDJA,oddSumDIA)| - |(b,jb,ib)|
\end{programbox}



@d multiply d matrices by appropriate s matrix and subtract
@{

sparseMult(rowDim,aOne,aOne,ao,jao,iao,
(dmatsA[timeOffset]),(dmatsJA[timeOffset]),(dmatsIA[timeOffset]),
b,jb,ib,nzmax,iw,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);
aSmallDouble=DBL_EPSILON;
dropSmallElements(rowDim,aOne,&aSmallDouble,nzmax,b,jb,ib,b,jb,ib,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);

/*actually want to subtract so mult elements by -1*/
for(j=0;j<ib[*rowDim]-1;j++)b[j]=(-1)*b[j];

sparseAdd(rowDim,aOne,
nzmax,iw,aOne,oddSumDA,oddSumDJA,oddSumDIA,
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
% \pstree[treemode=U]{\Tcircle{oddsumC}}{\pstree{\Treebox{APLB\_}}{\Tcircle{evensumC}\pstree{\Tcircle{b}}{\pstree{\Treebox{SUBMAT\_}}{\Tcircle{smats}}}}}
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

*firstColumn=(*numberOfEquations* *lagss)+1;
*lastColumn=*firstColumn + *rowDim-1;
extractSubmatrix(rowDim,aOne,aOne,rowDim,
firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,nr,nc,
oddSumCA,oddSumCJA,oddSumCIA);
*nonZeroNow=oddSumCIA[*rowDim]-oddSumCIA[0];

useMA50ID(cntl,icntl);
useMA50AD(rowDim,rowDim,nonZeroNow,nzmax,oddSumCA,oddSumCJA,jcn,oddSumCIA,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
/*wordybump(info[3]);*/
//pathNewtAssert(info[0]>=0);

@}


@d factorize matrix
@{
/* restore odd since ad is destructive*/
extractSubmatrix(rowDim,aOne,aOne,rowDim,
firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,nr,nc,
oddSumCA,oddSumCJA,jcn);
if(*ma50bdJob==1){
useMA50BD(rowDim,rowDim,nonZeroNow,ma50bdJob,
oddSumCA,oddSumCJA,jcn,
cntl,icntl,ip,oddSumCIA,np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump(info[3]);
//pathNewtAssert(info[0]>=0);
for(i=0;i<*rowDim+1;i++){ma50bdIq[i]=oddSumCIA[i];}
} else {
useMA50BD(rowDim,rowDim,nonZeroNow,ma50bdJob,
oddSumCA,oddSumCJA,jcn,
cntl,icntl,ip,ma50bdIq,np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
/*wordybump(info[3]);*/
//pathNewtAssert(info[0]>=0);
}
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
for(i=0;i<*balColumns;i++){


*lastColumn = *firstColumn=(1+i+ *rowDim +*numberOfEquations *(*lagss));

extractSubmatrix(rowDim,aOne,
aOne,rowDim,firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,
nr,nc,
b,jb,ib);

csrToDns(rowDim,aOne,b,jb,ib,
nsSumC,rowDim,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);
useMA50CD(rowDim,rowDim,icntl,ma50bdIq,np,trans,
lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
nsSumC,x,
w,info);
bump(info[3]);
//pathNewtAssert(info[0]>=0);

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
dropSmallElements(balColumns,aOne,&aSmallDouble,nzmax,tb,jtb,itb,tb,jtb,itb,ierr);
pathNewtAssert(*ierr == 0);
csrToCsc(balColumns,aOne,aOne,tb,jtb,itb,cmatsA[0],cmatsJA[0],cmatsIA[0]);

/*expand sum of d's*/
csrToDns(rowDim,aOne,oddSumDA,oddSumDJA,oddSumDIA,nsSumD,rowDim,ierr);
pathNewtAssert(*ierr == 0);
bump(*rowDim);
/*code should use info from previous call to set lfact
also can avoid calls to ma50ad once pattern settles down*/

useMA50CD(rowDim,rowDim,icntl,ma50bdIq,np,
trans,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,nsSumD,x,w,info);
//pathNewtAssert(info[0]>=0);
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
unsigned int timeOffset;
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
cPrunsigned intMatrixNonZero
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
*cColumns = *numberOfEquations * (1+(*leadss?*leadss:1));
*balColumns = *cColumns-*rowDim;
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
x = (double *)calloc(*rowDim * *numberOfEquations * (*leadss+1),sizeof(double));
trans = (unsigned int *) calloc(1,sizeof(unsigned int));
nsSumD = (double *)calloc(*rowDim,sizeof(double));
nsSumC = (double *)calloc(*rowDim /** *numberOfEquations * (*leadss+1+ *lagss)*/,sizeof(double));
aOne = (unsigned int *)calloc(1,sizeof(unsigned int));
*aOne = 1;
@}


@d nxtCDmats array variable allocations
@{

ib = (unsigned int *)calloc(*rowDim *(1+ *leadss) + 1,sizeof(unsigned int));
jb = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
b = (double *)calloc(*maxNumberHElements,sizeof(double));
itb = (unsigned int *)calloc(*cColumns + 1,sizeof(unsigned int));
jtb = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
tb = (double *)calloc(*maxNumberHElements,sizeof(double));


evenSumCIA = (unsigned int *)calloc(( *rowDim * (1+*leadss))+1,sizeof(unsigned int));
evenSumCJA = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
evenSumCA = (double *)calloc(*maxNumberHElements,sizeof(double));


/*larger than necessary now so that can use for transpose in csrcsc */
oddSumCIA = (unsigned int *)calloc(( *rowDim * (1+*leadss))+ 1,sizeof(unsigned int));
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
iw = (unsigned int *)calloc(*rowDim * ((1+*leadss) + *lagss + 3),sizeof(unsigned int));
w = (double *)calloc(*rowDim * ((1+*leadss) + *lagss + 3),sizeof(double));/*MLK*/@|
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
void oneStepBack(@<oneStepBack argument list@>)
{
@<oneStepBack variable definitions@>
@<oneStepBack variable allocations@>
/*cmat non zero then multiply else product is zero*/
if(cmatsIA[0][*rowDim]-cmatsIA[0][0]) {
  sparseMult(rowDim,aOne,aOne,cmatsA[0],cmatsJA[0],cmatsIA[0],
  yvecA[0+1 ],yvecJA[0+1 ],yvecIA[0+1 ],
  rcy,rcyj,rcyi,rowDim,iw,ierr);
pathNewtAssert(*ierr == 0);

aSmallDouble=DBL_EPSILON;

dropSmallElements(rowDim,aOne,&aSmallDouble,nzmax,rcy,rcyj,rcyi,rcy,rcyj,rcyi,ierr);
pathNewtAssert(*ierr == 0);

  for(i=0;i<rcyi[*rowDim]-rcyi[0];i++)rcy[i]=(-1)*rcy[i];
  sparseAdd(rowDim,aOne,
  nzmax,iw,aOne,
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
void (*func)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
void (*dfunc)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * params,
double * expansionPounsigned int,
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
void newNxtGuess(@<newNxtGuess argument list@>);
@}


@d nxtGuess argument list
@{
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads, unsigned int * capT,
double **fmats,unsigned int  **fmatsj,unsigned int  **fmatsi,
double **smats,unsigned int  **smatsj,unsigned int  **smatsi,
unsigned int *maxNumberHElements,
double * termConstr,unsigned int * termConstrj,unsigned int * termConstri,double * fp,double * intercept,
double * initialX,
/*double * shockVec,*/
double * updateDirection/*,
unsigned int * intControlParameters,double * doubleControlParameters*/
 @| termConstr fp initialX shockVec
theFunc theDrvFunc capT 
@}

@o stackC.h -d
@{
void constructFdrv(@<constructFdrv argument list@>);
@}

@d constructFdrv argument list
@{
unsigned int numberOfEquations,unsigned int lags, unsigned int leads,unsigned int pathLength,
double * xvec,double * params,
void (*vFunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
void (*vFuncDrv)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * termConstr,unsigned int * termConstrj,unsigned int * termConstri,
double * fixedPoint,double * intercept,double * linearizationPoint,
/*unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,*/
double * shockVec,
double * fvec,
double * fdrv,unsigned int * fdrvj,unsigned int * fdrvi,unsigned int ihomotopy,
unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo/*, double * doubleOutputInfo*/@}

@d chkDrv definition
@{
#include <math.h>
#define NEGLIGIBLEDOUBLE 1.0e-9
void chkDrv(unsigned int n, double * fdrv,unsigned int * fdrvj,unsigned int * fdrvi,
double * fvec,double * delxvec)
{
unsigned int i;
//unsigned int aOne=1;
double * fvals;
fvals = (double * ) calloc(n,sizeof(double));
#ifdef DEBUG 
printf("chkDrv:beginning\n");
#endif


sparseMatTimesVec(&n,delxvec,fvals,fdrv,fdrvj,fdrvi);
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

void constructFdrv(@<constructFdrv argument list@>)
{
double * deviations;
unsigned int * ignore;//double  dignore[1]={1.0};
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
fvec+numberOfEquations*lags,fvecj,fveci,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
addOneToFDrvEvals;
vFuncDrv(xvec,params,shockVec,
fdrv+soFar,fdrvj+soFar,fdrvi+numberOfEquations*lags,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
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
fvec+numberOfEquations*lags+i*numberOfEquations,fvecj,fveci,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
addOneToFDrvEvals;
vFuncDrv(xvec+i*numberOfEquations,params,zeroShockVec,
fdrv+soFar,fdrvj+soFar,fdrvi+numberOfEquations*lags+(i*numberOfEquations),homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
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
sparseMatTimesVec(&rowDim,termConstrj,termConstri,fvec+(numberOfEquations*(lags+pathLength)),deviations,
termConstr);
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
/*unsigned int * ma50bdIq,*/
double * ma50bdFact,
unsigned int * ma50bdIrnf,
unsigned int * ma50bdIptrl,
unsigned int * ma50bdIptru,
unsigned int * intControlParameters,double * doubleControlParameters
@}

@d nxtGuess definition
@{
#define sysDimSwitchLevel 30000
#define ma50DropThreshold (1e-8)

void newNxtGuess(@<newNxtGuess argument list@>)
{


//unsigned int i;
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
int * info;
double * rinfo;
unsigned int *lfact;
//double * fact;
//unsigned int *irnf;
//unsigned int * iptrl;
//unsigned int * iptru;
unsigned int *iw;
double * w;
unsigned int nzmax;
unsigned int  * aOne;
unsigned int nonZeroNow;
unsigned int  trans;
@}
@d nxtGuess definition
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
info = (int *)calloc(7,sizeof( int));
rinfo = (double *)calloc(1,sizeof(double));
lfact = (unsigned int *)calloc(1,sizeof(unsigned int));
*lfact = ( *maxNumberHElements);/*pessimistic setting for filling*/
aOne = (unsigned int *)calloc(1,sizeof(unsigned int));
*aOne=1;
@}
@d nxtGuess definition
@{
copyMatrix(sysDim,aOne,fdrv,fdrvj,fdrvi,aOne,
copychkfdrv,copychkfdrvj,copychkfdrvi);

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
wordybump((unsigned int)info[3]);
//pathNewtAssert(info[0]>=0);

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
wordybump((unsigned int)info[3]);
#ifdef DEBUG 
printf("\n ma50bd info\n");
for(i=0;i<7;i++)printf(" %d ",info[i]);
printf("\n ma50bd info\n");
#endif
//pathNewtAssert(info[0]>=0);
if(*ma50bdJob==1)*ma50bdJob=1;
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
wordybump((unsigned int)info[3]);
//pathNewtAssert(info[0]>=0);
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
wordybump((unsigned int)info[3]);
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
@d nxtGuess definition
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
*maxNumberHElements=maxElementsEncountered;
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
int iSigned;
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

cmats =(double **)calloc(*capT+(*lags+*leads)+1,sizeof(double *));
cmatsj =(unsigned int **)calloc(*capT+(*lags+*leads)+1,sizeof(unsigned int *));
cmatsi =(unsigned int **)calloc(*capT+(*lags+*leads)+1,sizeof(unsigned int *));
dmats =(double **)calloc(*capT+(*lags+*leads)+1,sizeof(double *));
dmatsj =(unsigned int **)calloc(*capT+(*lags+*leads)+1,sizeof(unsigned int *));
dmatsi =(unsigned int **)calloc(*capT+(*lags+*leads)+1,sizeof(unsigned int *));
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
@}

@d nxtGuess obtain sparse representation and compute next C and d
@{
nxtCDmats(numberOfEquations,lags,leads,
numberOfEquations,maxNumberHElements,
   smats[tNow],smatsj[tNow],smatsi[tNow],
   fmats[tNow],fmatsj[tNow],fmatsi[tNow],
   cmats+*lags+tNow,cmatsj+*lags+tNow,cmatsi+*lags+tNow,
   dmats+*lags+tNow,dmatsj+*lags+tNow,dmatsi+*lags+tNow,
ma50bdJob,
ma50bdIq,
ma50bdFact,
ma50bdIrnf,
ma50bdIptrl,
ma50bdIptru
);
bump(*maxNumberHElements);
@}
@d nxtGuess use terminal constraint
@{




copyMatrix(rowDim,aOne,termConstr,termConstrj,termConstri,aOne,
smats[*capT],smatsj[*capT],smatsi[*capT]);

/*xxxxxxxxx add code for deviations using gmat*/
for(i=0;i<*numberOfEquations* (*lags+ *leads);i++){
deviations[i]=initialX[*numberOfEquations* *capT+i]-fp[i+*numberOfEquations];}
sparseMatTimesVec(rowDim,smats[*capT],smatsj[*capT],smatsi[*capT],deviations,fullfvec);
for(i=0;i<*numberOfEquations* *leads;i++){fullfvec[i]=
fullfvec[i]-intercept[i];}
dnsToCsr(rowDim,aOne,rowDim,fullfvec,
aOne,
fmats[*capT],fmatsj[*capT],fmatsi[*capT],ierr);
pathNewtAssert(*ierr == 0);

#ifdef DEBUG 
/*
printf("nxtGuess:deviations computation\n");
printf("fullfvec={\n");
for(i=0;i<*rowDim-1;i++){printf("%f,",fullfvec[i]);}
printf("%f}\n",fullfvec[*rowDim]);
printf("deviations={\n");
for(i=0;i<*numberOfEquations* (*lags+ *leads)-1;i++){printf("%f,",deviations[i]);}
printf("%f}\n",deviations[*numberOfEquations* (*lags+ *leads)]);
*/
#endif


*ma50bdJob=1;
nxtCDmats(numberOfEquations,lags,leads,
rowDim,maxNumberHElements,
termConstr,termConstrj,termConstri,
fmats[*capT],fmatsj[*capT],fmatsi[*capT],
   cmats+*lags+*capT,cmatsj+*lags+*capT,cmatsi+*lags+*capT,
   dmats+*lags+*capT,dmatsj+*lags+*capT,dmatsi+*lags+*capT,
ma50bdJob,
ma50bdIq,
ma50bdFact,
ma50bdIrnf,
ma50bdIptrl,
ma50bdIptru
);
bump(*maxNumberHElements);
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
for(iSigned=*capT-1;iSigned>-1;iSigned--){
oneStepBack(numberOfEquations,
   ymats+*lags+iSigned,ymatsj+*lags+iSigned,ymatsi+*lags+iSigned,
   cmats+*lags+iSigned,cmatsj+*lags+iSigned,cmatsi+*lags+iSigned,
   dmats+*lags+iSigned,dmatsj+*lags+iSigned,dmatsi+*lags+iSigned);
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
/*double * initialX,*/double * updateDirection @| termConstr fp initialX 
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
//double *fullXvec;
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

cmats =(double **)calloc((*lags+*leads)+2,sizeof(double *));
cmatsj =(unsigned int **)calloc((*lags+*leads)+2,sizeof(unsigned int *));
cmatsi =(unsigned int **)calloc((*lags+*leads)+2,sizeof(unsigned int *));
dmats =(double **)calloc((*lags+*leads)+2,sizeof(double *));
dmatsj =(unsigned int **)calloc((*lags+*leads)+2,sizeof(unsigned int *));
dmatsi =(unsigned int **)calloc((*lags+*leads)+2,sizeof(unsigned int *));
gmats =(double *)calloc(*numberOfEquations,sizeof(double));
gmatsj =(unsigned int *)calloc(*numberOfEquations,sizeof(unsigned int));
gmatsi =(unsigned int *)calloc(*numberOfEquations,sizeof(unsigned int));
ymats =(double **)calloc((*lags+*leads)+2,sizeof(double *));
ymatsj =(unsigned int **)calloc((*lags+*leads)+2,sizeof(unsigned int *));
ymatsi =(unsigned int **)calloc((*lags+*leads)+2,sizeof(unsigned int *));
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
nxtCDmats(numberOfEquations,lags,leads,
numberOfEquations,maxNumberHElements,
   smats,smatsj,smatsi,
   fmats,fmatsj,fmatsi,
   cmats+*lags,cmatsj+*lags,cmatsi+*lags,
   dmats+*lags,dmatsj+*lags,dmatsi+*lags,
ma50bdJob,
ma50bdIq,
ma50bdFact,
ma50bdIrnf,
ma50bdIptrl,
ma50bdIptru
);
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
   dmats+*lags+1,dmatsj+*lags+1,dmatsi+*lags+1,
ma50bdJob,
ma50bdIq,
ma50bdFact,
ma50bdIrnf,
ma50bdIptrl,
ma50bdIptru
);

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
int iSigned;
for(iSigned=*lags;iSigned>-1;iSigned--){
oneStepBack(numberOfEquations,
   ymats+iSigned,ymatsj+iSigned,ymatsi+iSigned,
   cmats+iSigned,cmatsj+iSigned,cmatsi+iSigned,
   dmats+iSigned,dmatsj+iSigned,dmatsi+iSigned);
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


@d compPathError definition
@{
void compPathError(unsigned int * numberOfEquations,unsigned int * lags,unsigned int * leads,
void (* theFunction)(double*, double*, double*,double*,unsigned int*,unsigned int*/*,double *,double**/),
double * termConstr,unsigned int * termConstrj,unsigned int * termConstri,double * fp,double * intercept,
double * initialX,
double * shockVec,
unsigned int * capT,
/*unsigned int * maxNumberHElements,*/
double ** fmats,unsigned int ** fmatsj,unsigned int **fmatsi,
double ** smats,unsigned int ** smatsj,unsigned int **smatsi)
{
unsigned int maxElementsEncountered=0;
unsigned int * rowDim;
unsigned int * qColumns;
unsigned int tNow;
unsigned int * aOne;
unsigned int * aZero;
unsigned int * ierr;unsigned int i;
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

/*xxxxxxxxx add code for deviations using gmat*/
for(i=0;i<*numberOfEquations* (*lags+ *leads);i++){
deviations[i]=initialX[*numberOfEquations* *capT+i]-fp[i+*numberOfEquations];}
sparseMatTimesVec(rowDim,smats[*capT],smatsj[*capT],smatsi[*capT],deviations,fullfvec);
for(i=0;i<*numberOfEquations* *leads;i++){fullfvec[i]=
fullfvec[i]-intercept[i];}
dnsToCsr(rowDim,aOne,rowDim,fullfvec,
aOne,
fmats[*capT],fmatsj[*capT],fmatsi[*capT],ierr);
pathNewtAssert(*ierr == 0);
bump(fmatsi[*capT][*rowDim]-fmatsi[*capT][0]);
free(zeroShockVec);
}


@}


\section{Numerical Recipes Modifications of FPnewt}


\subsection{array allocation program}
\label{sec:allocarray}


@o stackC.c -d
@{
//#include "stochSims.h"
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

@<allocFPNewt signature@>
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
@<freeFPNewt signature@>
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
@<allocAltComputeAsymptoticQ signature@>
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


@<freeAltComputeAsymptoticQ signature@>
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


@<allocPathNewt signature@>
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
@}
@d freeFPNewt signature
@{
void freeFPNewt(unsigned int lags, unsigned int pathLength,
double ** genericFP,
double ** genericIntercept,
double***fmats,unsigned int***fmatsj,unsigned int***fmatsi,
double***smats,unsigned int***smatsj,unsigned int***smatsi)
@}
@o stackC.h -d
@{
@<freeFPNewt signature@>;
@<freeAltComputeAsymptoticQ signature@>;
@}

@d freeAltComputeAsymptoticQ signature
@{
void freeAltComputeAsymptoticQ(
double**AMqMatrix,unsigned int**AMqMatrixj,unsigned int**AMqMatrixi,
double**rootr,double**rooti)
@}
@d allocPathNewt signature
@{
void allocPathNewt(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int pathLength,unsigned int replications,unsigned int stochasticPathLength,
double**genericPath,
double**genericZeroPath,
double**genericEasyPath,
double**genericTargetPath
)@}

@d allocAltComputeAsymptoticQ signature
@{
void allocAltComputeAsymptoticQ(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int maxElements,double**AMqMatrix,unsigned int**AMqMatrixj,unsigned int**AMqMatrixi,
double** rootr,double**rooti)
@}

@d freePathNewt signature
@{
void freePathNewt(double ** genericPath)
@}

@o stackC.h -d
@{
@<allocAltComputeAsymptoticQ signature@>;
@<freePathNewt signature@>;
@<allocPathNewt signature@>;

@}

@d altComputeAsymptoticIMatrix signature
@{
void altComputeAsymptoticIMatrix(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
double * qMat,unsigned int * qMatj,unsigned int * qMati/*,
unsigned int * ierr*/
)
@}

@o stackC.h -d
@{
@<altComputeAsymptoticIMatrix signature@>;
@}

@o stackC.c -d
@{
@<freePathNewt signature@>
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

@d no longer needed
@{

@}

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

@o myNewt.c -d 
@{
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include "stochProto.h"
#include "useSparseAMA.h"
#include "stackC.h"
#include "stochSims.h"

/*taken from numerical recipes nrutil.h*/
@<define assert bump@>
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

@<FPnewt defines@>
@<FPnewt signature@>
{
@<FPnewt declarations@>
*check=0;done=0;
for(ihomotopy=0;(ihomotopy<numberAlphas)&&(!(*check));ihomotopy++){
#ifdef DEBUG 
printf("%d-th homotopy -- %f\n",ihomotopy,*(homotopyAlpha+ihomotopy));
#endif
@<evaluate func and find max element@>
for (its=1;its<=((maxitsInput)&&(!done));its++) {
@<get newton update@>
@<evaluate func update norm@>
@<check for convergence@>
	}
if(!done){
*check=1;
addOneToFailedQ;
 printf("MAXITS=%d exceeded in FPNewt\n",(maxitsInput));}
done=0;
}
if(ignoreFailQ){
if(*check==1){subOneFromFailedQ;}
*check=0;
printf("IGNORING FAILED CONVERGENCE!!!!!!!!\n");
}
FREERETURN
}
@}
@d FPnewt defines
@{
#include <math.h>
#define NRANSI


#define TOLMIN -1.0e-6

#define STPMX 100.0
#define GAMMA 1.0

unsigned int nn;
double *fvec;
#define FREERETURN {*homotopyAlpha=realAlpha; \
free(fvec);free(xold);free(shockVec);\
	free(p);free(g);free(aOne);free(ierr);\
	free(indx);return;}
#define PFREERETURN { assignRealizedTolf = testf;assignRealizedTolx=testx; \
free(fvec);free(xold);free(xoldls);free(xdel);\
free(deviations);free(fullfvec);\
	free(p);free(g);free(aOne);free(aZero);free(ierr);\
	free(indx);free(rowDim);free(qColumns);\
free(compfvec);free(compfvecj);free(compfveci);free(chkfdrv);free(chkfdrvj);free(chkfdrvi);free(zeroShockVec);\
*maxNumberElements = maxElementsEncountered;\
return;}
@}

@d FPnewt signature
@{
void FPnewt(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
void (*func)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
void (*dfunc)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * params,
double x[],double * linearizationPoint,/*unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,*/
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
unsigned int *check,
unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo)@}

@o stackC.h -d
@{
@<FPnewt signature@>;/*first one*/
@}

@d delete FPnewt declarations
@{
    unsigned int n/*;double * xdel;double  dignore[1]={1.0}*/;unsigned int done;
/*	double fmin(double x[]);*/unsigned int ihomotopy=1;double realAlpha;
/*	void lnsrch(unsigned int n,unsigned int np,unsigned int reps,
    double xold[], double   * fold, double g[], double p[],double * params,
		 double * shockVec,double * f, double stpmax, unsigned int *check, 
         void (*func)(double*,double*,double*,double*,unsigned int*,unsigned int*,double*,double*),double *x,
            unsigned int ihomotopy,double * linPt,unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo);*/
	void lubksb(double **a, unsigned int n, unsigned int *indx, double b[]);
	void ludcmp(double **a, unsigned int n, unsigned int *indx, double *d);
	unsigned int i,its/*,j*/,*indx,*aOne/*,*ndns*/,*ierr;
	double /*d,den,*/f,fold,stpmax,sum,temp,test,*g,*p,*xold,normSum;
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
	shockVec=(double *)calloc(*numberOfEquations,sizeof(double));
    for(i=0;i<*numberOfEquations;i++)shockVec[i]=0.0;
	nn=n;
/*save alphas but force fp to use orginal function*/
realAlpha=*homotopyAlpha;
/**homotopyAlpha=1;;*/
resetFailedQ;
resetNewtonSteps;
resetRealizedTolf;
resetRealizedTolx;
resetFEvals;
resetFDrvEvals;
resetShrinkSteps;
resetExpandSteps;
resetLnsrchSteps;

@}


@d evaluate func and find max element
@{


/* modification begin */
      func(x,params,shockVec,fmats[0],fmatsj[0],fmatsi[0],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
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
/*	if (test<0.01*(tolfInput)) FREERETURN*/
	for (sum=0.0,i=0;i<*numberOfEquations;i++) sum += SQR(x[i]);
	stpmax=(*lags+*leads+1)*STPMX*FMAX(sqrt(sum),(double)n);


        sparseMatTimesVec(numberOfEquations,
        smats[0],smatsj[0],smatsi[0],fvec,
        g+(*numberOfEquations * *lags));
        for(i=0;i>*numberOfEquations;i++){
          g[(*numberOfEquations* (*lags-1))+i]=fvec[(*numberOfEquations* *lags)+i];
          g[(*numberOfEquations* (*lags+1))+i]= -fvec[(*numberOfEquations* *lags)+i];
          }

@}

@d get newton update
@{
dfunc(x,params,shockVec,smats[0],smatsj[0],smatsi[0],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
		for (i=0;i<n;i++) xold[i]=x[i];
		/*modification begin*/
nxtFPGuess(numberOfEquations,lags,leads,
fmats[0],fmatsj[0],fmatsi[0],
smats[0],smatsj[0],smatsi[0],
maxNumberElements,/*x,*/p);




		lnsrch(n,*numberOfEquations,1,
        xold,&fold,g,p,params,shockVec,&f,stpmax,check,func,x,ihomotopy,
linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
/*intControlParameters,*/doubleControlParameters/*,
intOutputInfo, doubleOutputInfo*/);


		fold=f;


@}

@d evaluate func update norm
@{

      func(x,params,shockVec,fmats[0],fmatsj[0],fmatsi[0],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
        for (normSum=0.0,i=0;i<=fmatsi[0][*numberOfEquations]+fmatsi[0][0];i++) 
            normSum += SQR(fmats[0][i]);
        f= 0.5 * normSum * (*lags+*leads+1);
        csrToDns(numberOfEquations,
        aOne,fmats[0],fmatsj[0],fmatsi[0],fvec+(*numberOfEquations**lags),
        numberOfEquations,ierr);
pathNewtAssert(*ierr == 0);
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
		if (test < (tolfInput)) {
resetFailedQ;
			*check=0;
done=1;
		} else {
		test=0.0;
		for (i=0;i<n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < tolxInput) {
*check=0;done=1;}
}
@}


\subsection{pathNewt}
\label{sec:pathNewt} 




@d pathNewt signature
@{
void pathNewt(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
void (* fdjac)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * params,double * shockVec,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,double * linearizationPoint,
/*unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,*/
double x[],
unsigned int *check, double * lastDel,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * ma50bdJob,
/*unsigned int * ma50bdIq,*/
double * ma50bdFact,
unsigned int * ma50bdIrnf,
unsigned int * ma50bdIptrl,
unsigned int * ma50bdIptru
)
@}

@o stackC.h -d 
@{
@<pathNewt signature@>;
@}

@o myNewt.c -d 
@{

void exogenizeRow(unsigned int row, double * a,/* unsigned int * ja,*/ unsigned int * ia)
{
//unsigned int i; unsigned int strt;
a[ia[row-1]-1]=0.0;
}
void exogenizeDRow(unsigned int row, unsigned int col, double * a, unsigned int * ja, unsigned int * ia)
{
unsigned int i; unsigned int elems; unsigned int strt; unsigned int fin;
fin=ia[row]-1;strt=ia[row-1]-1;elems=fin-strt;
for(i=0;i<elems;i++){
if(ja[i+strt]==col){a[i+strt]=1.0;} else {a[i+strt]=0.0;}
}
}

@<pathNewt signature@>
{
@<pathNewt declarations@>
@<pathNewt initializations@>
*check=0;done=0;
for(ihomotopy=0;(ihomotopy<numberAlphas)&&(!(*check));ihomotopy++){
#ifdef DEBUG 
printf("%d-th homotopy -- %f\n",ihomotopy,*(homotopyAlpha+ihomotopy));
#endif
@<q terminal constraint computation@>

its=0;
@<pathNewt check for convergence @>
	for (its=1;(its<=(maxitsInput))&&(!done);its++) {
if(*check == 0){
@<pathNewt update path@>
@<q terminal constraint computation@>
@<another pn check for convergence@>
	}
}

if(!done){
*check=1;
addOneToFailedQ;
 printf("MAXITS=%d exceeded in pathNewt\n",(maxitsInput));}
done=0;
}
if(ignoreFailQ){
if(*check==1){/*subOneFromFailedQ;*/
printf("IGNORING FAILED CONVERGENCE!!!!!!!!\n");
}
*check=0;

}
PFREERETURN
}
@}
@o stackC.h -d
@{
@<homotopyPathNewt signature@>;/*first one*/
@}
@d delete homotopyPathNewt signature
@{
void homotopyPathNewt(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
unsigned int * pathLength,
void (* vecfunc)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
void (* fdjac),
double * params,
double x[],double * linearizationPoint,unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
unsigned int *check,
unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo)
@}
@o stackC.h -d
@{
@<FPnewt signature@>;/*second one*/
@}

@d FPnewt declarations
@{
    unsigned int n/*;double * xdel;double  dignore[1]={1.0}*/;unsigned int done;
/*	double fmin(double x[]);*/unsigned int ihomotopy=1;double realAlpha;
/*	void lnsrch(unsigned int n,unsigned int np,unsigned int reps,
    double xold[], double   * fold, double g[], double p[],double * params,
		 double * shockVec,double * f, double stpmax, unsigned int *check, 
         void (*func)(double*,double*,double*,double*,unsigned int*,unsigned int*,double*,double*),double *x,
            unsigned int ihomotopy,double * linPt,unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo);*/
	void lubksb(double **a, unsigned int n, unsigned int *indx, double b[]);
	void ludcmp(double **a, unsigned int n, unsigned int *indx, double *d);
	unsigned int i,its/*,j*/,*indx,*aOne/*,*ndns*/,*ierr;
	double /*d,den,*/f,fold,stpmax,sum,temp,test,*g,*p,*xold,normSum;
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
	shockVec=(double *)calloc(*numberOfEquations,sizeof(double));
    for(i=0;i<*numberOfEquations;i++)shockVec[i]=0.0;
	nn=n;
/*save alphas but force fp to use orginal function*/
realAlpha=*homotopyAlpha;
/**homotopyAlpha=1;;*/
resetFailedQ;
resetNewtonSteps;
resetRealizedTolf;
resetRealizedTolx;
resetFEvals;
resetFDrvEvals;
resetShrinkSteps;
resetExpandSteps;
resetLnsrchSteps;

@}


@d evaluate func and find max element
@{


/* modification begin */
      func(x,params,shockVec,fmats[0],fmatsj[0],fmatsi[0],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
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
/*	if (test<0.01*(tolfInput)) FREERETURN*/
	for (sum=0.0,i=0;i<*numberOfEquations;i++) sum += SQR(x[i]);
	stpmax=(*lags+*leads+1)*STPMX*FMAX(sqrt(sum),(double)n);


        sparseMatTimesVec(numberOfEquations,
        smats[0],smatsj[0],smatsi[0],fvec,
        g+(*numberOfEquations * *lags));
        for(i=0;i>*numberOfEquations;i++){
          g[(*numberOfEquations* (*lags-1))+i]=fvec[(*numberOfEquations* *lags)+i];
          g[(*numberOfEquations* (*lags+1))+i]= -fvec[(*numberOfEquations* *lags)+i];
          }

@}

@d get newton update
@{
dfunc(x,params,shockVec,smats[0],smatsj[0],smatsi[0],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
		for (i=0;i<n;i++) xold[i]=x[i];
		/*modification begin*/
nxtFPGuess(numberOfEquations,lags,leads,
fmats[0],fmatsj[0],fmatsi[0],
smats[0],smatsj[0],smatsi[0],
maxNumberElements,/*x,*/p);




		lnsrch(n,*numberOfEquations,1,
        xold,&fold,g,p,params,shockVec,&f,stpmax,check,func,x,ihomotopy,
linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
/*intControlParameters,*/doubleControlParameters/*,
intOutputInfo, doubleOutputInfo*/);


		fold=f;


@}

@d evaluate func update norm
@{

      func(x,params,shockVec,fmats[0],fmatsj[0],fmatsi[0],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
        for (normSum=0.0,i=0;i<=fmatsi[0][*numberOfEquations]+fmatsi[0][0];i++) 
            normSum += SQR(fmats[0][i]);
        f= 0.5 * normSum * (*lags+*leads+1);
        csrToDns(numberOfEquations,
        aOne,fmats[0],fmatsj[0],fmatsi[0],fvec+(*numberOfEquations**lags),
        numberOfEquations,ierr);
pathNewtAssert(*ierr == 0);
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
		if (test < (tolfInput)) {
resetFailedQ;
			*check=0;
done=1;
		} else {
		test=0.0;
		for (i=0;i<n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < tolxInput) {
*check=0;done=1;}
}
@}


\subsection{pathNewt}
\label{sec:pathNewt} 




@d delete pathNewt signature
@{
void pathNewt(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,unsigned int * pathLength,
void (* vecfunc)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
void (* fdjac)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * params,double * shockVec,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,double * linearizationPoint,unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
double x[],
unsigned int *check, double * lastDel,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * ma50bdJob,
unsigned int * ma50bdIq,
double * ma50bdFact,
unsigned int * ma50bdIrnf,
unsigned int * ma50bdIptrl,
unsigned int * ma50bdIptru
)
@}

@o stackC.h -d 
@{
@<pathNewt signature@>;
@}


@o stackC.h -d
@{
@<homotopyPathNewt signature@>;/*second one*/
@}
@d homotopyPathNewt signature
@{
void homotopyPathNewt(unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
unsigned int * pathLength,
void (* vecfunc)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
void (* fdjac)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * params,double * shockVec,
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,double * qMat,unsigned int * qMatj,unsigned int * qMati,
double * fixedPoint,double * intercept,double * linearizationPoint,
/*unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,*/
double easyX[],double targetX[],unsigned int * exogQ,double x[],
unsigned int *check,unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo,
unsigned int * ma50bdJob,
/*unsigned int * ma50bdIq,*/
double * ma50bdFact,
unsigned int * ma50bdIrnf,
unsigned int * ma50bdIptrl,
unsigned int * ma50bdIptru
)
@}
@o myNewt.c -d
@{

 void putBadHomotopy(
 unsigned int  numberOfEquations,unsigned int  lags,unsigned int leads,unsigned int pathLength,
 double * easyX,double * targetX,double * x,double * lastDel,double * fixedPt,
					double * shock, FILE * hFile);


@<homotopyPathNewt signature@>
{

FILE * debFile;
FILE * hFile;
unsigned int lclMaxElems=0;
unsigned int originalMaxElems;
unsigned int i;
unsigned int j;
unsigned int lagPart;
unsigned int vecLength;
double * zeroShock;
double * trialShock;
double * trialX;
double * lastDel;
originalMaxElems=*maxNumberElements;
lagPart=*numberOfEquations* *lags;
vecLength=*numberOfEquations*( *lags + *leads + *pathLength);
lastDel = (double *)calloc(vecLength,sizeof(double));
trialX = (double *)calloc(vecLength,sizeof(double));
zeroShock = (double *)calloc(vecLength,sizeof(double));
trialShock = (double *)calloc(vecLength,sizeof(double));
addOneToHomotopies;
*check=0;
for(i=0;i<numberBetas;i++){
/* do convex linear combination for shock and for x */
/*needs to keep endog at last value but update exog using a template*/
/* call pathNewt */

/*set endog values last x*/
switch(homotopyXGuess){
case useBigX:
for(j=0;j<vecLength;j++){
if(j<lagPart||exogQ[j%*numberOfEquations]){
trialX[j]= (1-homotopyBeta[i])*easyX[j]  +  homotopyBeta[i]*targetX[j];}
 else {trialX[j]= x[j];}
}
break;
case useBigEasy:
for(j=0;j<vecLength;j++){
if(j<lagPart||exogQ[j%*numberOfEquations]){
trialX[j]= (1-homotopyBeta[i])*easyX[j]  +  homotopyBeta[i]*targetX[j];}
 else {trialX[j]= x[j];}
}
}



if(*check ==0){
#ifdef DEBUG 
printf("for homotopyBeta[%d]=%f\n",i,homotopyBeta[i]);
#endif
for(j=0;j<*numberOfEquations;j++){
trialShock[j]= homotopyBeta[i]*shockVec[j];
}
*maxNumberElements=originalMaxElems;
pathNewt(numberOfEquations,lags, leads,
pathLength,
vecfunc,fdjac,params,trialShock,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,
fixedPoint,intercept,linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
trialX,
check,lastDel,intControlParameters,doubleControlParameters,
intOutputInfo,doubleOutputInfo,
ma50bdJob,
/*ma50bdIq,*/
ma50bdFact,
ma50bdIrnf,
ma50bdIptrl,
ma50bdIptru
);
if(*check==0)for(j=0;j<vecLength;j++){x[j]=trialX[j];} else {
printf("already failed earlier homotopy attempt, ignoring %e\n", homotopyBeta[i]);}
/*set endog values last x*/
switch(homotopyXGuess){
case useBigX:
for(j=0;j<vecLength;j++){
if(/*j<lagPart||*/exogQ[j%*numberOfEquations]){
trialX[j]= (1-homotopyBeta[i])*easyX[j]  +  homotopyBeta[i]*targetX[j];}
 else {trialX[j]= x[j];}
}
break;
case useBigEasy:
for(j=0;j<vecLength;j++){
trialX[j]= (1-homotopyBeta[i])*easyX[j]  +  homotopyBeta[i]*targetX[j];}

}

}
if(*check){
addOneToHomotopyFailures;
debFile=fopen("codeGenDebHFile","w");
fprintf(debFile,"homotopy failing\n");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
targetX,"hom$fail$targetX");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
easyX,"hom$fail$easyX");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
x,"hom$fail$x");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + *pathLength),
lastDel,"hom$fail$xdel");
fPrintMathDbl(debFile,*numberOfEquations* (*lags + *leads + 1),
targetX,"hom$fail$fixedPt");
fPrintMathDbl(debFile,*numberOfEquations,
trialShock,(const char *)"hom$fail$shock");
fclose(debFile);
}
if(lclMaxElems<*maxNumberElements)lclMaxElems=*maxNumberElements;

}

if(*check){
printf("homotopy failing, writing to codeGenHFile\n");
hFile=fopen("codeGenHFile","w");
putBadHomotopy(*numberOfEquations,*lags,*leads,*pathLength,
easyX,targetX,x,lastDel,fixedPoint,shockVec,hFile);
fclose(hFile);

}

*maxNumberElements=lclMaxElems;
free(lastDel);
free(trialX);
free(zeroShock);
free(trialShock);

}

@}


@d q terminal constraint computation
@{

@<compute fmats for path@>

@<apply q terminal constraint@>
{/*use shock first period only*/
addOneToFDrvEvals;
      fdjac(x,params,shockVec,
      smats[0],smatsj[0],smatsi[0],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
}
for(tNow=1;tNow<*pathLength;tNow++) {
addOneToFDrvEvals;
      fdjac(x+(tNow* *numberOfEquations),params,zeroShockVec,
      smats[tNow],smatsj[tNow],smatsi[tNow],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
}


@}





@d compute fmats for path
@{
normSum=0.0;
{/*use shock for first period only*/
addOneToFEvals;
      vecfunc(x,params,shockVec,
      fmats[0],fmatsj[0],fmatsi[0],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
        for (i=0;
        i<=fmatsi[0][*numberOfEquations]-fmatsi[0][0];i++)
            normSum += SQR(fmats[0][i]);
}

for(tNow=1;tNow<*pathLength;tNow++) {
addOneToFEvals;
      vecfunc(x+(tNow* *numberOfEquations),params,zeroShockVec,
      fmats[tNow],fmatsj[tNow],fmatsi[tNow],homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/);
        for (i=0;
        i<=fmatsi[tNow][*numberOfEquations]-fmatsi[tNow][0];i++)
            normSum += SQR(fmats[0][i]);
}


@}


@d apply q terminal constraint
@{
if(*leads>0){
/*
copyMatrix(rowDim,aOne,qMat,qMatj,qMati,
aOne,
smats[*pathLength],smatsj[*pathLength],smatsi[*pathLength]);
*/

/*xxxxxxxxx add code for deviations using gmat*/
for(i=0;i<*numberOfEquations* (*lags+ *leads);i++){
deviations[i]=x[*numberOfEquations* *pathLength+i]-fixedPoint[i+*numberOfEquations];}
/*sparseMatTimesVec(rowDim,smats[*pathLength],smatsj[*pathLength],smatsi[*pathLength],deviations,fullfvec);
*/
sparseMatTimesVec(rowDim,qMat,qMatj,qMati,deviations,fullfvec);
for(i=0;i<*numberOfEquations* *leads;i++){fullfvec[i]=
fullfvec[i]-intercept[i];}
dnsToCsr(rowDim,aOne,rowDim,fullfvec,
aOne,
fmats[*pathLength],fmatsj[*pathLength],fmatsi[*pathLength],ierr);
pathNewtAssert(*ierr == 0);
bump(fmatsi[*pathLength][*rowDim]-fmatsi[*pathLength][0]);
for(i=0;i<*rowDim;i++){
  fvec[(*pathLength * *numberOfEquations)+i]=fullfvec[i];
            normSum += SQR(fullfvec[i]);}

        f= 0.5 * normSum;
#ifdef DEBUG 
/*printf("pathNewt:deviations computation\n");
printf("fullfvec={\n");
for(i=0;i<*rowDim-1;i++){printf("%f,",fullfvec[i]);}
printf("%f}\n",fullfvec[*rowDim]);
printf("deviations={\n");
for(i=0;i<*numberOfEquations* (*lags+ *leads)-1;i++){printf("%f,",deviations[i]);}
printf("%f}\n",deviations[*numberOfEquations* (*lags+ *leads)]);
*/
#endif

}
for(tNow=0;tNow<*pathLength;tNow++) {
        csrToDns(numberOfEquations,
        aOne,fmats[tNow],fmatsj[tNow],fmatsi[tNow],
        fvec+((*lags+tNow) * *numberOfEquations),numberOfEquations,ierr);
pathNewtAssert(*ierr == 0);
bump(*numberOfEquations);
        }
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
*qRows=0;
maxHElementsForSparseAim=maxHElements;
(void)sparseAMA(&maxHElementsForSparseAim,DISCRETE_TIME,*numberOfEquations,
*numberOfEquations*(*lags+1+*leads),*leads,
smats[0],smatsj[0],smatsi[0],
newH,newHj,newHi,
auxInit,qRows,
qMat,qMatj,qMati,
&essential,rootr,rooti,&returnCode);
*maxNumberElements=maxHElementsForSparseAim;
@|NEGLIGIBLEDOUBLE AIMROOTBOUND
@}

\subsection{computeAsymptoticQMatrix}
\label{sec:computeAsymptoticQMatrix}


@d create asymptotic linearization for Anderson-Moore algorithm
@{
dfunc(genericModelFP,params,shockVec,
smats[0],smatsj[0],smatsi[0],homotopyAlpha+ihomotopy,genericModelFP/*,exogRows,exogCols,exogenizeQ*/);






/**hColumns=(*lags+*leads+1)* *numberOfEquations;*/

@}


@o myNewt.c -d
@{

void putQ(unsigned int aux,unsigned int qrows,double * q,unsigned int * jq,unsigned int * iq,
unsigned int nroots,double * rootr,double * rooti,FILE * qFile){

  unsigned int ii;
  fprintf(qFile,"%10d %10d\n",aux,qrows);
  for(ii=0;ii<qrows+1;ii++){
	fprintf(qFile,"%12d ",iq[ii]);}
  for(ii=0;ii<iq[qrows]-iq[0];ii++){
	fprintf(qFile,"%12d ",jq[ii]);}
  for(ii=0;ii<iq[qrows]-iq[0];ii++){
	fprintf(qFile,"%30.15e ",q[ii]);}
  fprintf(qFile,"%10d \n",nroots);
  for(ii=0;ii<nroots;ii++){
	fprintf(qFile,"%30.15e ",rootr[ii]);}
  for(ii=0;ii<nroots;ii++){
	fprintf(qFile,"%30.15e ",rooti[ii]);}

}
void getQ(unsigned int * aux,unsigned int * qrows,double * q,unsigned int * jq,unsigned int * iq,
unsigned int * nroots,double * rootr,double * rooti,FILE * qFile){

  unsigned int ii;
  fscanf(qFile,"%10u %10u",aux,qrows);
  for(ii=0;ii<*qrows+1;ii++){
	fscanf(qFile,"%12u ",jq+ii);}
  for(ii=0;ii<iq[*qrows]-iq[0];ii++){
	fscanf(qFile,"%le ",q+ii);}
  fscanf(qFile,"%10u",nroots);
  for(ii=0;ii<*nroots;ii++){
	fscanf(qFile,"%le ",rootr+ii);}
  for(ii=0;ii<*nroots;ii++){
	fscanf(qFile,"%le ",rooti+ii);}
}
/*vecLength corresponds to actual number of shocks -- typically less than the numberOfEquations. routine pads vector with zeros to make the shockVec numberOfEquations long  numShocks is how many shocks vec to get*/
@<getShocks signature@>
{

  unsigned int ii=0;unsigned int jj=0;unsigned int kk=0;unsigned int ll=0;double ignore[1];
ll=lowLim;
for(ii=0;ii<lowLim*vecLength;ii++){
fscanf(sFile,"%le",ignore);
}
ii=0;
while((ii <(numShocks)*numberOfEquations)&&(ll<=upLim)){
while((jj<vecLength)&&(fscanf(sFile,"%le",shockArray+ii) != EOF))
{ii++;jj++;}
jj=0;ll++;
for(kk=0;kk<numberOfEquations-vecLength;kk++){shockArray[ii]=0;ii++;}
}
}
@}
@o stackC.h -d
@(
@<getShocks signature@>;
@<getData signature@>;

@}

@d getShocks signature
@{
void getShocks(unsigned int  numShocks,unsigned int  vecLength,unsigned int numberOfEquations,
double * shockArray,FILE * sFile,unsigned int lowLim,unsigned int upLim)
@}
@d getData signature
@{
void getData(unsigned int  numPeriods,unsigned int  vecLength,unsigned int numberOfEquations,
double * shockArray,FILE * sFile,unsigned int lowLim,unsigned int upLim)
@}
@d putData signature
@{
void putData(unsigned int numberOfEquations,double * q,FILE * qFile)
@}
@o stackC.h -d
@{
@<putData signature@>;
@}
@d computeAsymptoticQMatrix signature
@{
void computeAsymptoticQMatrix(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
/*void (*func)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),*/
void (*dfunc)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * params,double * shockVec,
double genericModelFP[],/*unsigned int *exogRows,unsigned int *exogCols,unsigned int *exogenizeQ,unsigned int * pthLngth,*/
/*double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,*/
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
/*unsigned int * maxNumberElements,*/
/*double * AMqMatrix,*/
/*unsigned int * ierr,*/unsigned int ihomotopy,
/*unsigned int * intControlParameters,*/double * doubleControlParameters/*,
unsigned int * intOutputInfo, double * doubleOutputInfo*/
)
@}

@o myNewt.c -d
@{

 

@<putData signature@>
{

  unsigned int ii;
  for(ii=0;ii<numberOfEquations;ii++){
	fprintf(qFile,"%30.15e ",q[ii]);}
}

@<getData signature@>{

  unsigned int ii=0;unsigned int jj=0;unsigned int kk=0;unsigned int ll=0;double ignore[1];
ll=lowLim;
for(ii=0;ii<lowLim*vecLength;ii++){
fscanf(sFile,"%le",ignore);
}
ii=0;
while((ii <(numPeriods)*numberOfEquations)&&(ll<=upLim)){
while((jj<vecLength)&&(fscanf(sFile,"%le",shockArray+ii) != EOF))
{ii++;jj++;}
jj=0;ll++;
for(kk=0;kk<numberOfEquations-vecLength;kk++){shockArray[ii]=0;ii++;}
}
}



void putBadHomotopy(
unsigned int  numberOfEquations,unsigned int  lags,unsigned int leads,unsigned int pathLength,
double * easyX,double * targetX,double * x,double * lastDel,double * fixedPt,
					double * shock, FILE * hFile)
{

unsigned int ii;unsigned int bigVecLength;
bigVecLength=numberOfEquations*(lags+leads+pathLength);
  fprintf(hFile,"%10d %10d %10d %10d\n",
  numberOfEquations,lags,leads,pathLength);
  for(ii=0;ii<bigVecLength;ii++){
	fprintf(hFile,"%30.15e ",easyX[ii]);}
  for(ii=0;ii<bigVecLength;ii++){
	fprintf(hFile,"%30.15e ",targetX[ii]);}
  for(ii=0;ii<bigVecLength;ii++){
	fprintf(hFile,"%30.15e ",x[ii]);}
  for(ii=0;ii<bigVecLength;ii++){
	fprintf(hFile,"%30.15e ",lastDel[ii]);}
  for(ii=0;ii<numberOfEquations*(lags+leads+1);ii++){
	fprintf(hFile,"%30.15e ",fixedPt[ii]);}
  for(ii=0;ii<numberOfEquations;ii++){
	fprintf(hFile,"%30.15e ",shock[ii]);}
}

void getBadHomotopy(
unsigned int  * numberOfEquations,unsigned int  * lags,unsigned int * leads,unsigned int * pathLength,
double * easyX,double * targetX,double * x,double * lastDel,double * fixedPt,
					double * shock, FILE * hFile)
{

unsigned int ii;unsigned int bigVecLength;
  fscanf(hFile,"%10u %10u %10u %10u",
  numberOfEquations,lags,leads,pathLength);
bigVecLength=*numberOfEquations*(*lags+*leads+*pathLength);
  for(ii=0;ii<bigVecLength;ii++){
	fscanf(hFile,"%le ",easyX+ii);}
  for(ii=0;ii<bigVecLength;ii++){
	fscanf(hFile,"%le ",targetX+ii);}
  for(ii=0;ii<bigVecLength;ii++){
	fscanf(hFile,"%le ",x+ii);}
  for(ii=0;ii<bigVecLength;ii++){
	fscanf(hFile,"%le ",lastDel+ii);}
  for(ii=0;ii<*numberOfEquations*(*lags+*leads+1);ii++){
	fscanf(hFile,"%le ",fixedPt+ii);}
  for(ii=0;ii<*numberOfEquations;ii++){
	fscanf(hFile,"%le ",shock+ii);}
}




@<computeAsymptoticQMatrix signature@>
{
@<computeAsymptoticQMatrix variable declarations@>
@<computeAsymptoticQMatrix variable allocations@>
//  *ierr=0;
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
double * asymptoticLinearization;//double  dignore[1]={1.0};
double * cond;
double * epsi;
unsigned int * inform;
unsigned int * iq;
unsigned int * itsbad;
unsigned int * nbig;
unsigned int * nexa;
unsigned int * nnum;
unsigned int * nroot;
//double * zeroVector;
unsigned int * qColumns;
double * uprbnd;
//double * wts;
//double * err;
unsigned int * hColumns;
//unsigned int i;
@}
@d allocFPNewt signature
@{
void allocFPNewt(unsigned int numberOfEquations,unsigned int lags,unsigned int leads,
unsigned int pathLength,unsigned int maxElements,
double ** genericFP,
double ** genericIntercept,
double***fmats,unsigned int***fmatsj,unsigned int***fmatsi,
double***smats,unsigned int***smatsj,unsigned int***smatsi)
@}

@o stackC.h -d
@{
@<allocFPNewt signature@>;
@<altComputeAsymptoticQMatrix signature@>;
@}
@o myNewt.c -d
@{
@<altComputeAsymptoticIMatrix signature@>
{
unsigned int j;
/*construct identity matrices*/
for(j=0;j<*numberOfEquations* *leads;j++){
qMat[j]=1;
qMatj[j]=(*numberOfEquations * (*lags))+j+1;
qMati[j]=j+1;
}
qMati[*numberOfEquations* *leads]= (*numberOfEquations * *leads)+1;
}


void altComputeAsymptoticDMatrix(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
double * qMat,unsigned int * qMatj,unsigned int * qMati/*,
unsigned int * ierr*/
)
{
unsigned int j;
/*construct identity matrices*/
for(j=0;j<*numberOfEquations* *leads;j++){
qMat[2*j]=1;
qMatj[2*j]=(*numberOfEquations * (*lags))+j+1;
qMat[(2*j)+1]=-1;
qMatj[(2*j)+1]=(*numberOfEquations * (*lags-1))+j+1;
qMati[j]=(2*j)+1;
}
qMati[*numberOfEquations* *leads]= 2*(*numberOfEquations * *leads)+1;
}

@<altComputeAsymptoticQMatrix signature@>
{
@<altComputeAsymptoticQMatrix variable declarations@>
@<altComputeAsymptoticQMatrix variable allocations@>
  *ierr=0;

@<create asymptotic linearization for Anderson-Moore algorithm@>
@<alt apply the Anderson-Moore algorithm  to get AMqMatrix@>

@<altComputeAsymptoticQMatrix variable deallocations@>
}
@}
@d altComputeAsymptoticQMatrix signature
@{
void altComputeAsymptoticQMatrix(
unsigned int * numberOfEquations,unsigned int * lags, unsigned int * leads,
/*void (*func)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),*/
void (*dfunc)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
double * params,double * shockVec,
double genericModelFP[],/*unsigned int *exogRows,unsigned int *exogCols,unsigned int *exogenizeQ,*//*unsigned int * pthLngth,*//*
double ** fmats, unsigned int ** fmatsj, unsigned int ** fmatsi,*/
double ** smats, unsigned int ** smatsj, unsigned int ** smatsi,
unsigned int * maxNumberElements,
double * qMat,unsigned int * qMatj,unsigned int * qMati,unsigned int * auxInit,unsigned int * qRows,
double * rootr, double * rooti,
unsigned int * ierr,unsigned int ihomotopy/*,*/
/*unsigned int * intControlParameters,double * doubleControlParameters,*/
/*unsigned int * intOutputInfo, double * doubleOutputInfo*/
)
@}

@d altComputeAsymptoticQMatrix variable deallocations
@{
free(newH);
free(newHj);
free(newHi);
@}
@d altComputeAsymptoticQMatrix variable allocations
@{
maxHElements=*maxNumberElements;
newH=(double *)calloc(maxHElements,sizeof(double));
newHj=(unsigned int *)calloc(maxHElements,sizeof(unsigned int));
newHi=(unsigned int *)calloc(maxHElements,sizeof(unsigned int));

@}
@d altComputeAsymptoticQMatrix variable declarations
@{
/*void * hookForPeter;*//*troll needs to have this for their calloc in sparseAMA*/
//double  dignore[1]={1.0};
unsigned int maxHElements=50000;
unsigned int maxHElementsForSparseAim;
unsigned int essential;
unsigned int returnCode;
//double * zeroVector;
//unsigned int * qColumns;
double * newH;
unsigned int * newHj;
unsigned int * newHi;
//double * uprbnd;
//double * wts;
//double * err;
//unsigned int * hColumns;
//unsigned int i;
@}


\subsection{Prepare Linear Terminator}
\label{sec:lunsigned interm}
@d linearTerminator definition
@{
void linearTerminator(
@<linearTerminator argument list@>
)
{
@<linearTerminator declarations@>
@<linearTerminator allocations@>

/*compute linearization*/
dfunc(linearizationPoint,params,shockVec,hMat,hMatj,hMati,homotopyAlpha+ihomotopy,linearizationPoint,exogRows,exogCols,exogenizeQ);

/*apply anderson-moore algorithm*/
maxHElementsForSparseAim=*maxNumberElements;
sparseAMA(&maxHElementsForSparseAim,DISCRETE_TIME,*numberOfEquations,
*numberOfEquations*(*lags+1+*leads),*leads,
hMat,hMatj,hMati,
hzMat,hzMatj,hzMati,/*will be overwritten by subsequent call to get hzMat*/
auxInit,qRows,
qMat,qMatj,qMati,
&essential,rootr,rooti,ierr,hookForPeter);


/*compute reduced form b*/
qCols=*numberOfEquations*(*lags+*leads);
obtainSparseReducedForm(maxNumberElements, 
*qRows,qCols,qMat,qMatj,qMati,
bMat,bMatj,bMati);

/*compute phiInv*/
/*sparseAMA destroys hMat*/
dfunc(linearizationPoint,params,shockVec,hMat,hMatj,hMati,homotopyAlpha+ihomotopy,linearizationPoint,exogRows,exogCols,exogenizeQ);
computePhiInv(maxNumberElements, 
*numberOfEquations,*numberOfEquations*(*lags+*leads+1),
hMat,hMatj,hMati,
*numberOfEquations*(*leads),*numberOfEquations*(*lags),
bMat,bMatj,bMati,
phiInvMat,phiInvMatj,phiInvMati
);

/*compute fmat varthetac varthetazstar*/

exfunc(params,linearizationPoint,hzMat,hzMatj,hzMati);

computeF(maxNumberElements, 
*numberOfEquations,*numberOfEquations*(*lags+*leads+1),
hMat,hMatj,hMati,
hzMat,hzMatj,hzMati,
*numberOfEquations*(*leads),*numberOfEquations*(*lags),
bMat,bMatj,bMati,
phiInvMat,phiInvMatj,phiInvMati,
*numberExogenous,
psimat,psimatj,psimati,
upsilonMat,upsilonMatj,upsilonMati,
fmat,fmatj,fmati,
varthetaZstar,varthetaZstarj,varthetaZstari,
impact,impactj,impacti,
varthetaC,varthetaCj,varthetaCi
);



/*func(linearizationPoint,params,zeroShock,cstar,cstarj,cstari,homotopyAlpha+ihomotopy,linearizationPoint,exogRows,exogCols,exogenizeQ);*/


@<linearTerminator frees@>

}
@}

@d linearTerminator declarations
@{
unsigned int i;double  dignore[1]={1.0};
double * zeroShock;unsigned int ihomotopy=0;
unsigned int  essential;
unsigned int maxHElementsForSparseAim;
unsigned int qCols;
double * psimat;
unsigned int * psimatj;
unsigned int * psimati;

@}

@d linearTerminator allocations
@{
zeroShock=(double *)calloc(*numberOfEquations,sizeof(double));
psimat=(double *)calloc(*maxNumberElements,sizeof(double));
psimatj=(unsigned int *)calloc(*maxNumberElements,sizeof(unsigned int));
psimati=(unsigned int *)calloc(*numberOfEquations+1,sizeof(unsigned int));
@}

@d linearTerminator frees
@{
free(zeroShock);
free(psimat);
free(psimatj);
free(psimati);
@}

@d linearTerminator argument list
@{
unsigned int * numberOfEquations,
unsigned int * lags, unsigned int * leads,unsigned int * numberExogenous,
double * linearizationPoint,unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
void (*func)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
void (*dfunc)(double *,double*, double*, double*,unsigned  int*,unsigned  int*,double *, double *),
void(*exfunc)(),double * params,double * shockVec,
double * upsilonMat, unsigned int * upsilonMatj, unsigned int * upsilonMati,
unsigned int * maxNumberElements,
double * hMat,unsigned int * hMatj,unsigned int * hMati,
double * hzMat,unsigned int * hzMatj,unsigned int * hzMati,
double * qMat,unsigned int * qMatj,unsigned int * qMati,
unsigned int * auxInit,unsigned int * qRows,
double * rootr, double * rooti,
double * bMat,unsigned int * bMatj,unsigned int * bMati,
double * phiInvMat,unsigned int * phiInvMatj,unsigned int * phiInvMati,
double * fmat, unsigned int * fmatj, unsigned int * fmati,
double * varthetaZstar, unsigned int * varthetaZstarj, unsigned int * varthetaZstari,
double * impact, unsigned int * impactj, unsigned int * impacti,
double * varthetaC, unsigned int * varthetaCj, unsigned int * varthetaCi,
unsigned int * ierr,
unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo
@}
\subsection{computeF}
\label{sec:computeF}

@d kron definition
@{
void kron(
unsigned int leftRows,unsigned int leftCols,
double * leftMat, unsigned int * leftMatj, unsigned int * leftMati,
unsigned int rightRows,unsigned int rightCols,
double * rightMat, unsigned int * rightMatj, unsigned int * rightMati,
unsigned int * maxNumberElements,
unsigned int * resultRows,unsigned int  * resultCols,
double * resultMat, unsigned int * resultMatj, unsigned int * resultMati)
{
unsigned int i; unsigned int j;unsigned int elements;unsigned int soFar=0;
unsigned int leftRowNow; unsigned int rightRowNow;
unsigned int leftColNow; unsigned int rightColNow;
/*check if enough space*/
resultMati[1]=resultMati[0]=1;
if((elements=
(leftMati[leftRows]-leftMati[0])*
(rightMati[rightRows]-rightMati[0]))>*maxNumberElements){ 
  *maxNumberElements=elements;
  *resultRows=0;*resultCols=0;} else {
/*for each row in left matrix*/
for(leftRowNow=0;leftRowNow<leftRows;leftRowNow++){
/*for each row in right matrix*/
for(rightRowNow=0;rightRowNow<rightRows;rightRowNow++){
/*for each element in row in left matrix*/
for(i=leftMati[leftRowNow]-1;i<leftMati[leftRowNow+1]-1;i++){
/*for each element in row in right matrix*/
for(j=rightMati[rightRowNow]-1;j<rightMati[rightRowNow+1]-1;j++){
resultMat[soFar]=leftMat[i]*rightMat[j];
resultMatj[soFar]=(leftMatj[i]-1)*rightCols+rightMatj[j];
soFar++;
resultMati[leftRowNow*rightRows+rightRowNow+1]=soFar+1;
}
resultMati[leftRowNow*rightRows+rightRowNow+1]=soFar+1;
}
resultMati[leftRowNow*rightRows+rightRowNow+1]=soFar+1;
}
resultMati[leftRowNow*rightRows+rightRowNow]=soFar+1;
}
resultMati[leftRows*rightRows]=soFar+1;
*resultRows=leftRows*rightRows;*resultCols=leftCols*rightCols;
}
}

@}
@d computeUnsigned Intercept definition
@{
void computeUnsigned Intercept(@<computeUnsigned Intercept argument list@>)
{
@<computeUnsigned Intercept declarations@>
@<computeUnsigned Intercept allocations@>
if(leads=1){
/*initialize unsigned intercept to zero*/
for(i=0;i<hrows*leads;i++){unsigned intercept[i]=0;}
@<compute zdeviation at time t@>
} else { printf("multiple leads not implemented yet !!!!!!!!!!!!!!!!!!!!\n");}
@<computeUnsigned Intercept frees@>
}
@}
@d compute zdeviation at time t
@{
        sparseMatTimesVec(&numberExogenous,
        selectZmat,selectZmatj,selectZmati,linearizationPt,zLin);
        sparseMatTimesVec(&numberExogenous,
        selectZmat,selectZmatj,selectZmati,newPt,zNew);
for(i=0;i<numberExogenous;i++){zDeviation[i]=zNew[i]-zLin[i];}
        sparseMatTimesVec(&hrows,
        impact,impactj,impacti,zDeviation,unsigned intercept);
csrToDns(&hrows,&aOne,cstar,cstarj,cstari,cXstar,&hrows,&ierr);
pathNewtAssert(ierr == 0);
bump(hrows);
        sparseMatTimesVec(&hrows,
        varthetaC,varthetaCj,varthetaCi,cXstar,unsigned intercept);

@}

@d computeUnsigned Intercept declarations
@{
unsigned int maxElementsEncountered=0;/*local to obtain Sparse still need to pass back to caller*/
unsigned int i;unsigned int ierr;unsigned int aOne=1;
double * zDeviation;
double * zLin;
double * zNew;
double * cXstar;
@}
@d computeUnsigned Intercept allocations
@{
zDeviation = calloc(numberExogenous,sizeof(double));
zLin = calloc(numberExogenous,sizeof(double));
zNew = calloc(numberExogenous,sizeof(double));
cXstar = calloc(hrows,sizeof(double));
@}
@d computeUnsigned Intercept frees
@{
free(zDeviation);
free(zLin);
free(zNew);
free(cXstar);
@}

@d computeUnsigned Intercept argument list
@{
unsigned int * maxNumberHElements, 
unsigned int hrows, unsigned int lags, unsigned int leads,
unsigned int numberExogenous,
 double * hzMat, unsigned int * hzMatj, unsigned int * hzMati,
 double * selectZmat, unsigned int * selectZmatj, unsigned int * selectZmati,
 double * vartheta, unsigned int * varthetaj, unsigned int * varthetai,
 double * impact, unsigned int * impactj, unsigned int * impacti,
 double * varthetaC, unsigned int * varthetaCj, unsigned int * varthetaCi,
double * cstar,unsigned int * cstarj,unsigned int * cstari,
double * linearizationPt,
double * newPt,
double * intercept
@}


@d computeF definition
@{
void kron(
unsigned int leftRows,unsigned int leftCols,
double * leftMat, unsigned int * leftMatj, unsigned int * leftMati,
unsigned int rightRows,unsigned int rightCols,
double * rightMat, unsigned int * rightMatj, unsigned int * rightMati,
unsigned int * maxNumberElements,
unsigned int * resultRows,unsigned int  * resultCols,
double * resultMat, unsigned int * resultMatj, unsigned int * resultMati);

computeF(@<computeF argument list@>)
{
@<computeF declarations@>
@<computeF allocations@>
@<factor phiInv@>
if(hcols-(bcols+2*hrows)>0){
printf("not implemented yet!!!!!!!!!!!!!!!!!!!!!!!!!!********************");
@<extract longb@>
} else {
@<compute factored phiInv hplus@>
@<compute factored phiInv imat@>
@<compute factored phiInv psi@>
@<transpose upsilon compute kron subtract@>
@<factor kron thing@>
@<compute vec vartheta@>
@<compute exogenous impact@>
@<imat compute kron subtract@>
@<factor kron thing@>
@<compute vec varthetaC@>
}

@<computeF deallocations@>
}
@}
@d compute factored phiInv psi
@{








for(i=0;i<numberExogenous;i++){

firstColumn=( numberExogenous*(bcols/hrows))+1+i;
lastColumn=( numberExogenous*(bcols/hrows))+1+i;
extractSubmatrix(&hrows,aOne,
aOne,&hrows,&firstColumn,&lastColumn,
hzMat,hzMatj,hzMati,
&nr,&nc,
psimat,psimatj,psimati);


csrToDns(&hrows,aOne,psimat,psimatj,psimati,
columnVals,&hrows,&ierr);
pathNewtAssert(ierr == 0);
bump(hrows);

    trans = 1;
useMA50CD(&hrows,&hrows,icntl,lclphiInvmati,np,&trans,
lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
columnVals,phipsi+(i*hrows),
w,info);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
/*change sign*/


}
for(i=0;i<hrows*numberExogenous;i++)phipsi[i]=(-1)*phipsi[i];/*since the kron part is negative*/
@}

@d compute factored phiInv imat
@{

/*set  columVals is vec of identitymatrix*/
for(i=0;i<hrows;i++){
for(j=0;j<hrows;j++){columnVals[j]=0;}
columnVals[i]=1;
    trans = 1;
useMA50CD(&hrows,&hrows,icntl,lclphiInvmati,np,&trans,
lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
columnVals,thePhi+(i*hrows),
w,info);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
/*change sign*/
}

for(i=0;i<hrows*hrows;i++)thePhi[i]=(-1)*thePhi[i];/*since the kron part is negative*/
@}

@d transpose upsilon compute kron subtract
@{
/*compute kron product*/
csrToCsc(&numberExogenous,aOne,aOne,upsilonmat,upsilonmatj,upsilonmati,
tupsilon,tupsilonj,tupsiloni);

kron(numberExogenous, 
numberExogenous,
tupsilon,tupsilonj,tupsiloni,
hrows,hrows,
fmat,fmatj,fmati,
maxNumberHElements,
&resultRows,&resultCols,
resultMat,resultMatj,resultMati);

/*compute negative of I - tup kron f*/
aplsca_(&resultRows,
resultMat,resultMatj,resultMati,
aMinusOne,iw);

@}

@d imat compute kron subtract
@{

/*set tupsilon to identitymatrix*/

for(i=0;i<hrows;i++){
tupsilon[i]=1;
tupsiloni[i]=tupsilonj[i]=i+1;
}
tupsiloni[hrows]=hrows+1;



/*compute kron product*/
kron(hrows,
hrows,
tupsilon,tupsilonj,tupsiloni,
hrows,hrows,
fmat,fmatj,fmati,
maxNumberHElements,
&resultRows,&resultCols,
resultMat,resultMatj,resultMati);

/*compute negative of I - tup kron f*/
aplsca_(&resultRows,
resultMat,resultMatj,resultMati,
aMinusOne,iw);

@}
@d computeF argument list
@{
unsigned int * maxNumberHElements, 
unsigned int hrows, unsigned int hcols,
 double * hmat, unsigned int * hmatj, unsigned int * hmati,
 double * hzMat, unsigned int * hzMatj, unsigned int * hzMati,
unsigned int brows, unsigned int bcols,
 double * bmat, unsigned int * bmatj, unsigned int * bmati,
 double * phiInvmat, unsigned int * phiInvmatj, unsigned int * phiInvmati,
unsigned int numberExogenous,
 double * psimat, unsigned int * psimatj, unsigned int * psimati,
 double * upsilonmat, unsigned int * upsilonmatj, unsigned int * upsilonmati,
 double * fmat,unsigned int * fmatj, unsigned int * fmati,
 double * vartheta, unsigned int * varthetaj, unsigned int * varthetai,
 double * impact, unsigned int * impactj, unsigned int * impacti,
 double * varthetaC, unsigned int * varthetaCj, unsigned int * varthetaCi
@}

@d computeF declarations
@{
unsigned int maxElementsEncountered=0;/*local to obtain Sparse still need to pass back to caller*/
 double * tupsilon; unsigned int * tupsilonj; unsigned int * tupsiloni;
 double * tfmat; unsigned int * tfmatj; unsigned int * tfmati;
unsigned int resultRows;unsigned int resultCols;
double * resultMat;
unsigned int * resultMatj;
unsigned int * resultMati;
unsigned int nzmaxLeft;
unsigned int lastColumn;unsigned int firstColumn;unsigned int nr;unsigned int nc;unsigned int ierr;
double * columnVals;
double * columnResult;
double * thePhi;
double * phipsi;
double * longb;
unsigned int * longbj;
unsigned int * longbi;
unsigned int i;
unsigned int j;
unsigned int soFar=0;
@}
@d computeF allocations
@{
impactRows=resultRows= (hcols+numberExogenous)*(hcols+numberExogenous);
columnVals = (double *)calloc(hrows*(hrows+numberExogenous+1),sizeof(double));
columnResult = (double *)calloc(resultRows,sizeof(double));
resultMat = (double *)calloc(*maxNumberHElements,sizeof(double));
resultMatj = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
resultMati = (unsigned int *)calloc((hrows+numberExogenous)*hrows+1,sizeof(unsigned int));/*too big to accomodate both uses*/
thePhi = (double *)calloc(hrows*(hrows+numberExogenous),sizeof(double));
phipsi = (double *)calloc(hrows*(hrows+numberExogenous),sizeof(double));
longb = (double *)calloc(*maxNumberHElements,sizeof(double));
longbj = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
longbi = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
tfmat = (double *)calloc(*maxNumberHElements,sizeof(double));
tfmatj = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
tfmati = (unsigned int *)calloc(resultRows+1,sizeof(unsigned int));
tupsilon = (double *)calloc(*maxNumberHElements,sizeof(double));
tupsilonj = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
tupsiloni = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));

@}

@d computeF deallocations
@{
free(resultMat);
free(resultMatj);
free(resultMati);
free(aMinusOne);
free(columnVals);
free(columnResult);
free(thePhi);
free(phipsi);
free(longb);
free(longbj);
free(longbi);
free(tfmat);
free(tfmatj);
free(tfmati);
free(tupsilon);
free(tupsilonj);
free(tupsiloni);
@}

@d computeF declarations
@{
unsigned int impactRows;
double aSmallDouble=DBL_EPSILON;
unsigned int * ma50bdIq;
double * ma50bdFact;
unsigned int * ma50bdIrnf;
unsigned int * ma50bdIptrl;
unsigned int * ma50bdIptru;
unsigned int * ma50bdJob;
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
unsigned int  * notAOne;
unsigned int  * aOne;
double  * aMinusOne;
unsigned int nonZeroNow;
unsigned int  trans;
double * lclphiInvmat;
unsigned int * lclphiInvmatj;
unsigned int * lclphiInvmati;
@}

@d computeF allocations
@{
ma50bdIptru = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
ma50bdIptrl = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
ma50bdIrnf = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
ma50bdFact = (double *)calloc(*maxNumberHElements,sizeof(double));
ma50bdIq = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
ma50bdJob = (unsigned int *)calloc(1,sizeof(unsigned int));
jcn = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
cntl= (double *)calloc(5,sizeof(double));
icntl= (unsigned int *)calloc(9,sizeof(unsigned int));
ip = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
np = (unsigned int *)calloc(1,sizeof(unsigned int));
jfirst = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
lenr = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
lastr = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
nextr = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
w = (double *)calloc(resultRows,sizeof(double));
iw = (unsigned int *)calloc(3*hrows*(hrows+numberExogenous+1),sizeof(unsigned int));/*extra space to cover both uses of iw*/
ifirst = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
lenc = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
lastc = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
nextc = (unsigned int *)calloc(resultRows,sizeof(unsigned int));
info = (unsigned int *)calloc(7,sizeof(unsigned int));
rinfo = (double *)calloc(1,sizeof(double));
lfact = (unsigned int *)calloc(1,sizeof(unsigned int));
*lfact = ( *maxNumberHElements);/*pessimistic setting for filling*/
aOne = (unsigned int *)calloc(1,sizeof(unsigned int));
*aOne=1;
notAOne = (unsigned int *)calloc(1,sizeof(unsigned int));
aMinusOne = (double *)calloc(1,sizeof(double));
*aMinusOne=-1.;
lclphiInvmat= (double * )calloc(*maxNumberHElements,sizeof(double));
lclphiInvmatj= (unsigned int * )calloc(*maxNumberHElements,sizeof(unsigned int));
lclphiInvmati= (unsigned int * )calloc((1+numberExogenous+hrows)*resultRows+1,sizeof(unsigned int));

@}
@d computeF deallocations
@{
free(ma50bdJob);
free(ma50bdIq);
free(ma50bdFact);
free(ma50bdIrnf);
free(ma50bdIptrl);
free(ma50bdIptru);

free(lclphiInvmat);
free(lclphiInvmatj);
free(lclphiInvmati);
free(notAOne);
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
@}

@d factor phiInv
@{


useMA50ID(cntl,icntl);
ma50bdJob[0] =1;


nzmax=*maxNumberHElements;
nonZeroNow=phiInvmati[hrows]-phiInvmati[0];

copyMatrix(&hrows,aOne,phiInvmat,phiInvmatj,phiInvmati,aOne,
lclphiInvmat,lclphiInvmatj,lclphiInvmati);

useMA50AD(&hrows,&hrows,&nonZeroNow,
&nzmax,lclphiInvmat,lclphiInvmatj,jcn,lclphiInvmati,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
if(*ma50bdJob!=2){
useMA50BD(&hrows,&hrows,&nonZeroNow,ma50bdJob,
phiInvmat,phiInvmatj,phiInvmati,
cntl,icntl,ip,lclphiInvmati,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
if(*ma50bdJob=1)*ma50bdJob=1;/* if it was 1 promote 
                             to 3(ie conservative alternative)*/
} else {
useMA50BD(&hrows,&hrows,&nonZeroNow,ma50bdJob,
phiInvmat,phiInvmatj,phiInvmati,
cntl,icntl,ip,lclphiInvmati,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
if(info[0]<-7){
/* small pivot values,
reset to 3*/
#ifdef DEBUG 
printf("small pivots, resetting to 3\n");
#endif
*ma50bdJob=1;
useMA50BD(&hrows,&hrows,&nonZeroNow,ma50bdJob,
phiInvmat,phiInvmatj,phiInvmati,
cntl,icntl,ip,lclphiInvmati,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
}
}
@}

@d factor kron thing
@{
useMA50ID(cntl,icntl);
ma50bdJob[0] =1;


nzmax=*maxNumberHElements;
nonZeroNow=resultMati[resultRows]-resultMati[0];

copyMatrix(&resultRows,aOne,resultMat,resultMatj,resultMati,aOne,
lclphiInvmat,lclphiInvmatj,lclphiInvmati);

useMA50AD(&resultRows,&resultRows,&nonZeroNow,
&nzmax,lclphiInvmat,lclphiInvmatj,jcn,lclphiInvmati,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
if(*ma50bdJob!=2){
useMA50BD(&resultRows,&resultRows,&nonZeroNow,ma50bdJob,
resultMat,resultMatj,resultMati,
cntl,icntl,ip,lclphiInvmati,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
if(*ma50bdJob=1)*ma50bdJob=1;/* if it was 1 promote 
                             to 3(ie conservative alternative)*/
} else {
useMA50BD(&resultRows,&resultRows,&nonZeroNow,ma50bdJob,
resultMat,resultMatj,resultMati,
cntl,icntl,ip,lclphiInvmati,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
if(info[0]<-7){
/* small pivot values,
reset to 3*/
#ifdef DEBUG 
printf("small pivots, resetting to 3\n");
#endif
*ma50bdJob=1;
useMA50BD(&resultRows,&resultRows,&nonZeroNow,ma50bdJob,
resultMat,resultMatj,resultMati,
cntl,icntl,ip,lclphiInvmati,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
}
}
@}
@d compute factored phiInv hplus
@{


for(i=0;i<hrows;i++){


lastColumn = firstColumn=(hcols-hrows+1+i);
for(j=0;j<hrows;j++){longbi[j]=1;}/*initialize to zero matrix*/
extractSubmatrix(&hrows,aOne,
aOne,&hrows,&firstColumn,&lastColumn,
hmat,hmatj,hmati,
&nr,&nc,
longb,longbj,longbi);

csrToDns(&hrows,aOne,longb,longbj,longbi,
columnVals,&hrows,&ierr);
pathNewtAssert(ierr == 0);
bump(longbi[hrows]-longbi[0]);
    trans = 1;
useMA50CD(&hrows,&hrows,icntl,lclphiInvmati,np,&trans,
lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
columnVals,columnResult,
w,info);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);

/*
useMA50CD(&hrows,&hrows,icntl,qrmati,np,&trans,
lfact,fact,irnf,iptrl,iptru,
nsSumC,x,
w,info);
wordybump((unsigned int)info[3]);
pathNewtAssert(info[0]>=0);
*/
nzmaxLeft= nzmax-soFar-1;
dnsToCsr(aOne,&hrows,&nzmaxLeft,columnResult,
aOne,tfmat+(tfmati[i]-1),tfmatj+(tfmati[i]-1),tfmati+i,
&ierr);
pathNewtAssert(ierr == 0);
tfmati[i+1]=tfmati[i+1]+soFar;
tfmati[i]=tfmati[i]+soFar;
soFar=tfmati[i+1]-1;
}


dropSmallElements(&hrows,aOne,&aSmallDouble,&nzmax,tfmat,tfmatj,tfmati,tfmat,tfmatj,tfmati,&ierr);
pathNewtAssert(ierr == 0);
csrcsc2_(&hrows,&hrows,aOne,aOne,tfmat,tfmatj,tfmati,fmat,fmatj,fmati);
/*change sign*/
for(i=0;i<fmati[hrows]-fmati[0];i++)fmat[i]=(-1)*fmat[i];

@}
@d compute vec vartheta
@{
    trans = 1;
useMA50CD(&resultRows,&resultRows,icntl,lclphiInvmati,np,&trans,
lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
phipsi,columnResult,
w,info);
pathNewtAssert(info[0]>=0);
dnsToCsr(&hrows,&numberExogenous,&nzmaxLeft,columnResult,
&hrows,vartheta,varthetaj,varthetai,
&ierr);
wordybump((unsigned int)info[3]);
pathNewtAssert(ierr == 0);
@}

@d compute vec varthetaC
@{
    trans = 1;
useMA50CD(&resultRows,&resultRows,icntl,lclphiInvmati,np,&trans,
lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
thePhi,columnResult,
w,info);
pathNewtAssert(info[0]>=0);
dnsToCsr(&hrows,&hrows,&nzmaxLeft,columnResult,
&hrows,varthetaC,varthetaCj,varthetaCi,
&ierr);
wordybump((unsigned int)info[3]);
pathNewtAssert(ierr == 0);
@}

@d compute exogenous impact
@{
copyMatrix(&hrows,aOne,vartheta,varthetaj,varthetai,aOne,
impact,impactj,impacti);
/*compute vartheta * upsilon*/
sparseMult(&hrows,&numberExogenous,aOne,vartheta,varthetaj,varthetai,
upsilonmat,upsilonmatj,upsilonmati,tfmat,tfmatj,tfmati,maxNumberHElements,iw,&ierr);
pathNewtAssert(ierr == 0);
*notAOne=impacti[hrows]-impacti[0]+1;
/*post mult by selector matrix*/

copyMatrix(&hrows,aOne,tfmat,tfmatj,tfmati,notAOne,
impact,impactj,impacti+hrows);

@}

@d extract longb
@{
for(i=0;i<hrows;i++){
longbi[i]=soFar+1;
longb[soFar]=1;
longbj[soFar]=i;
soFar++;
/*check to see if more than one lead if not were done with the row*/
}
@}

\subsection{computePhiInv}
\label{sec:computePhiInv}





@d computePhiInv definition
@{
computePhiInv(@<computePhiInv argument list@>)
{
@<computePhiInv declarations@>
@<computePhiInv allocations@>
@< extract hzero@>
@< extract hplus@>
@< extract bright@>
@< compute hplus times bright@>
@< add hzero@>
@<computePhiInv deallocations@>
#ifdef DEBUG 
printf("leaving computePhiInv\n");
#endif
}
@}

@d computePhiInv argument list
@{
unsigned int * maxNumberHElements, 
unsigned int hrows, unsigned int hcols,
 double * hmat, unsigned int * hmatj, unsigned int * hmati,
unsigned int brows, unsigned int bcols,
 double * bmat, unsigned int * bmatj, unsigned int * bmati,
 double * phiInvmat, unsigned int * phiInvmatj, unsigned int * phiInvmati
@}
@d computePhiInv declarations
@{
double * hzero;unsigned int *    hzeroj;unsigned int * hzeroi;
double * hplus;unsigned int *    hplusj;unsigned int * hplusi;
double * bright;unsigned int *    brightj;unsigned int * brighti;
double * hb;unsigned int * hbj;unsigned int * hbi;
unsigned int * iw;
unsigned int ierr[1];
unsigned int nr[1];
unsigned int nc[1];
unsigned int firstColumn[1];
unsigned int lastColumn[1];
unsigned int aOne[1]={1};
@}
@d computePhiInv allocations
@{
hzero = (double *)calloc(*maxNumberHElements,sizeof(double));
hzeroj = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
hzeroi = (unsigned int *)calloc(hrows+1,sizeof(unsigned int));
hplus = (double *)calloc(*maxNumberHElements,sizeof(double));
hplusj = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
hplusi = (unsigned int *)calloc(hrows+1,sizeof(unsigned int));
bright = (double *)calloc(*maxNumberHElements,sizeof(double));
brightj = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
brighti = (unsigned int *)calloc(brows+1,sizeof(unsigned int));
hb = (double *)calloc(*maxNumberHElements,sizeof(double));
hbj = (unsigned int *)calloc(*maxNumberHElements,sizeof(unsigned int));
hbi = (unsigned int *)calloc(hrows+1,sizeof(unsigned int));
iw = (unsigned int *)calloc(hrows,sizeof(unsigned int));
@}
@d extract hzero
@{
*firstColumn=bcols+1;
*lastColumn=bcols+hrows;
#ifdef DEBUG 
printf("in computePhiInv prunsigned inting hmat\n");
cPrunsigned intSparse(hrows,hmat,hmatj,hmati);
#endif
extractSubmatrix(&hrows,aOne,aOne,&hrows,firstColumn,lastColumn,
hmat,hmatj,hmati,nr,nc,hzero,hzeroj,hzeroi);
#ifdef DEBUG 
printf("in computePhiInv prunsigned int hzero\n");
cPrunsigned intSparse(hrows,hzero,hzeroj,hzeroi);
#endif
@}
@d extract hplus
@{
*firstColumn=bcols+hrows+1;
*lastColumn=hcols;
extractSubmatrix(&hrows,aOne,aOne,&hrows,firstColumn,lastColumn,
hmat,hmatj,hmati,nr,nc,hplus,hplusj,hplusi);
#ifdef DEBUG 
printf("in computePhiInv prunsigned inting hplus\n");
cPrunsigned intSparse(hrows,hplus,hplusj,hplusi);
#endif
@}
@d extract bright
@{
*firstColumn=bcols-hrows+1;
*lastColumn=bcols;
#ifdef DEBUG 
printf("in computePhiInv prunsigned inting bmat\n");
cPrunsigned intSparse(brows,bmat,bmatj,bmati);
#endif
extractSubmatrix(&brows,aOne,aOne,&brows,firstColumn,lastColumn,
bmat,bmatj,bmati,nr,nc,bright,brightj,brighti);
#ifdef DEBUG 
printf("in computePhiInv prunsigned inting bright\n");
cPrunsigned intSparse(hrows,bright,brightj,brighti);
#endif

@}

@d compute hplus times bright
@{
sparseMult(&hrows,&hrows,aOne,hplus,hplusj,hplusi,bright,brightj,brighti,
                       hb,hbj,hbi,maxNumberHElements,iw,ierr);
pathNewtAssert(ierr == 0);
#ifdef DEBUG 
printf("in computePhiInv prunsigned inting hb\n");
cPrunsigned intSparse(hrows,hb,hbj,hbi);
#endif
@}
@d add hzero
@{
sparseAdd(&hrows,&hrows,maxNumberHElements,iw,aOne,hzero,hzeroj,hzeroi,hb,hbj,hbi,
          phiInvmat,phiInvmatj,phiInvmati,ierr);
pathNewtAssert(ierr == 0);
#ifdef DEBUG 
printf("in computePhiInv prunsigned inting phiInvmat\n");
cPrunsigned intSparse(hrows,phiInvmat,phiInvmatj,phiInvmati);
#endif

@}


@d computePhiInv deallocations
@{
free(hzero);
free(hzeroj);
free(hzeroi);
free(hplus);
free(hplusj);
free(hplusi);
free(bright);
free(brightj);
free(brighti);
free(hb);
free(hbj);
free(hbi);
free(iw);
@}


@d another pn check for convergence
@{
#ifdef DEBUG 
printf("pathNewt:another testing for convergence\n");
#endif


@<pathNewt check for convergence@>


@}


@d pathNewt check for convergence
@{

		test=testx=testf=0.0;testfloc=99;
		for (i=0;i<n;i++)
			if (fabs(fvec[i]) > test) {testf=test=fabs(fvec[i]);testfloc=i;}
printf("pathNewt:checking test<TOLF return,max abs(fvec[i])=%e for i=%d\n",test,testfloc);/*
badArgs(i,fvec[i],chkfdrvj,chkfdrvi,0,*numberOfEquations); */
#ifdef DEBUG 
printf("pathNewt:checking test<TOLF return,max abs(fvec[i])=%e for i = %d\n",test,testfloc);
#endif
		if (test < (tolfInput)) {
resetFailedQ;
			*check=0;
			done=1;
printf("pathNewt:test<TOLF return,its=%d,max abs(fvec[i])=%e for i = %d\n",its,test,testfloc);
#ifdef DEBUG 
printf("pathNewt:test<TOLF return,its=%d,max abs(fvec[i])=%e for i = %d\n",its,test,testfloc);
#endif
		} else { 
		testx=test=0.0;
		for (i=0;i<n;i++) {
			temp=(fabs(x[i]-xold[i]))/(FMIN(fabs(x[i]),fabs(xold[i]))+GAMMA);
			if (temp > test) (testx=test=temp);
		}
printf("pathNewt:checking test<TOLX return, max(abs(x[i]-xold[i])/(min(x[i],xold[i])+gamma))=%e\n",test);
		if ((its>0)&&(test < tolxInput)) {
resetFailedQ;
			*check=0;
done=1;
printf("pathNewt:test<TOLX return,its=%d, max(abs(x[i]-xold[i])/(min(abs(x[i])),abs(xold[i]))+gamma)=%e\n",its,test);
printf("pathNewt:with TOLF max abs(fvec[i])=%e at i=%d\n",testf,testfloc);
#ifdef DEBUG 
printf("pathNewt:test<TOLX return,its=%d, max(abs(x[i]-xold[i])/(min(abs(x[i])),abs(xold[i]))+gamma)=%e\n",its,test);
#endif
}}

@}
@d pathNewt check for convergence 
@{
	for (sum=0.0,i=0;i<n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
		fold=f;
		for (i=0;i<n;i++) xold[i]=x[i];

@}


@d pathNewt update path
@{
		for (i=0;i<n;i++) xoldls[i]=x[i];
addOneToNewtonSteps;
if(useStackQ){
printf("using stack\n");
nxtGuess(numberOfEquations,lags,leads,pathLength,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,qMat,qMatj,qMati,fixedPoint,intercept,x,/*shockVec,*/xdel/*,
intControlParameters,doubleControlParameters*/);
bump(*maxNumberElements);
} else {
/* code to replace stack with sparse inversion*/
printf("not using stack\n");
constructFdrv(*numberOfEquations,*lags,*leads,*pathLength,
x,params,vecfunc,fdjac,qMat,qMatj,qMati,fixedPoint,intercept,linearizationPoint,/*exogRows,exogCols,exogenizeQ,*/
shockVec,
compfvec,chkfdrv,chkfdrvj,chkfdrvi,ihomotopy,
intControlParameters, doubleControlParameters,
intOutputInfo/*, doubleOutputInfo*/);

multMax=  (*pathLength)* *maxNumberElements+forQ;/* try 10 for multmax*/
newNxtGuess(&n,&multMax,
compfvec,chkfdrv,chkfdrvj,chkfdrvi,xdel,
ma50bdJob,
/*ma50bdIq,*/
ma50bdFact,
ma50bdIrnf,
ma50bdIptrl,
ma50bdIptru,
intControlParameters,doubleControlParameters
);
bump((((multMax-forQ)/(*pathLength))));
}
#ifdef DEBUG 
/*
printf("xdel=\n");
cPrunsigned intMatrixNonZero(1,n,xdel,1e-8);
printf("x=\n");
cPrunsigned intMatrixNonZero(1,n,x,1e-8);
*/


/* need to construct if using stack*/
if(useStackQ){
constructFdrv(*numberOfEquations,*lags,*leads,*pathLength,
x,params,vecfunc,fdjac,qMat,qMatj,qMati,fixedPoint,unsigned intercept,linearizationPoint,exogRows,exogCols,exogenizeQ,
shockVec,
compfvec,chkfdrv,chkfdrvj,chkfdrvi,ihomotopy,
intControlParameters, doubleControlParameters,
intOutputInfo, doubleOutputInfo);
}
chkDrv(*numberOfEquations*(*lags+*leads+*pathLength),
chkfdrv,chkfdrvj,chkfdrvi,compfvec,xdel);

#endif
		for (i=0;i<n;i++) {x[i]=xold[i]-xdel[i];lastDel[i]=xdel[i]; }
{/*use shock for first only*/
addOneToFEvals;
((*vecfunc)(x,params,shockVec,compfvec,compfvecj,compfveci,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/));}

		for(i=1;i<*pathLength;i++){
addOneToFEvals;
((*vecfunc)(x+(i* *numberOfEquations),params,zeroShockVec,compfvec+i* *numberOfEquations,compfvecj,compfveci,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/));}


if(expandFactorInput>0){
for(icnt=0;!lclValidVectorVerbose(*pathLength* *numberOfEquations,compfvec, *numberOfEquations,x+(i* *numberOfEquations),chkfdrvj,chkfdrvi,
 *numberOfEquations * *lags )&&icnt<5;icnt++){
addOneToExpandSteps;
printf("expanding...");
for (i=0;i<n;i++) {xdel[i]= xdel[i] *(expandFactorInput);}
		for (i=0;i<n;i++) x[i]=xold[i]-xdel[i];

{/* use shock for first period only*/
addOneToFEvals;
((*vecfunc)(x,params,shockVec,compfvec,compfvecj,compfveci,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/));}

		for(i=1;i<*pathLength;i++){
addOneToFEvals;
((*vecfunc)(x+(i* *numberOfEquations),params,zeroShockVec,compfvec+i* *numberOfEquations,compfvecj,compfveci,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/));}


}}
if(shrinkFactorInput>0){
for(icnt=0;!lclValidVectorVerbose(*pathLength* *numberOfEquations,compfvec,*numberOfEquations,x+(i* *numberOfEquations),chkfdrvj,chkfdrvi,
 *numberOfEquations * *lags)&&icnt<20;icnt++){
addOneToShrinkSteps;
printf("shrinking...");
for (i=0;i<n;i++) {xdel[i]= xdel[i] *(shrinkFactorInput);}
		for (i=0;i<n;i++) x[i]=xold[i]-xdel[i];


{/*use shock for first period only*/
addOneToFEvals;
((*vecfunc)(x,params,shockVec,compfvec,compfvecj,compfveci,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/));}
		for(i=1;i<*pathLength;i++){
addOneToFEvals;
((*vecfunc)(x+(i* *numberOfEquations),params,shockVec,compfvec+i* *numberOfEquations,compfvecj,compfveci,homotopyAlpha+ihomotopy,linearizationPoint/*,exogRows,exogCols,exogenizeQ*/));}


}}
if(icnt>20)printf("shrink failed to get rid of NaN and other non  numbers!!!!!!!!!\n");

		for (i=0;i<n;i++) x[i]=xold[i];

if(*pathLength)  {
if(useLnsrchQ){
@<compute gradient for lnsrch@>
		lnsrch(n,*numberOfEquations,*pathLength,
        xoldls,&fold,g,xdel,params,shockVec,&f,stpmax,check,vecfunc,x,ihomotopy,linearizationPoint,/*exogRows,exogCols,exogenizeQ,intControlParameters,*/doubleControlParameters/*,
intOutputInfo, doubleOutputInfo*/);
} else {for(i=*numberOfEquations* *lags;i<n;i++)x[i]=x[i]-xdel[i];
}
        } else {
for(i=*numberOfEquations* *lags;i<n;i++)x[i]=x[i]-xdel[i];
}


@}
@d compute gradient for lnsrch
@{
if(useStackQ && useLnsrchQ){
/* need to do this if using stack and lnsrch
constructFdrv(*numberOfEquations,*lags,*leads,*pathLength,
x,params,vecfunc,fdjac,qMat,qMatj,qMati,fixedPoint,unsigned intercept,linearizationPoint,exogRows,exogCols,exogenizeQ,
shockVec,
compfvec,chkfdrv,chkfdrvj,chkfdrvi,ihomotopy,
intControlParameters, doubleControlParameters,
intOutputInfo, doubleOutputInfo);
*/
}
        sparseMatTimesVec(&n,chkfdrv,chkfdrvj,chkfdrvi,compfvec,
        g);
for(i=0;i<*pathLength;i++){
sparseMatTimesVec(rowDim,smats[i],
smatsj[i],smatsi[i],fvec+(i * *numberOfEquations),
g+(*numberOfEquations * (*lags+i)));}

@}


@d pathNewt initializations
@{
/*
resetFailedQ;
resetNewtonSteps;
resetRealizedTolf;
resetRealizedTolx;
resetFEvals;
resetFDrvEvals;
resetShrinkSteps;
resetExpandSteps;
resetLnsrchSteps;
*/
  *check=0;
    n=*numberOfEquations*(*lags+*leads+*pathLength);
	rowDim=(unsigned int *)calloc(1,sizeof(unsigned int));
	*rowDim=*numberOfEquations**leads;
    aOne=(unsigned int *)calloc(1,sizeof(unsigned int));
    *aOne=1;
    aZero=(unsigned int *)calloc(1,sizeof(unsigned int));
    *aZero=0;
    qColumns=(unsigned int *)calloc(1,sizeof(unsigned int));
    *qColumns=*numberOfEquations*(*leads+*lags+1);
    deviations=(double *)calloc(*qColumns,sizeof(double));
    zeroShockVec=(double *)calloc(*numberOfEquations,sizeof(double));
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
forQ=qMati[*numberOfEquations * *leads]-1;
printf("setting forQ to %d\n",forQ);
compfvec = (double * ) calloc(n,sizeof(double));
compfvecj = (unsigned int * ) calloc(n,sizeof(unsigned int));
compfveci = (unsigned int * ) calloc(n+1,sizeof(unsigned int));
chkfdrv = (double * )calloc(*pathLength * *maxNumberElements+
/**pathLength* *pathLength * *maxNumberElements+*/forQ,sizeof(double));
chkfdrvj = (unsigned int * )calloc(*pathLength * *maxNumberElements+
/**pathLength* *pathLength * *maxNumberElements+*/forQ,sizeof(unsigned int));
chkfdrvi=(unsigned int * ) calloc(n+1,sizeof(unsigned int));


@}
@d pathNewt declarations
@{    unsigned int  n; unsigned int done;/*unsigned unsigned int ihomotopy=1;*/
unsigned int ihomotopy=0;double * zeroShockVec;
//#include "stochSims.h"
unsigned int forQ/*;double  dignore[1]={1.0}*/;
unsigned int maxElementsEncountered=0;
unsigned int icnt;
unsigned int  multMax;
    unsigned int tNow;unsigned int * rowDim;unsigned int * qColumns;
    double * deviations;double * fullfvec;
/*	double fmin(double x[]);*/
/* declared below
	void lnsrch(unsigned int  n, unsigned int np,unsigned int reps,double xold[], double  *fold, double g[], double p[],double * params,double * shockVec,
		 double * f, double stpmax, unsigned int *check, 
void (*func)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),double * x,
            unsigned int ihomotopy,double * linPt,unsigned int * exogRows, unsigned int * exogCols, unsigned int * exogenizeQ,
unsigned int * intControlParameters,double * doubleControlParameters,
unsigned int * intOutputInfo, double * doubleOutputInfo);
*/

	unsigned int i,its/*,j*/,*indx,*aOne,*aZero/*,*ndns*/,*ierr,testfloc;
	double /*d,den,*/f,fold,stpmax,sum,temp,test,testx,testf,*g,*p,*xold,*xoldls,*xdel,normSum;
double * compfvec;
unsigned int * compfvecj;
unsigned int * compfveci;
double * chkfdrv;
unsigned int * chkfdrvj;
unsigned int * chkfdrvi;

@}

p
\subsection{lnsrch}
\label{sec:lnsrch}

@d lnsrch argument list
@{
unsigned int  n,unsigned int np,unsigned int reps,
double xold[], double * fold, double g[], double p[], 
		 double * params,double * shockVec,double * f,double stpmax, unsigned int *check,
void (*func)(double*, double*, double*,double*,unsigned int*,unsigned int*,double *,double*),
double * x,
            unsigned int ihomotopy,double * linPt,/*unsigned int * exogRows,*/
/*unsigned int * exogCols,unsigned int * exogenizeQ,*/
/*unsigned int * intControlParameters,*/double * doubleControlParameters/*,
unsigned int * intOutputInfo, double * doubleOutputInfo*/
@}
@o stackC.h -d
@{
#include <unistd.h>
//pid_t getpid(void);


void lnsrch(@<lnsrch argument list@>);
@}
@o myNewt.c -d
@{
#include <math.h>
#define NRANSI


extern void dgemm_(char * noTransp1,char * noTransp2,
unsigned int * neq,unsigned int * aOne,unsigned int *BMatrixColumns1,
double *aFloatOne,
double *__asymptoticBmatrices,unsigned int *BMatrixRows1,
double *deviationsFromPeriodicPath,unsigned int *BMatrixColumns2,
double *aZero,
double *pathNow,
unsigned int *BMatrixRows2);



#define TRUE 1
unsigned int lclValidVector(unsigned int numRows,double * vec)
{
unsigned int i,allFiniteNumbers;
unsigned int printed=0;
      allFiniteNumbers=TRUE;
      for(i=0;i<numRows;i++){
        allFiniteNumbers=(isfinite(vec[i])&&allFiniteNumbers);
        if(!isfinite(vec[i])) {
if(printed<2) {
printf("<lclValidVectorNew:problem with i=%d,vec[i]=%e,with numRows=%d>",i,vec[i],numRows);printed++;} else 
{printf(".");}
}
}
printf("\n");
      return(allFiniteNumbers);
}
@}
@o stackC.h -d
@{
unsigned int lclValidVectorVerbose(@<lclValidVectorVerbose argument list@>);
@}
@d lclValidVectorVerbose argument list
@{
unsigned int numRows,double * vec,unsigned int  forMod,double * x,unsigned int * ja,unsigned int * ia,unsigned int skipRows
@}
@o myNewt.c -d
@{unsigned int lclValidVectorVerbose(@<lclValidVectorVerbose argument list@>)
{
unsigned int i,allFiniteNumbers;
unsigned int printed=0;
      allFiniteNumbers=TRUE;
      for(i=0;i<numRows;i++){
        allFiniteNumbers=(isfinite(vec[i])&&allFiniteNumbers);
        if(!isfinite(vec[i])) {
if(printed<1) {
printf("<lclValidVectorVerbose:problem with eqn i=%d in block=%d,vec[i]=%e,with numRows=%d>",i%forMod+1,i/forMod,vec[i],numRows);printed++;
badArgs(i,x,ja,ia,skipRows,forMod);} else 
{printf(".");}
}
}
printf("\n");
      return(allFiniteNumbers);
}
@}
@d badArgs signature
@{
void badArgs(unsigned int ii,double * x,unsigned int * ja,unsigned int * ia,unsigned int skipRows,unsigned int forMod)
@}
@o stackC.h -d
@{
@<badArgs signature@>;
@}
@o myNewt.c -d
@{


@<badArgs signature@>{
unsigned int kk;unsigned int fir;unsigned int lst;//unsigned int elems;
fir=ia[ii+skipRows]-1;
lst=ia[ii+skipRows+1]-1;
//elems=lst-fir;
for(kk=fir;kk<lst;kk++){
printf("x[%d,%d]=%e ",((ja[kk]-1)/forMod)-1,((ja[kk]-1)%forMod)+1,x[ja[kk]-1]);
}
}

extern double  dnrm2_();
void lnsrch(@<lnsrch argument list@>)
{
@< lnsrch declarations@>
@< lnsrch preloop@>
@< lnsrch loop@>
@< lnsrch postloop@>
free(aOne);free(aZero);free(aTwo);free(fvec);free(fvecj);free(fveci);
}

#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software '>9m_L31.. */

@}

@d lnsrch declarations
@{
	unsigned int i;double * zeroShockVec;
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
	zeroShockVec=(double *)calloc(n,sizeof(double));
	xorig=(double *)calloc(n,sizeof(double));
	fvecj=(unsigned int *)calloc(n,sizeof(double));

@}
@d lnsrch preloop
@{


/*
    fold = dnrm2_(&np,foldp,aOne);
*/

*fold=0;
{/* use shock for first only*/
addOneToFEvals;
(*func)(xold,params,shockVec,fvec,fvecj,fveci,homotopyAlpha+ihomotopy,linPt/*,exogRows,exogCols,exogenizeQ*/);
        useCNRMS(&np,aTwo,fvec,fvecj,fveci,xorig);
        *fold += xorig[0];
        }
		for(i=1;i<reps;i++){
addOneToFEvals;
(*func)(xold+(i*np),params,zeroShockVec,fvec,fvecj,fveci,homotopyAlpha+ihomotopy,linPt/*,exogRows,exogCols,exogenizeQ*/);
        useCNRMS(&np,aTwo,fvec,fvecj,fveci,xorig);
        *fold += xorig[0];
        }

      *fold *= 0.5;

/*
csrToDns(&n,aOne,fvec,fvecj,fveci,xorig,&n,ierr);
pathNewtAssert(ierr == 0);

    dgemm_(transp,noTransp,
           aOne,aOne,&n,aDoubleOne,
           xorig,&n,
           p,&n,
           aDoubleZero,dir,aOne);
if(*dir<0)
*/

for (i=0;i<n;i++) p[i]= (-p[i]);

resetFailedQ;
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
		if (temp > test) (test=temp);
	}
	if(test>tolxInput){alamin=1e-5*tolxInput/test;} else {alamin=(alaminInput);}
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
			free(zeroShockVec);
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


*f=0;
{/*use shock first period only*/
addOneToFEvals;
(*func)(x,params,shockVec,fvec,fvecj,fveci,homotopyAlpha+ihomotopy,linPt/*,exogRows,exogCols,exogenizeQ*/);
        useCNRMS(&np,aTwo,fvec,fvecj,fveci,xorig);
        *f += xorig[0];
        }
		for(i=1;i<reps;i++){
addOneToFEvals;
(*func)(x+(i*np),params,zeroShockVec,fvec,fvecj,fveci,homotopyAlpha+ihomotopy,linPt/*,exogRows,exogCols,exogenizeQ*/);
        useCNRMS(&np,aTwo,fvec,fvecj,fveci,xorig);
        *f += xorig[0];
        }

      *f *= 0.5;

#define DEBUG 1
#ifdef DEBUG 
printf("lnsrch:f %f\n",*f);
printf("lnsrch:fold %f\n",*fold);
printf("lnsrch:alam %f\n",alam);
printf("lnsrch:alamin %f\n",alamin);
#endif

		if (alam < NEGLIGIBLEDOUBLE ) {
		  			for (i=0;i<n;i++) x[i]=xold[i]+p[i];
			 /**check=1;*/
@<lnsrch free storage@>
printf("lnsrch:alam<NEGLIGIBLEDOUBLE, just do full newton step\n");
#ifdef DEBUG 
printf("lnsrch:alam<NEGLIGIBLEDOUBLE, just do full newton step\n");
#endif
			return;
		} else 	if (alam < alamin) {
		  			for (i=0;i<n;i++) x[i]=xold[i]+p[i];
printf("lnsrch:alam<alamin, just do full newton step\n");
			/* *check=1;*/
@<lnsrch free storage@>
#ifdef DEBUG 
printf("lnsrch:alam<alamin, just do full newton step\n");
#endif
			return;
		} else if (*f < *fold+(alfInput)*alam*slope) {
@<lnsrch free storage@>
#ifdef DEBUG 
printf("lnsrch:*f < *fold+(alfInput)*alam*slope return\n,slope=%e,alam=%e\n,*fold=%e,f=%e\n",slope,alam,*fold,*f);
#endif
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


@o validArrayElements.c
@{
@<define assert bump@>
unsigned int validArrayElements(unsigned int elements,double * mata)
{
unsigned int result,allPositive,allFiniteNumbers;
result=(elements >0);
      allFiniteNumbers=TRUE;
      for(i=0;i<elements;i++){
        allFiniteNumbers=(isfinite(mata[i])&&allFiniteNumbers);}
      result=(result &&  allFiniteNumbers);
      return(result);
}
@}




\section{Development Notes}
\label{sec:devnotes}


\subsection{Known Bugs}
\label{sec:knownbugs}


\subsection{Things to Do}
\label{sec:thingstodo}

\begin{itemize}
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
theSparseFunction, theSparseFunctionDerivative, exogRows,exogCols,exogenizeQ,paramVector, 
fp, 
fmats, fmatsj, fmatsi, 
smats, smatsj, smatsi, 
maxNumberElements, 
linearConstrainMatrix, 
ierr
);
)
\item[{\bf  pathNewt}] pathNewt(NEQS, NLAGS, NLEADS, pathLength, 
theSparseFunction, theSparseFunctionDerivative, paramVector, shockVector,
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
FIU: file descriptor 26: <reserved for Purify unsigned internal use>
FIU: file descriptor 27: <reserved for Purify unsigned internal use>

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
        1889 libunsigned internal_stubs.so.1 (shared code)
         208 libunsigned internal_stubs.so.1 (private data)


\end{verbatim}

\subsection{Diagnostic Prunsigned inting}
\label{sec:diag}

Use emacs to enter:
\begin{verbatim}
gdb testStack -x gdbCommands
\end{verbatim}


\subsection{SPARSEKIT2}
\label{sec:sparsekit2}

\subsubsection{dump}
@o stackC.h -d
@{
void dump_();
@}


\begin{verbatim}
c----------------------------------------------------------------------- 
      subroutine dump (i1,i2,values,a,ja,ia,iout)
      unsigned integer i1, i2, ia(*), ja(*), iout
      real*8 a(*) 
      logical values 
c-----------------------------------------------------------------------
c outputs rows i1 through i2 of a sparse matrix stored in CSR format 
c (or columns i1 through i2 of a matrix stored in CSC format) in a file, 
c one (column) row at a time in a nice readable format. 
c This is a simple routine which is useful for debugging. 
c-----------------------------------------------------------------------
c on entry:
c---------
c i1    = first row (column) to prunsigned int out
c i2    = last row (column) to prunsigned int out 
c values= logical. indicates whether or not to prunsigned int real values.
c         if value = .false. only the pattern will be output.
c a,
c ja, 
c ia    =  matrix in CSR format (or CSC format) 
c iout  = logical unit number for output.
c---------- 
c the output file iout will have written in it the rows or columns 
c of the matrix in one of two possible formats (depending on the max 
c number of elements per row. The values are output with only 
c two digits of accuracy (D9.2). )

\end{verbatim}

\subsubsection{amub}
@o stackC.h -d
@{
//void sparseMult();
@}

\begin{verbatim}

/*
c----------------------------------------------------------------------c
       subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                  c,jc,ic,nzmax,iw,ierr) 
      real*8 a(*), b(*), c(*) 
      unsigned integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix by matrix product C = A B 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = unsigned integer. The row dimension of A = row dimension of C
c ncol  = unsigned integer. The column dimension of B = column dimension of C
c job   = unsigned integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = unsigned integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic    = resulting matrix C in compressed sparse row sparse format.
c           
c ierr  = unsigned integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = unsigned integer work array of length equal to the number of
c         columns in A.
c Note: 
c-------
c   The row dimension of B is not needed. However there is no checking 
c   on the condition that ncol(A) = nrow(B). 
c
c----------------------------------------------------------------------- 

*/
\end{verbatim}

\subsubsection{amux}

@o stackC.h -d
@{
//void sparseMatTimesVec();











@}


\begin{verbatim}

/*
c----------------------------------------------------------------------c
      subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      unsigned integer n, ja(*), ia(*)
c-----------------------------------------------------------------------
c         A times a vector
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
*/
\end{verbatim}




\subsubsection{aplb}
@o stackC.h -d
@{
//void sparseAdd();
@}


\begin{verbatim}

/*
c-----------------------------------------------------------------------
      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      unsigned integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix sum  C = A+B. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = unsigned integer. The row dimension of A and B
c ncol  = unsigned integer. The column dimension of A and B.
c job   = unsigned integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = unsigned integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic    = resulting matrix C in compressed sparse row sparse format.
c       
c ierr  = unsigned integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = unsigned integer work array of length equal to the number of
c         columns in A.
c
c-----------------------------------------------------------------------
*/
\end{verbatim}


\subsubsection{csrcsc}

@o stackC.h -d
@{
//void csrcsc2_();
@}

\begin{verbatim}
c-----------------------------------------------------------------------
      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      unsigned integer ia(n+1),iao(n2+1),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c Rectangular version.  n is number of rows of CSR matrix,
c                       n2 (input) is number of columns of CSC matrix.
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n = number of rows of CSR matrix.
c n2    = number of columns of CSC matrix.
c job   = unsigned integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this unsigned into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
c     for any other normal usage, enter ipos=1.
c a = real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja    = unsigned integer array of length nnz containing the column positions
c     of the corresponding elements in a.
c ia    = unsigned integer of size n+1. ia(k) contains the position in a, ja of
c     the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao    = real array of size nzz containing the "a" part of the transpose
c jao   = unsigned integer array of size nnz containing the column indices.
c iao   = unsigned integer array of size n+1 containing the "ia" index array of
c     the transpose. 
c
c----------------------------------------------------------------------- 
\end{verbatim}
@o stackC.h -d
@{
//void csrToCsc();
@}


\subsubsection{csrcsc}
\begin{verbatim}
c----------------------------------------------------------------------- 
      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      unsigned integer ia(n+1),iao(n+1),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n	= dimension of A.
c job	= unsigned integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this unsigned into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= unsigned integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= unsigned integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= unsigned integer array of size nnz containing the column indices.
c iao	= unsigned integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 

\end{verbatim}

\subsubsection{csrdns}
@o stackC.h -d
@{
//void csrToDns();
@}


\begin{verbatim}

/*c----------------------------------------------------------------------c
      subroutine csrdns(nrow,ncol,a,ja,ia,dns,ndns,ierr) 
      real*8 dns(ndns,*),a(*)
      unsigned integer ja(*),ia(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row    to    Dense 
c-----------------------------------------------------------------------
c
c converts a row-stored sparse matrix unsigned into a densely stored one
c
c On entry:
c---------- 
c
c nrow  = row-dimension of a
c ncol  = column dimension of a
c a, 
c ja, 
c ia    = input matrix in compressed sparse row format. 
c         (a=value array, ja=column array, ia=pounsigned inter array)
c dns   = array where to store dense matrix
c ndns  = first dimension of array dns 
c
c on return: 
c----------- 
c dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
c 
c ierr  = unsigned integer error indicator. 
c         ierr .eq. 0  means normal return
c         ierr .eq. i  means that the code has stopped when processing
c         row number i, because it found a column number .gt. ncol.
c 
c----------------------------------------------------------------------- 
*/
\end{verbatim}


\subsubsection{dnscsr}
@o stackC.h -d
@{
//void dnsToCsr();
@}


\begin{verbatim}

/*
c----------------------------------------------------------------------- 
      subroutine dnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)
      real*8 dns(ndns,*),a(*)
      unsigned integer ia(*),ja(*)
c-----------------------------------------------------------------------
c Dense     to    Compressed Row Sparse 
c----------------------------------------------------------------------- 
c
c converts a densely stored matrix unsigned into a row orientied
c compactly sparse matrix. ( reverse of csrdns )
c Note: this routine does not check whether an element 
c is small. It considers that a(i,j) is zero if it is exactly
c equal to zero: see test below.
c-----------------------------------------------------------------------
c on entry:
c---------
c
c nrow  = row-dimension of a
c ncol  = column dimension of a
c nzmax = maximum number of nonzero elements allowed. This
c         should be set to be the lengths of the arrays a and ja.
c dns   = input nrow x ncol (dense) matrix.
c ndns  = first dimension of dns. 
c
c on return:
c---------- 
c 
c a, ja, ia = value, column, pounsigned inter  arrays for output matrix 
c
c ierr  = unsigned integer error indicator: 
c         ierr .eq. 0 means normal retur
c         ierr .eq. i means that the the code stopped while
c         processing row number i, because there was no space left in
c         a, and ja (as defined by parameter nzmax).
c----------------------------------------------------------------------- 

*/
\end{verbatim}


\subsubsection{submat}
@o stackC.h -d
@{
//void extractSubmatrix();
@}

\begin{verbatim}

/*
      subroutine submat (n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)
      unsigned integer n,job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
      real*8 a(*),ao(*) 
c-----------------------------------------------------------------------
c extracts the submatrix A(i1:i2,j1:j2) and puts the result in 
c matrix ao,iao,jao
c---- In place: ao,jao,iao may be the same as a,ja,ia.
c-------------- 
c on input
c---------
c n = row dimension of the matrix 
c i1,i2 = two unsigned integers with i2 .ge. i1 indicating the range of rows to be
c          extracted. 
c j1,j2 = two unsigned integers with j2 .ge. j1 indicating the range of columns 
c         to be extracted.
c         * There is no checking whether the input values for i1, i2, j1,
c           j2 are between 1 and n. 
c a,
c ja,
c ia    = matrix in compressed sparse row format. 
c
c job   = job indicator: if job .ne. 1 then the real values in a are NOT
c         extracted, only the column indices (i.e. data structure) are.
c         otherwise values as well as column indices are extracted...
c         
c on output
c-------------- 
c nr    = number of rows of submatrix 
c nc    = number of columns of submatrix 
c     * if either of nr or nc is nonpositive the code will quit.
c
c ao,
c jao,iao = extracted matrix in general sparse format with jao containing
c   the column indices,and iao being the pounsigned inter to the beginning 
c   of the row,in arrays a,ja.
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
*/
\end{verbatim}



\subsection{HARWELL Routines}
\label{sec:harwell}


\subsubsection{MA50AD}

@o stackC.h -d
@{
//void useMA50AD();
@}

\begin{verbatim}
      SUBROUTINE MA50AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP,JFIRST,
     +                  LENR,LASTR,NEXTR,IW,IFIRST,LENC,LASTC,NEXTC,
     +                  INFO,RINFO)

C MA50A/AD chooses a pivot sequence using a Markowitz criterion with
C     threshold pivoting.

C If  the user requires a more convenient data unsigned interface then the MA48
C     package should be used. The MA48 subroutines call the MA50
C     subroutines after checking the user's input data and optionally
C     permute the matrix to block triangular form.

      UNSIGNED INTEGER M,N,NE,LA
      DOUBLE PRECISION A(LA)
      UNSIGNED INTEGER IRN(LA),JCN(LA),IQ(N)
      DOUBLE PRECISION CNTL(4)
      UNSIGNED INTEGER ICNTL(7),IP(M),NP,JFIRST(M),LENR(M),LASTR(M),NEXTR(M),
     +        IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(7)
      DOUBLE PRECISION RINFO

C M is an unsigned integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an unsigned integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an unsigned integer variable that must be set to the number of entries
C      in the input matrix. It is not altered by the subroutine.
C LA is an unsigned integer variable that must be set to the size of A, IRN, and
C      JCN. It is not altered by the subroutine.
C A is an array that holds the input matrix on entry and is used as
C      workspace.
C IRN  is an unsigned integer array.  Entries 1 to NE must be set to the
C      row indices of the corresponding entries in A.  IRN is used
C      as workspace and holds the row indices of the reduced matrix.
C JCN  is an unsigned integer array that need not be set by the user. It is
C      used to hold the column indices of entries in the reduced
C      matrix.
C IQ is an unsigned integer array of length N. On entry, it holds pounsigned inters
C      to column starts. During execution, IQ(j) holds the position of
C      the start of column j of the reduced matrix or -IQ(j) holds the
C      column index in the permuted matrix of column j. On exit, IQ(j)
C      holds the index of the column that is in position j of the
C      permuted matrix.


C CNTL must be set by the user as follows and is not altered.
C     CNTL(1)  Full matrix processing will be used if the density of
C       the reduced matrix is MIN(CNTL(1),1.0) or more.
C     CNTL(2) determines the balance between pivoting for sparsity and
C       for stability, values near zero emphasizing sparsity and values
C       near one emphasizing stability. Each pivot must have absolute
C       value at least CNTL(2) times the greatest absolute value in the
C       same column of the reduced matrix.
C     CNTL(3) If this is set to a positive value, any entry of the
C       reduced matrix whose modulus is less than CNTL(3) will be
C       dropped.
C     CNTL(4)  Any entry of the reduced matrix whose modulus is less
C       than or equal to CNTL(4) will be regarded as zero from the
C        pounsigned int of view of rank.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 3, plus all parameters on entry and exit.
C     ICNTL(4) If set to a positive value, the pivot search is limited
C       to ICNTL(4) columns (Zlatev strategy). This may result in
C       different fill-in and execution time. If ICNTL(4) is positive,
C       the workspace arrays LASTR and NEXTR are not referenced.
C     ICNTL(5) The block size to be used for full-matrix processing.
C     ICNTL(6) The last ICNTL(6) columns of A must be the last
C       ICNTL(6) columns of the permuted matrix. A value outside the
C       range 1 to N-1 is treated as zero.
C     ICNTL(7) If given the value 1, pivots are limited to
C       the main diagonal, which may lead to a premature switch to full
C       processing if no suitable diagonal entries are available.
C       If given the value 2, IFIRST must be set so that IFIRST(i) is
C       the column in position i of the permuted matrix and IP must
C       be set so that IP(i) < IP(j) if row i is recommended to
C       precede row j in the pivot sequence.
C IP is an unsigned integer array of length M that need not be set on entry
C      unless ICNTL(7)=2 (see ICNTL(7) for details of this case).
C      During execution, IP(i) holds the position of the start of row i
C      of the reduced matrix or -IP(i) holds the row index in the
C      permuted matrix of row i. Before exit, IP(i) is made positive.
C NP is an unsigned integer variable. It need not be set on entry. On exit,
C     it will be set to the number of columns to be processed in
C     packed storage.
C JFIRST is an unsigned integer workarray of length M. JFIRST(i) is the
C      first column of the reduced matrix to have i entries or is
C      zero if no column has i entries.
C LENR is an unsigned integer workarray of length M that is used to hold the
C      numbers of entries in the rows of the reduced matrix.
C LASTR is an unsigned integer workarray of length M, used only if ICNTL(4) = 0.
C      For rows in the reduced matrix, LASTR(i) indicates the previous
C      row to i with the same number of entries. LASTR(i) is zero if
C      no such row exists.
C NEXTR is an unsigned integer workarray of length M, used only if ICNTL(4) = 0
C      or ICNTL(7)=2. If ICNTL(4)=0, for rows in the reduced matrix,
C      NEXTR(i) indicates the next row to i with the same number of
C      entries; and if row i is the last in the chain, NEXTR is
C      equal to zero. If ICNTL(7)=2, NEXTR is a copy of the value of
C      IP on entry.
C IW is an unsigned integer array of length M used as workspace and is used to
C     assist the detection of duplicate entries and the sparse SAXPY
C     operations. It is reset to zero each time round the main loop.
C IFIRST is an unsigned integer array of length N, used only if ICNTL(4) = 0
C      or ICNTL(7)=2. If ICNTL(4) = 0, it is a workarray; IFIRST(i)
C      pounsigned ints to the first row of the reduced matrix to have i entries
C      or is zero if no row has i entries. If ICNTL(7)=2, IFIRST
C      must be set on entry (see ICNTL(7) for details of this case).
C LENC is an unsigned integer workarray of length N that is used to hold
C      the numbers of entries in the columns of the reduced matrix.
C LASTC is an unsigned integer workarray of length N.  For columns in the reduced
C      matrix, LASTC(j) indicates the previous column to j with the same
C      number of entries.  If column j is the first in the chain,
C      LASTC(j) is equal to zero.
C NEXTC is an unsigned integer workarray of length N.  For columns in the reduced
C      matrix, NEXTC(j) indicates the next column to j with the same
C      number of entries.  If column j is the last column in the chain,
C      NEXTC(j) is zero.
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1):
C       0  Successful entry.
C      -1  M < 1 or N < 1.
C      -2  NE < 1.
C      -3  Insufficient space.
C      -4  Duplicated entries.
C      -5  Faulty column permutation in IFIRST when ICNTL(7)=2.
C      -6  ICNTL(4) not equal to 1 when ICNTL(7)=2.
C      +1  Rank deficient.
C      +2  Premature switch to full processing because of failure to
C          find a stable diagonal pivot (ICNTL(7)>=1 case only).
C      +3  Both of these warnings.
C    INFO(2) Number of compresses of the arrays.
C    INFO(3) Minimum LA recommended to analyse matrix.
C    INFO(4) Minimum LFACT required to factorize matrix.
C    INFO(5) Upper bound on the rank of the matrix.
C    INFO(6) Number of entries dropped from the data structure.
C    INFO(7) Number of rows processed in full storage.
C RINFO need not be set on entry. On exit, it holds the number of
C    floating-pounsigned int operations needed for the factorization.


\end{verbatim}


\subsubsection{CNRMS}
\begin{verbatim}
c----------------------------------------------------------------------- 
      subroutine cnrms   (nrow, nrm, a, ja, ia, diag) 
      real*8 a(*), diag(nrow) 
      unsigned integer ja(*), ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each column of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= unsigned integer. The row dimension of A
c
c nrm   = unsigned integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------

\end{verbatim}

\subsubsection{MA50BD}

@o stackC.h -d
@{
//void useMA50BD();
@}

\begin{verbatim}
      
SUBROUTINE MA50BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP,IQ,NP,
     +                  LFACT,FACT,IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)
C MA50B/BD factorizes the matrix in AA/IRNA/IPTRA as P L U Q where
C     P and Q are permutations, L is lower triangular, and U is unit
C     upper triangular. The prior information that it uses depends on
C     the value of the parameter JOB.
C
      UNSIGNED INTEGER M,N,NE,JOB
      DOUBLE PRECISION AA(NE)
      UNSIGNED INTEGER IRNA(NE),IPTRA(N)
      DOUBLE PRECISION CNTL(4)
      UNSIGNED INTEGER ICNTL(7),IP(M),IQ(N),NP,LFACT
      DOUBLE PRECISION FACT(LFACT)
      UNSIGNED INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION W(M)
      UNSIGNED INTEGER IW(M+2*N),INFO(7)
      DOUBLE PRECISION RINFO
C
C M is an unsigned integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an unsigned integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an unsigned integer variable that must be set to the number of entries
C      in the input matrix.  It is not altered by the subroutine.
C JOB is an unsigned integer variable that must be set to the value 1, 2, or 3.
C     If JOB is equal to 1 and any of the first NP recommended pivots
C      fails to satisfy the threshold pivot tolerance, the row is
C      unsigned interchanged with the earliest row in the recommended sequence
C      that does satisfy the tolerance. Normal row unsigned interchanges are
C      performed in the last N-NP columns.
C     If JOB is equal to 2, then M, N, NE, IRNA, IPTRA, IP, IQ,
C      LFACT, NP, IRNF, IPTRL, and IPTRU must be unchanged since a
C      JOB=1 entry for the same matrix pattern and no unsigned interchanges are
C      performed among the first NP pivots; if ICNTL(6) > 0, the first
C      N-ICNTL(6) columns of AA must also be unchanged.
C     If JOB is equal to 3, ICNTL(6) must be in the range 1 to N-1.
C      The effect is as for JOB=2 except that unsigned interchanges are
C      performed.
C     JOB is not altered by the subroutine.
C AA is an array that holds the entries of the matrix and
C      is not altered.
C IRNA is an unsigned integer array of length NE that must be set to hold the
C      row indices of the corresponding entries in AA. It is not
C      altered.
C IPTRA is an unsigned integer array that holds the positions of the starts of
C      the columns of AA. It is not altered by the subroutine.
C CNTL  must be set by the user as follows and is not altered.
C     CNTL(2) determines the balance between pivoting for sparsity and
C       for stability, values near zero emphasizing sparsity and values
C       near one emphasizing stability.
C     CNTL(3) If this is set to a positive value, any entry whose
C       modulus is less than CNTL(3) will be dropped from the factors.
C       The factorization will then require less storage but will be
C       inaccurate.
C     CNTL(4)  Any entry of the reduced matrix whose modulus is less
C       than or equal to CNTL(4) will be regarded as zero from the
C        pounsigned int of view of rank.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 2, plus all parameters on entry and exit.
C     ICNTL(5) The block size to be used for full-matrix processing.
C       If <=0, the BLAS1 version is used.
C       If =1, the BLAS2 version is used.
C     ICNTL(6) If N > ICNTL(6) > 0, only the columns of A that
C       correspond to the last ICNTL(6) columns of the permuted matrix
C       may change prior to an entry with JOB > 1.
C IP is an unsigned integer array. If JOB=1, it must be set so that IP(I) < IP(J)
C      if row I is recommended to precede row J in the pivot sequence.
C      If JOB>1, it need not be set. If JOB=1 or JOB=3, IP(I) is set
C      to -K when row I is chosen for pivot K and IP is eventually
C      reset to recommend the chosen pivot sequence to a subsequent
C      JOB=1 entry. If JOB=2, IP is not be referenced.
C IQ is an unsigned integer array that must be set so that either IQ(J) is the
C      column in position J in the pivot sequence, J=1,2,...,N,
C      or IQ(1)=0 and the columns are taken in natural order.
C      It is not altered by the subroutine.
C NP is an unsigned integer variable that holds the number of columns to be
C      processed in packed storage. It is not altered by the subroutine.
C LFACT is an unsigned integer variable set to the size of FACT and IRNF.
C      It is not altered by the subroutine.
C FACT is an array that need not be set on a JOB=1 entry and must be
C      unchanged since the previous entry if JOB>1. On return, FACT(1)
C      holds the value of CNTL(3) used, FACT(2) will holds the value
C      of CNTL(4) used, FACT(3:IPTRL(N)) holds the packed part of L/U
C      by columns, and the full part of L/U is held by columns
C      immediately afterwards. U has unit diagonal entries, which are
C      not stored. In each column of the packed part, the entries of
C      U precede the entries of L; also the diagonal entries of L
C      head each column of L and are reciprocated.
C IRNF is an unsigned integer array of length LFACT that need not be set on
C      a JOB=1 entry and must be unchanged since the previous entry
C      if JOB>1. On exit, IRNF(1) holds the number of dropped entries,
C      IRNF(2) holds the number of rows MF in full storage,
C      IRNF(3:IPTRL(N)) holds the row numbers of the packed part
C      of L/U, IRNF(IPTRL(N)+1:IPTRL(N)+MF) holds the row indices
C      of the full part of L/U, and IRNF(IPTRL(N)+MF+I), I=1,2,..,N-NP
C      holds the vector IPIV output by MA50GD.
C      If JOB=2, IRNF will be unaltered.
C IPTRL is an unsigned integer array that need not be set on a JOB=1 entry and
C     must be unchanged since the previous entry if JOB>1.
C     For J = 1,..., NP, IPTRL(J) holds the position in
C     FACT and IRNF of the end of column J of L.
C     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
C IPTRU is an unsigned integer array that need not be set on a JOB=1 entry and
C     must be unchanged since the previous entry if JOB>1.
C     For J = 1,..., N, IPTRU(J) holds the position in
C     FACT and IRNF of the end of the packed part of column J of U.
C W is an array of length M used as workspace for holding
C      the expanded form of a sparse vector.
C IW is an unsigned integer array of length M+2*N used as workspace.
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1) A negative value will indicate an error return and a
C       positive value a warning. Possible nonzero values are:
C      -1  M < 1 or N < 1.
C      -2  NE < 0.
C      -3  Insufficient space.
C      -4  There are duplicated entries.
C      -5  JOB < 1, 3 when ICNTL(6)=0, or > 3.
C      -6  JOB = 2, but entries were dropped in the corresponding JOB=1
C          entry.
C      -7  NP < 0 or NP > N.
C     -(7+K) Pivot too small in column K when JOB=2.
C      +1  Rank deficient.
C    INFO(4) Minimum storage required to factorize matrix (or
C            recommended value for LFACT if INFO(1) = -3.
C    INFO(5) Computed rank of the matrix.
C    INFO(6) Number of entries dropped from the data structure.
C    INFO(7) Number of rows processed in full storage.
C RINFO need not be set on entry. On exit, it holds the number of
C    floating-pounsigned int operations performed.
\end{verbatim}

\subsubsection{MA50CD}

@o stackC.h -d
@{
//void useMA50CD();
@}

\begin{verbatim}

      SUBROUTINE MA50CD(M,N,ICNTL,IQ,NP,TRANS,LFACT,FACT,IRNF,IPTRL,
     +                  IPTRU,B,X,W,INFO)
C MA50C/CD uses the factorization produced by
C     MA50B/BD to solve A x = b or (A trans) x = b.
C
      UNSIGNED INTEGER M,N,ICNTL(7),IQ(N),NP
      LOGICAL TRANS
      UNSIGNED INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      UNSIGNED INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION B(*),X(*),W(*)
      UNSIGNED INTEGER INFO(7)
C
C M  is an unsigned integer variable set to the number of rows.
C     It is not altered by the subroutine.
C N  is an unsigned integer variable set to the number of columns.
C     It is not altered by the subroutine.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 2, plus all parameters on entry and exit.
C     ICNTL(5) must be set to control the level of BLAS used:
C       0 Level 1 BLAS.
C      >0 Level 2 BLAS.
C IQ is an unsigned integer array holding the permutation Q.
C     It is not altered by the subroutine.
C NP is an unsigned integer variable that must be unchanged since calling
C     MA50B/BD. It holds the number of rows and columns in packed
C     storage. It is not altered by the subroutine.
C TRANS a logical variable thatmust be set to .TRUE. if (A trans)x = b
C     is to be solved and to .FALSE. if A x = b is to be solved.
C     TRANS is not altered by the subroutine.
C LFACT is an unsigned integer variable set to the size of FACT and IRNF.
C     It is not altered by the subroutine.
C FACT is an array that must be unchanged since calling MA50B/BD. It
C     holds the packed part of L/U by columns, and the full part of L/U
C     by columns. U has unit diagonal entries, which are not stored, and
C     the signs of the off-diagonal entries are inverted.  In the packed
C     part, the entries of U precede the entries of L; also the diagonal
C     entries of L head each column of L and are reciprocated.
C     FACT is not altered by the subroutine.
C IRNF is an unsigned integer array that must be unchanged since calling
C     MA50B/BD. It holds the row numbers of the packed part of L/U, and
C     the row numbers of the full part of L/U.
C     It is not altered by the subroutine.
C IPTRL is an unsigned integer array that must be unchanged since calling
C     MA50B/BD. For J = 1,..., NP, IPTRL(J) holds the position in
C     FACT and IRNF of the end of column J of L.
C     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
C     It is not altered by the subroutine.
C IPTRU is an unsigned integer array that must be unchanged since calling
C     MA50B/BD. For J = 1,..., N, IPTRU(J) holds the position in
C     FACT and IRNF of the end of the packed part of column J of U.
C     It is not altered by the subroutine.
C B is an array that must be set to the vector b.
C     It is not altered.
C X is an array that need not be set on entry. On return, it holds the
C    solution x.
C W is a work array of length max(M,N).
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1) A nonzero value will indicate an error return. Possible
C      nonzero values are:
C      -1  M < 1 or N < 1
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




%\bibliographystyle{authordate2}
\bibliographystyle{plain}
\bibliography{files,anderson}
\end{document}


