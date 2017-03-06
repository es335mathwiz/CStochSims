%	$Id: stackC.w,v 1.39 2000/03/28 21:05:12 m1gsa00 Exp m1gsa00 $	
\documentclass{article}
%\documentclass[html]{article}
@i /msu/home/m1gsa00/miscProjects/math/miscLatexPkg.tex
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
#include <float.h>
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


@o stack.h -d
@{
void nxtCDmats(@<nxtCDmats argument list@>);
@}
Function uses SPARSEKIT's CSR format.
Since HARWELL functions expects CSC,
MA50CD operates on the transpose.

@d nxtCDmats definition
@{
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
int * numberOfEquations, int * lagss, int * leadss,
int * rowDim,
int * maxNumberHElements,
double * smatsA,int * smatsJA,int * smatsIA,
double * fvecA,int * fvecJA,int * fvecIA,
double ** cmatsA,int **cmatsJA,int **cmatsIA,
double **dmatsA,int **dmatsJA,int **dmatsIA
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
submat_(rowDim,aOne,aOne,rowDim,firstColumn,lastColumn,
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
amub_(rowDim,cColumns,aOne,ao,jao,iao,
(cmatsA[timeOffset]),(cmatsJA[timeOffset]),(cmatsIA[timeOffset]),
b,jb,ib,nzmax,iw,ierr);
pathNewtAssert(*ierr == 0);
bump((cmatsIA[timeOffset])[*rowDim]-(cmatsIA[timeOffset])[0]);
aSmallDouble=DBL_EPSILON;
filter_(rowDim,aOne,&aSmallDouble,b,jb,ib,b,jb,ib,nzmax,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);
/*actually want to subtract so mult elements by -1 also need to shift right*/
for(j=0;j<ib[*rowDim]-1;j++)
{b[j]=(-1)*b[j];jb[j]=jb[j]+(*numberOfEquations*(timeOffset+*lagss+1));};

aplb_(rowDim,cColumns,aOne,oddSumCA,oddSumCJA,oddSumCIA,
b,jb,ib,evenSumCA,evenSumCJA,evenSumCIA,
nzmax,iw,ierr);
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

amub_(rowDim,aOne,aOne,ao,jao,iao,
(dmatsA[timeOffset]),(dmatsJA[timeOffset]),(dmatsIA[timeOffset]),
b,jb,ib,nzmax,iw,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);
aSmallDouble=DBL_EPSILON;
filter_(rowDim,aOne,&aSmallDouble,b,jb,ib,b,jb,ib,nzmax,ierr);
pathNewtAssert(*ierr == 0);
bump(ib[*rowDim]-ib[0]);

/*actually want to subtract so mult elements by -1*/
for(j=0;j<ib[*rowDim]-1;j++)b[j]=(-1)*b[j];

aplb_(rowDim,aOne,aOne,oddSumDA,oddSumDJA,oddSumDIA,
b,jb,ib,
evenSumDA,evenSumDJA,evenSumDIA,
nzmax,iw,ierr);
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
submat_(rowDim,aOne,aOne,rowDim,
firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,nr,nc,
oddSumCA,oddSumCJA,oddSumCIA);
*nonZeroNow=oddSumCIA[*rowDim]-oddSumCIA[0];

ma50id_(cntl,icntl);
ma50ad_(rowDim,rowDim,nonZeroNow,nzmax,oddSumCA,oddSumCJA,jcn,oddSumCIA,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
/*wordybump(info[3]);*/
pathNewtAssert(info[0]>=0);
/* restore odd since ad is destructive*/
submat_(rowDim,aOne,aOne,rowDim,
firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,nr,nc,
oddSumCA,oddSumCJA,jcn);

@}


@d factorize matrix
@{
ma50bd_(rowDim,rowDim,nonZeroNow,aOne,
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

submat_(rowDim,aOne,
aOne,rowDim,firstColumn,lastColumn,
evenSumCA,evenSumCJA,evenSumCIA,
nr,nc,
b,jb,ib);

csrdns_(rowDim,aOne,b,jb,ib,
nsSumC,rowDim,ierr);

ma50cd_(rowDim,rowDim,icntl,oddSumCIA,np,trans,
lfact,fact,irnf,iptrl,iptru,
nsSumC,x,
w,info);
bump(info[3]);
pathNewtAssert(info[0]>=0);

*nzmaxLeft= *nzmax-cmatsExtent-1;
dnscsr_(aOne,rowDim,nzmaxLeft,x,
aOne,tb+(itb[i]-1),jtb+(itb[i]-1),itb+i,
ierr);
pathNewtAssert(*ierr == 0);

itb[i+1]=itb[i+1]+cmatsExtent;
itb[i]=itb[i]+cmatsExtent;
cmatsExtent=itb[i+1]-1;
}
bump(cmatsExtent);

aSmallDouble=DBL_EPSILON;
filter_(cColumns,aOne,&aSmallDouble,tb,jtb,itb,tb,jtb,itb,nzmax,ierr);
csrcsc_(cColumns,aOne,aOne,tb,jtb,itb,cmatsA[0],cmatsJA[0],cmatsIA[0]);

/*expand sum of d's*/
csrdns_(rowDim,aOne,oddSumDA,oddSumDJA,oddSumDIA,nsSumD,rowDim,ierr);
pathNewtAssert(*ierr == 0);
bump(*rowDim);
/*code should use info from previous call to set lfact
also can avoid calls to ma50ad once pattern settles down*/

ma50cd_(rowDim,rowDim,icntl,oddSumCIA,np,
trans,lfact,fact,irnf,iptrl,iptru,nsSumD,x,w,info);

dnscsr_(rowDim,aOne,rowDim,x,
rowDim,dmatsA[0],dmatsJA[0],dmatsIA[0],ierr);
/*wordybump(info[3]);*/
pathNewtAssert(*ierr == 0);
@}

\subsection{nxtCDmats Variable Declarations}
\label{sec:nxtcddeclarations}


@d nxtCDmats variable declarations
@{
int maxElementsEncountered=0;
/*void * calloc(unsigned amt,int size);*/
double * evenSumCA;int * evenSumCJA;int * evenSumCIA;
double * evenSumDA;int * evenSumDJA;int * evenSumDIA;
double * oddSumCA;int * oddSumCJA;int * oddSumCIA;
double * oddSumDA;int * oddSumDJA;int * oddSumDIA;
double  *ao;int *jao,*iao;
double *b;int *jb,*ib;
double *tb;int *jtb,*itb;
double *tmp;int *jtmp,*itmp;
int *firstColumn,*lastColumn,*nr,*nc;
int *iw,*ierr,*nzmax,*nonZeroNow;int cmatsExtent;int *nzmaxLeft;
int i;
int j;
int timeOffset;
int *jcn;
double * cntl;
int * icntl;
int * ip ;
int * np;
int * jfirst;
int * lenr;
int * lastr;
int * nextr;
int * ifirst;
int * lenc;
int * lastc;
int * nextc;
int * info;
double * rinfo;
int *lfact;
double * fact;
int *irnf;
int * iptrl;
int * iptru;
int * aOne;
double aSmallDouble;
double * w;
double * x;
int * trans;
int * cColumns;
int * balColumns;
int *hColumns;
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
firstColumn = (int *)calloc(1,sizeof(int));
lastColumn = (int *)calloc(1,sizeof(int));
nr = (int *)calloc(1,sizeof(int));
nc = (int *)calloc(1,sizeof(int));
ierr = (int *)calloc(1,sizeof(int));
nzmax = (int *)calloc(1,sizeof(int));
nzmaxLeft = (int *)calloc(1,sizeof(int));
hColumns = (int *)calloc(1,sizeof(int));
nonZeroNow = (int *)calloc(1,sizeof(int));
cColumns = (int *)calloc(1,sizeof(int));
balColumns = (int *)calloc(1,sizeof(int));
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
jcn = (int *)calloc(*maxNumberHElements,sizeof(int));
cntl= (double *)calloc(5,sizeof(double));
icntl= (int *)calloc(9,sizeof(int));
ip = (int *)calloc(*rowDim,sizeof(int));
np = (int *)calloc(1,sizeof(int));
jfirst = (int *)calloc(*rowDim,sizeof(int));
lenr = (int *)calloc(*rowDim,sizeof(int));
lastr = (int *)calloc(*rowDim,sizeof(int));
nextr = (int *)calloc(*rowDim,sizeof(int));
ifirst = (int *)calloc(*rowDim,sizeof(int));
lenc = (int *)calloc(*rowDim,sizeof(int));
lastc = (int *)calloc(*rowDim,sizeof(int));
nextc = (int *)calloc(*rowDim,sizeof(int));
info = (int *)calloc(7,sizeof(int));
rinfo = (double *)calloc(1,sizeof(double));
/* ma50bd*/
lfact =(int *)calloc(1,sizeof(int));
*lfact = ( *maxNumberHElements);/*pessimistic setting for filling*/
fact = (double *)calloc(*lfact,sizeof(double));
irnf = (int *)calloc(*lfact,sizeof(int));
iptrl = (int *)calloc(*rowDim,sizeof(int));
iptru = (int *)calloc(*rowDim,sizeof(int));
x = (double *)calloc(*rowDim * *numberOfEquations * (*leadss+1),sizeof(double));
trans = (int *) calloc(1,sizeof(int));
nsSumD = (double *)calloc(*rowDim,sizeof(double));
nsSumC = (double *)calloc(*rowDim /** *numberOfEquations * (*leadss+1+ *lagss)*/,sizeof(double));
aOne = (int *)calloc(1,sizeof(int));
*aOne = 1;
@}


@d nxtCDmats array variable allocations
@{

ib = (int *)calloc(*rowDim * *leadss + 1,sizeof(int));
jb = (int *)calloc(*maxNumberHElements,sizeof(int));
b = (double *)calloc(*maxNumberHElements,sizeof(double));
itb = (int *)calloc(*cColumns + 1,sizeof(int));
jtb = (int *)calloc(*maxNumberHElements,sizeof(int));
tb = (double *)calloc(*maxNumberHElements,sizeof(double));


evenSumCIA = (int *)calloc(( *rowDim * *leadss)+1,sizeof(int));
evenSumCJA = (int *)calloc(*maxNumberHElements,sizeof(int));
evenSumCA = (double *)calloc(*maxNumberHElements,sizeof(double));


/*larger than necessary now so that can use for transpose in csrcsc */
oddSumCIA = (int *)calloc(( *rowDim * *leadss)+ 1,sizeof(int));
oddSumCJA = (int *)calloc(*maxNumberHElements,sizeof(int));
oddSumCA = (double *)calloc(*maxNumberHElements,sizeof(double));


evenSumDIA = (int *)calloc(*rowDim+1,sizeof(int));/*MLK*/
evenSumDJA = (int *)calloc(*rowDim,sizeof(int));
evenSumDA = (double *)calloc(*rowDim,sizeof(double));


oddSumDIA = (int *)calloc(*rowDim+1,sizeof(int));/*MLK*/
oddSumDJA = (int *)calloc(*rowDim,sizeof(int));
oddSumDA = (double *)calloc(*rowDim,sizeof(double));


iao = (int *)calloc(*rowDim+1,sizeof(int));
jao = (int *)calloc(*maxNumberHElements,sizeof(int));
ao = (double *)calloc(*maxNumberHElements,sizeof(double));

/*work array needs elements equal to number of columns of matrix*/
iw = (int *)calloc(*rowDim * (*leadss + *lagss + 3),sizeof(int));
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
cfree(firstColumn);
cfree(lastColumn);
cfree(nr);
cfree(nc);
cfree(ierr);
cfree(nzmax);
cfree(nzmaxLeft);
cfree(nonZeroNow);
cfree(hColumns);
cfree(cColumns);
cfree(balColumns);
@}


@d nxtCDmats array variable deallocations
@{
cfree(ib);
cfree(jb);
cfree(b);
cfree(itb);
cfree(jtb);
cfree(tb);
cfree(evenSumCIA);
cfree(evenSumCJA);
cfree(evenSumCA);
cfree(oddSumCIA);
cfree(oddSumCJA);
cfree(oddSumCA);
cfree(evenSumDIA);
cfree(evenSumDJA);
cfree(evenSumDA);
cfree(oddSumDIA);
cfree(oddSumDJA);
cfree(oddSumDA);
cfree(iao);
cfree(jao);
cfree(ao);
cfree(iw);
cfree(w);
@}

@d ma50xx array variable deallocations
@{

/*for ma50ad*/
cfree(jcn);
cfree(cntl);
cfree(icntl);
cfree(ip);
cfree(np);
cfree(jfirst);
cfree(lenr);
cfree(lastr);
cfree(nextr);
cfree(ifirst);
cfree(lenc);
cfree(lastc);
cfree(nextc);
cfree(info);
cfree(rinfo);
/* ma50bd*/
cfree(lfact);
cfree(fact);
cfree(irnf);
cfree(iptrl);
cfree(iptru);
cfree(x);
cfree(trans);
cfree(nsSumD);
cfree(nsSumC);
cfree(aOne);
@}


\subsection{backSub Definition}
\label{sec:backSub}
@d oneStepBack argument list
@{
int * rowDim,
double ** yvecA,int ** yvecJA,int ** yvecIA,
double ** cmatsA,int **cmatsJA,int **cmatsIA,
double **dmatsA,int **dmatsJA,int **dmatsIA
 @| numberOfEquations lagss leadss 
smatsA smatsJA smatsIA 
yvecA yvecJA yvecIA 
cmatsA cmatsJA cmatsIA
dmatsA dmatsJA dmatsIA
@}



\subsection{oneStepBack Definition}
\label{sec:oneStepBack}

@o stack.h -d
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
  amub_(rowDim,aOne,aOne,cmatsA[0],cmatsJA[0],cmatsIA[0],
  yvecA[0+1 ],yvecJA[0+1 ],yvecIA[0+1 ],
  rcy,rcyj,rcyi,rowDim,iw,ierr);
pathNewtAssert(*ierr == 0);

aSmallDouble=DBL_EPSILON;

filter_(rowDim,aOne,&aSmallDouble,rcy,rcyj,rcyi,rcy,rcyj,rcyi,nzmax,ierr);
pathNewtAssert(*ierr == 0);

  for(i=0;i<rcyi[*rowDim]-rcyi[0];i++)rcy[i]=(-1)*rcy[i];
  aplb_(rowDim,aOne,aOne,
  dmatsA[0],dmatsJA[0],dmatsIA[0],
  rcy,rcyj,rcyi,
  yvecA[(0) ],yvecJA[(0) ],yvecIA[(0)],
  nzmax,iw,ierr);
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
int * aOne;
double *rcy;
int *rcyj ,*rcyi;
int *ierr;
int i;
int *nzmax;
int * iw;
@}
@d oneStepBack variable allocations
@{
nzmax=(int *)calloc(1,sizeof(int));
aOne=(int *)calloc(1,sizeof(int));
ierr=(int *)calloc(1,sizeof(int));
*aOne=1;
*nzmax=*rowDim;
rcy= (double *)calloc(*rowDim,sizeof(double));
rcyj=(int *) calloc(*rowDim,sizeof(int));
rcyi=(int *)calloc(*rowDim+1,sizeof(int));
iw = (int *)calloc(1,sizeof(int));
@}
@d oneStepBack variable deallocations
@{
cfree(ierr);
cfree(aOne);
cfree(rcy);
cfree(rcyj);
cfree(rcyi);
cfree(iw);
cfree(nzmax);
@}

\subsection{terminalConditions}
\label{sec:term}
@d computeR argument list
@{
int * numberOfEquations,int * lags, int * leads,
int * maxNumberElements,
void (* func)(),void (* dfunc)(),double * params,
double * expansionPoint,
double * qMat,int * qMatj,int * qMati,
double * impact, int * impactj, int * impacti,
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
amux_(&rowDim,expansionPoint,rvec,termConstr,termConstrj,termConstri);
amux_(&rowDim,expansionPoint,impactPart1+(numberOfEuations*lags),
impactr,impactrj,impactri);
amux_(&rowDim,impactPart1,impactPart2,
impactr,impactrj,impactri);
for(i=0;i<numberOfEquations*leads;i++){rvec[i]=rvec[i]-impactPart2[i];
}
cfree(impactPart1);
cfree(impactPart2);
}
@}




\subsection{nxtGuess Definition}
\label{sec:nxtGuess}


@o stack.h -d
@{
void nxtGuess(@<nxtGuess argument list@>);
void newNxtGuess(@<nxtGuess argument list@>);
@}


@d nxtGuess argument list
@{
int * numberOfEquations,int * lags, int * leads, int * capT,
double **fmats,int  **fmatsj,int  **fmatsi,
double **smats,int  **smatsj,int  **smatsi,
int *maxNumberHElements,
double * termConstr,int * termConstrj,int * termConstri,double * fp,
double * initialX,
double * shockVec,
double * updateDirection @| termConstr fp initialX shockVec
theFunc theDrvFunc capT 
@}

@d chkDrv definition
@{
#include <math.h>
#define NEGLIGIBLEDOUBLE 1.0e-9
chkDrv(int n, double * fdrv,int * fdrvj,int * fdrvi,
double * fvec,double * delxvec)
{
int i;
int aOne=1;
double * fvals;
fvals = (double * ) calloc(n,sizeof(double));
#ifdef DEBUG 
printf("chkDrv:beginning\n");
#endif


amux_(&n,delxvec,fvals,fdrv,fdrvj,fdrvi);
for(i=0;i<n;i++){
if(fabs(fvals[i]-fvec[i])>NEGLIGIBLEDOUBLE){
#ifdef DEBUG 
printf("chkDrv:discrepancy for %d,(%e,%e)\n",i,fvals[i],fvec[i]);
#endif
}
}
cfree(fvals);
#ifdef DEBUG 
printf("chkDrv:done\n");
#endif

}
#define addOneToFEvals (intOutputInfo[2])++
#define homotopyAlpha (doubleControlParameters+10)
#define addOneToFDrvEvals (intOutputInfo[3])++
#define maxNumberStackElements intControlParameters[5]
constructFdrv(int numberOfEquations,int lags, int leads,int pathLength,
double * xvec,double * params,void (* vFunc)(),void (* vFuncDrv)(),
double * termConstr,int * termConstrj,int * termConstri,
double * fixedPoint,double * intercept,double * linearizationPoint,int * exogRows, int * exogCols, int * exogenizeQ,
double * shockVec,
double * fvec,
double * fdrv,int * fdrvj,int * fdrvi,int ihomotopy,
int * intControlParameters,double * doubleControlParameters,
int * intOutputInfo, double * doubleOutputInfo)
{
double * deviations;
int * ignore;double  dignore[1]={1.0};
int rowDim;int * fvecj;int * fveci;
int i;int j;int soFar;double * zeroShockVec;
ignore = (int *)calloc(numberOfEquations+1,sizeof(int));
fvecj = (int *)calloc(numberOfEquations+1,sizeof(int));
fveci = (int *)calloc(numberOfEquations+1,sizeof(int));
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
amux_(&rowDim,deviations,fvec+(numberOfEquations*(lags+pathLength)),
termConstr,termConstrj,termConstri);
for(i=0;i<numberOfEquations*leads;i++){fvec[numberOfEquations*(lags+pathLength)+i]=fvec[numberOfEquations*(lags+pathLength)+i]-intercept[i];}
cfree(ignore);
cfree(fvecj);
cfree(fveci);
cfree(deviations);
cfree(zeroShockVec);
}
@}

@d newNxtGuess argument list
@{
int * sysDim,
int * maxNumberHElements,
double * fvec,
double * fdrv,int * fdrvj,int * fdrvi,
double * xdel,
int * ma50bdJob,
int * ma50bdIq,
double * ma50bdFact,
int * ma50bdIrnf,
int * ma50bdIptrl,
int * ma50bdIptru,
int * intControlParameters,double * doubleControlParameters
@}

@d newNxtGuess definition
@{
#define sysDimSwitchLevel 30000
#define ma50DropThreshold (1e-8)

void newNxtGuess(@<newNxtGuess argument list@>)
{


int i;
int maxElementsEncountered=0;
double * copychkfdrv;int * copychkfdrvj;int * copychkfdrvi;
int * jcn;
double * cntl;
int * icntl;
int * ip ;
int * np;
int * jfirst;
int * lenr;
int * lastr;
int * nextr;
int * ifirst;
int * lenc;
int * lastc;
int * nextc;
int * info;
double * rinfo;
int *lfact;
double * fact;
int *irnf;
int * iptrl;
int * iptru;
int *iw;
double * w;
int nzmax;
int  * aOne;
int nonZeroNow;
int  trans;
@}
@d newNxtGuess definition
@{

copychkfdrv = (double * )calloc(*maxNumberHElements,sizeof(double));
copychkfdrvj = (int * )calloc(*maxNumberHElements,sizeof(int));
copychkfdrvi=(int * ) calloc(*sysDim+1,sizeof(int));
jcn = (int *)calloc(*maxNumberHElements,sizeof(int));
cntl= (double *)calloc(5,sizeof(double));
icntl= (int *)calloc(9,sizeof(int));
ip = (int *)calloc(*sysDim,sizeof(int));
np = (int *)calloc(1,sizeof(int));
jfirst = (int *)calloc(*sysDim,sizeof(int));
lenr = (int *)calloc(*sysDim,sizeof(int));
lastr = (int *)calloc(*sysDim,sizeof(int));
nextr = (int *)calloc(*sysDim,sizeof(int));
w = (double *)calloc(*sysDim,sizeof(double));
iw = (int *)calloc(3**sysDim,sizeof(int));
ifirst = (int *)calloc(*sysDim,sizeof(int));
lenc = (int *)calloc(*sysDim,sizeof(int));
lastc = (int *)calloc(*sysDim,sizeof(int));
nextc = (int *)calloc(*sysDim,sizeof(int));
info = (int *)calloc(7,sizeof(int));
rinfo = (double *)calloc(1,sizeof(double));
lfact = (int *)calloc(1,sizeof(int));
*lfact = ( *maxNumberHElements);/*pessimistic setting for filling*/
aOne = (int *)calloc(1,sizeof(int));
*aOne=1;
@}
@d newNxtGuess definition
@{
copmat_(sysDim,fdrv,fdrvj,fdrvi,
copychkfdrv,copychkfdrvj,copychkfdrvi,aOne,aOne);

ma50id_(cntl,icntl);
cntl[1]=ma50Balance;
cntl[2]=ma50DropEntry;
cntl[3]=ma50DropCol;
icntl[3]=ma50PivotSearch;
/*if(*sysDim>sysDimSwitchLevel){cntl[2]=ma50DropThreshold;}*/
nzmax=*maxNumberHElements;
nonZeroNow=copychkfdrvi[*sysDim]-copychkfdrvi[0];
ma50ad_(sysDim,sysDim,&nonZeroNow,
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
ma50bd_(sysDim,sysDim,&nonZeroNow,ma50bdJob,
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
ma50bd_(sysDim,sysDim,&nonZeroNow,ma50bdJob,
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
ma50bd_(sysDim,sysDim,&nonZeroNow,ma50bdJob,
fdrv,fdrvj,fdrvi,
cntl,icntl,ip,copychkfdrvi,
np,lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
w,iw,info,rinfo);
wordybump(info[3]);
pathNewtAssert(info[0]>=0);
}
}
    trans = 1;
ma50cd_(sysDim,sysDim,icntl,copychkfdrvi,np,&trans,
lfact,ma50bdFact,ma50bdIrnf,ma50bdIptrl,ma50bdIptru,
fvec,xdel,
w,info);
*maxNumberHElements=maxElementsEncountered;
pathNewtAssert(info[0]>=0);
@}
@d newNxtGuess definition
@{
cfree(aOne);
cfree(jcn );
cfree(cntl);
cfree(icntl);
cfree(ip );
cfree(np );
cfree(jfirst );
cfree(lenr );
cfree(lastr );
cfree(nextr );
cfree(w);
cfree(iw);
cfree(ifirst );
cfree(lenc );
cfree(lastc );
cfree(nextc );
cfree(info );
cfree(rinfo );
cfree(lfact);
cfree(copychkfdrv);cfree(copychkfdrvj);cfree(copychkfdrvi);
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
int maxElementsEncountered=0;
int * ma50bdJob;
int * ma50bdIq;
double * ma50bdFact;
int * ma50bdIrnf;
int * ma50bdIptrl;
int * ma50bdIptru;

int tNow;
double * deviations;
double **cmats;int  **cmatsj;int  **cmatsi;
double **dmats;int **dmatsj;int **dmatsi;
double **ymats;int  **ymatsj;int  **ymatsi;
double *gmats;int  *gmatsj;int  *gmatsi;
int i,j;
int *hColumns;
int *qColumns;
int *rowDim;
double *fullfvec;
double *fulldfvec;
double *fullXvec;
int * aOne;
int * ierr;
@}


@d nxtGuess variable initializations
@{
aOne=(int *)calloc(1,sizeof(int));
*aOne=1;
ierr=(int *)calloc(1,sizeof(int));
deviations=(double *)calloc(*numberOfEquations*(*lags+*leads),sizeof(double));
hColumns=(int *)calloc(1,sizeof(int));
qColumns=(int *)calloc(1,sizeof(int));
rowDim=(int *)calloc(1,sizeof(int));
*hColumns=*numberOfEquations*(*lags+1+*leads);
*qColumns=*numberOfEquations*(*lags+*leads);
*rowDim=*numberOfEquations* (*leads? *leads:1);

fullfvec=(double *)calloc(*numberOfEquations**leads,sizeof(double));
fulldfvec=(double *)calloc(*numberOfEquations*(*lags+*leads+ 1),sizeof(double));
fullXvec = (double *)calloc(*numberOfEquations *  (*lags+*leads+*capT),sizeof(double));
for(i=0;i<(*lags+*leads+*capT) * *numberOfEquations ;i++){fullXvec[i]=1.0;};

cmats =(double **)calloc(*capT+(*lags+*leads)+1,sizeof(double *));
cmatsj =(int **)calloc(*capT+(*lags+*leads)+1,sizeof(int *));
cmatsi =(int **)calloc(*capT+(*lags+*leads)+1,sizeof(int *));
dmats =(double **)calloc(*capT+(*lags+*leads)+1,sizeof(double *));
dmatsj =(int **)calloc(*capT+(*lags+*leads)+1,sizeof(int *));
dmatsi =(int **)calloc(*capT+(*lags+*leads)+1,sizeof(int *));
gmats =(double *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(double));
gmatsj =(int *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(int));
gmatsi =(int *)calloc(*numberOfEquations*(*leads?*leads:1)+1,sizeof(int));
ymats =(double **)calloc(*capT+(*leads+*lags)+1,sizeof(double *));
ymatsj =(int **)calloc(*capT+(*leads+*lags)+1,sizeof(int *));
ymatsi =(int **)calloc(*capT+(*leads+*lags)+1,sizeof(int *));
for(i=0;i<*capT+(*lags+*leads)+1;i++){
ymats[i] =(double *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(double));
ymatsj[i] =(int *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(int));
ymatsi[i] =(int *)calloc(*numberOfEquations*(*leads?*leads:1)+1,sizeof(int));
cmats[i] =(double *)calloc(*maxNumberHElements,sizeof(double));
cmatsj[i] =(int *)calloc(*maxNumberHElements,sizeof(int));
cmatsi[i] =(int *)calloc(*numberOfEquations*(*leads?*leads:1)+1,sizeof(int));
dmats[i] =(double *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(double));
dmatsj[i] =(int *)calloc(*numberOfEquations*(*leads?*leads:1),sizeof(int));
dmatsi[i] =(int *)calloc(*numberOfEquations*(*leads?*leads:1)+1,sizeof(int));
}

@}
@d nxtGuess storage allocations
@{
ma50bdIptru = (int *)calloc(*numberOfEquations* (*leads?*leads:1),sizeof(int));
ma50bdIptrl = (int *)calloc(*numberOfEquations* (*leads?*leads:1),sizeof(int));
ma50bdIrnf = (int *)calloc(*maxNumberHElements,sizeof(int));
ma50bdFact = (double *)calloc(*maxNumberHElements,sizeof(double));
ma50bdIq = (int *)calloc(*maxNumberHElements,sizeof(int));
ma50bdJob = (int *)calloc(1,sizeof(int));
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
csrdns_(numberOfEquations,aOne,ymats[i+*lags],ymatsj[i+*lags],ymatsi[i+*lags],
updateDirection+((*lags + i) * *numberOfEquations),numberOfEquations,ierr);
pathNewtAssert(*ierr == 0);
bump(ymatsi[i+*lags][*numberOfEquations]-ymatsi[i+*lags][0]);
}
if(*leads>0){
csrdns_(rowDim,aOne,ymats[*capT+*lags],ymatsj[*capT+*lags],ymatsi[*capT+*lags],
updateDirection+((*lags + *capT) * *numberOfEquations),rowDim,ierr);
pathNewtAssert(*ierr == 0);
bump(ymatsi[*capT+*lags][*numberOfEquations]-ymatsi[*capT+*lags][0]);
}

@}
@d nxtGuess storage deallocations
@{
cfree(ma50bdIptru);
cfree(ma50bdIptrl);
cfree(ma50bdIrnf);
cfree(ma50bdFact);
cfree(ma50bdIq);
cfree(ma50bdJob);
cfree(deviations);
cfree(hColumns);
cfree(qColumns);
cfree(rowDim);
cfree(fullXvec);
cfree(fullfvec);
cfree(fulldfvec);
cfree(ierr);
cfree(aOne);

for(i=0;i<*capT+(*lags+*leads)+1;i++){
cfree(cmats[i]);
cfree(cmatsj[i]);
cfree(cmatsi[i]);
cfree(dmats[i]);
cfree(dmatsj[i]);
cfree(dmatsi[i]);
cfree(ymats[i]);
cfree(ymatsj[i]);
cfree(ymatsi[i]);
}
cfree(cmats);
cfree(cmatsj);
cfree(cmatsi);
cfree(dmats);
cfree(dmatsj);
cfree(dmatsi);
cfree(gmats);
cfree(gmatsj);
cfree(gmatsi);
cfree(ymats);
cfree(ymatsj);
cfree(ymatsi);

@}

\subsection{nxtFPGuess Definition}
\label{sec:nxtFPGuess}


@o stack.h -d
@{
void nxtFPGuess(@<nxtFPGuess argument list@>);
@}



@d nxtFPGuess argument list
@{
int * numberOfEquations,int * lags, int * leads,
double * fmats,int * fmatsj,int * fmatsi,
double * smats,int * smatsj,int * smatsi,
int *maxNumberHElements,
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
int * ma50bdJob;
int * ma50bdIq;
double * ma50bdFact;
int * ma50bdIrnf;
int * ma50bdIptrl;
int * ma50bdIptru;

int j;
double **cmats;int  **cmatsj;int  **cmatsi;
double **dmats;int **dmatsj;int **dmatsi;
double **ymats;int  **ymatsj;int  **ymatsi;
double *gmats;int  *gmatsj;int  *gmatsi;
int i;
int *hColumns;
int *qColumns;
int *rowDim;
double *fullfvec;
double *fulldfvec;
double *fullXvec;
int * aOne;
int * ierr;
@}


@d nxtFPGuess variable initializations
@{
aOne=(int *)calloc(1,sizeof(int));
*aOne=1;
ierr=(int *)calloc(1,sizeof(int));
hColumns=(int *)calloc(1,sizeof(int));
qColumns=(int *)calloc(1,sizeof(int));
rowDim=(int *)calloc(1,sizeof(int));
*hColumns=*numberOfEquations*(*lags+1+*leads);
*qColumns=*numberOfEquations*(*lags+*leads);
*rowDim=*numberOfEquations* *leads;

fullfvec=(double *)calloc(*numberOfEquations**leads,sizeof(double));
fulldfvec=(double *)calloc(*numberOfEquations*(*lags+*leads+ 1),sizeof(double));

cmats =(double **)calloc((*lags+*leads)+2,sizeof(double *));
cmatsj =(int **)calloc((*lags+*leads)+2,sizeof(int *));
cmatsi =(int **)calloc((*lags+*leads)+2,sizeof(int *));
dmats =(double **)calloc((*lags+*leads)+2,sizeof(double *));
dmatsj =(int **)calloc((*lags+*leads)+2,sizeof(int *));
dmatsi =(int **)calloc((*lags+*leads)+2,sizeof(int *));
gmats =(double *)calloc(*numberOfEquations,sizeof(double));
gmatsj =(int *)calloc(*numberOfEquations,sizeof(int));
gmatsi =(int *)calloc(*numberOfEquations,sizeof(int));
ymats =(double **)calloc((*lags+*leads)+2,sizeof(double *));
ymatsj =(int **)calloc((*lags+*leads)+2,sizeof(int *));
ymatsi =(int **)calloc((*lags+*leads)+2,sizeof(int *));
for(i=0;i<(*lags+*leads)+2;i++){
ymats[i] =(double *)calloc(*numberOfEquations* *leads,sizeof(double));
ymatsj[i] =(int *)calloc(*numberOfEquations* *leads,sizeof(int));
ymatsi[i] =(int *)calloc(*numberOfEquations* *leads+1,sizeof(int));
cmats[i] =(double *)calloc(*maxNumberHElements,sizeof(double));
cmatsj[i] =(int *)calloc(*maxNumberHElements,sizeof(int));
cmatsi[i] =(int *)calloc((*numberOfEquations  * *leads) +1,sizeof(int));
dmats[i] =(double *)calloc(*numberOfEquations* *leads,sizeof(double));
dmatsj[i] =(int *)calloc(*numberOfEquations* *leads,sizeof(int));
dmatsi[i] =(int *)calloc(*numberOfEquations* *leads+1,sizeof(int));
}

@}
@d nxtFPGuess storage allocations
@{
ma50bdIptru = (int *)calloc(*numberOfEquations* *leads,sizeof(int));
ma50bdIptrl = (int *)calloc(*numberOfEquations* *leads,sizeof(int));
ma50bdIrnf = (int *)calloc(*maxNumberHElements,sizeof(int));
ma50bdFact = (double *)calloc(*maxNumberHElements,sizeof(double));
ma50bdIq = (int *)calloc(*maxNumberHElements,sizeof(int));
ma50bdJob = (int *)calloc(1,sizeof(int));
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
dnscsr_(numberOfEquations,aOne,numberOfEquations,fullfvec,
numberOfEquations,
fmats,fmatsj,fmatsi,
ierr);
pathNewtAssert(*ierr == 0);
bump(fmatsi[*numberOfEquations]-fmatsi[0]);
dnscsr_(numberOfEquations,hColumns,maxNumberHElements,
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
csrdns_(numberOfEquations,aOne,ymats[i],ymatsj[i],ymatsi[i],
updateDirection+(( i) * *numberOfEquations),numberOfEquations,ierr);
pathNewtAssert(*ierr == 0);


}
csrdns_(rowDim,aOne,ymats[*lags+1],ymatsj[*lags+1],ymatsi[*lags+1],
updateDirection+((*lags+1 ) * *numberOfEquations),rowDim,ierr);
pathNewtAssert(*ierr == 0);




@}
@d nxtFPGuess storage deallocations
@{
cfree(ma50bdIptru);
cfree(ma50bdIptrl);
cfree(ma50bdIrnf);
cfree(ma50bdFact);
cfree(ma50bdIq);
cfree(ma50bdJob);
cfree(hColumns);
cfree(qColumns);
cfree(rowDim);

cfree(fullfvec);
cfree(fulldfvec);
cfree(ierr);
cfree(aOne);
for(i=0;i<(*lags+*leads)+2;i++){
cfree(cmats[i]);
cfree(cmatsj[i]);
cfree(cmatsi[i]);
cfree(dmats[i]);
cfree(dmatsj[i]);
cfree(dmatsi[i]);
cfree(ymats[i]);
cfree(ymatsj[i]);
cfree(ymatsi[i]);
}
cfree(cmats);
cfree(cmatsj);
cfree(cmatsi);
cfree(dmats);
cfree(dmatsj);
cfree(dmatsi);
cfree(gmats);
cfree(gmatsj);
cfree(gmatsi);
cfree(ymats);
cfree(ymatsj);
cfree(ymatsi);

@}



\subsection{compPathError Definition}
\label{sec:compPathError}


@d compPathError definition
@{
compPathError(int * numberOfEquations,int * lags,int * leads,
void (* theFunction)(),
double * termConstr,int * termConstrj,int * termConstri,double * fp,
double * initialX,
double * shockVec,
int * capT,
int * maxNumberHElements,
double ** fmats,int ** fmatsj,int **fmatsi,
double ** smats,int ** smatsj,int **smatsi)
{
int maxElementsEncountered=0;
int * rowDim;
int * qColumns;
int tNow;
int * aOne;
int * aZero;
int * ierr;int i;
double * deviations;
double * zeroShockVec;
double * fullfvec;
ierr = (int *) calloc(1,sizeof(int));
aOne= (int *) calloc(1,sizeof(int));
*aOne=1;
aZero= (int *) calloc(1,sizeof(int));
*aZero=1;
rowDim= (int *) calloc(1,sizeof(int));
qColumns= (int *) calloc(1,sizeof(int));
*rowDim=*numberOfEquations* *leads;
*qColumns=*numberOfEquations* (*leads+*lags);
deviations = (double *) calloc(*numberOfEquations* (*lags+*leads+1),sizeof(double));
fullfvec = (double *) calloc(*numberOfEquations,sizeof(double));

for(tNow=0;tNow<*capT;tNow++) {
 theFunction(initialX+((tNow/*+(*lags-1)*/) * *numberOfEquations),initialX,
shockVec,
  fmats[tNow],fmatsj[tNow],fmatsi[tNow]);

}


copmat_(rowDim,termConstr,termConstrj,termConstri,
smats[*capT],smatsj[*capT],smatsi[*capT],&aOne,&aOne);

/*xxxxxxxxx add code for deviations using gmat*/
for(i=0;i<*numberOfEquations* (*lags+ *leads);i++){
deviations[i]=initialX[*numberOfEquations* *capT+i]-fp[i];}
amux_(rowDim,deviations,fullfvec,smats[*capT],smatsj[*capT],smatsi[*capT]);
dnscsr_(rowDim,aOne,rowDim,fullfvec,
aOne,
fmats[*capT],fmatsj[*capT],fmatsi[*capT],ierr);
}
@}


\subsection{B matrix computation}
\label{sec:bmat}

@d applySparseReducedForm argument list
@{
int rowDim,int colDim,
double * initialX, double * fp,
double * bmat,int * bmatj,int * bmati,
double * resultX
@}
@d applySparseReducedForm declarations
@{
double * deviations;
int i;
@}
@d applySparseReducedForm allocations
@{
deviations = (double *) calloc(colDim,sizeof(double));
@}
@d applySparseReducedForm frees
@{
cfree(deviations);
@}
@d applySparseReducedForm
@{
void applySparseReducedForm(@<applySparseReducedForm argument list@>)
{
@<applySparseReducedForm declarations@>
@<applySparseReducedForm allocations@>
for(i=0;i<colDim;i++){deviations[i]=initialX[i]-fp[i%rowDim];}
amux_(&rowDim,deviations,resultX,bmat,bmatj,bmati);
for(i=0;i<rowDim;i++){resultX[i]=resultX[i]+fp[i%rowDim];}
@<applySparseReducedForm frees@>
}
@}

@d obtainSparseReducedForm argument list
@{
int * maxNumberHElements, 
int qrows, int qcols, double * qmat, int * qmatj, int * qmati,
 double * bmat, int * bmatj, int * bmati
@}

@d obtainSparseReducedForm allocations
@{

jcn = (int *)calloc(*maxNumberHElements,sizeof(int));
cntl= (double *)calloc(5,sizeof(double));
icntl= (int *)calloc(9,sizeof(int));
ip = (int *)calloc(qrows,sizeof(int));
np = (int *)calloc(1,sizeof(int));
jfirst = (int *)calloc(qrows,sizeof(int));
lenr = (int *)calloc(qrows,sizeof(int));
lastr = (int *)calloc(qrows,sizeof(int));
nextr = (int *)calloc(qrows,sizeof(int));
w = (double *)calloc(qrows,sizeof(double));
iw = (int *)calloc(3*qrows,sizeof(int));
ifirst = (int *)calloc(qrows,sizeof(int));
lenc = (int *)calloc(qrows,sizeof(int));
lastc = (int *)calloc(qrows,sizeof(int));
nextc = (int *)calloc(qrows,sizeof(int));
info = (int *)calloc(7,sizeof(int));
rinfo = (double *)calloc(1,sizeof(double));
/* ma50bd*/

qrmat = (double *) calloc(*maxNumberHElements,sizeof(double));
qrmatj = (int *) calloc(*maxNumberHElements,sizeof(double));
qrmati = (int *) calloc(qrows+1,sizeof(double));
tb = (double *) calloc(*maxNumberHElements,sizeof(double));
jtb = (int *) calloc(*maxNumberHElements,sizeof(double));
itb = (int *) calloc(qcols-qrows+1,sizeof(double));
b = (double *) calloc(*maxNumberHElements,sizeof(double));
jb = (int *) calloc(*maxNumberHElements,sizeof(double));
ib = (int *) calloc(qcols-qrows+1,sizeof(double));


lfact =(int *)calloc(1,sizeof(int));
*lfact = ( *maxNumberHElements);/*pessimistic setting for filling*/
fact = (double *)calloc(*lfact,sizeof(double));
irnf = (int *)calloc(*lfact,sizeof(int));
iptrl = (int *)calloc(qrows,sizeof(int));
iptru = (int *)calloc(qrows,sizeof(int));
x = (double *)calloc(  qcols,sizeof(double));
nsSumC = (double *)calloc(qrows ,sizeof(double));
@}
@d obtainSparseReducedForm frees
@{
cfree(w);
cfree(iw);
cfree(b);
cfree(jb);
cfree(ib);
cfree(tb);
cfree(jtb);
cfree(itb);
cfree(jcn );
cfree(cntl);
cfree(icntl);
cfree(ip );
cfree(np );
cfree(jfirst );
cfree(lenr );
cfree(lastr );
cfree(nextr );
cfree(ifirst );
cfree(lenc );
cfree(lastc );
cfree(nextc );
cfree(info );
cfree(rinfo );
cfree(/* ma50bd*/

qrmat );
cfree(qrmatj );
cfree(qrmati );
cfree(lfact );
cfree(fact );
cfree(irnf );
cfree(iptrl );
cfree(iptru );
cfree(x );
cfree(nsSumC );
@}


@d obtainSparseReducedForm declarations
@{
void * calloc(unsigned amt,int size);
double * nsSumC;int ierr;double * x;
int nzmaxLeft;double aSmallDouble;double nsSumD;
int cmatsExtent;int i;int cColumns;
double *b;int *jb,*ib;
double *tb;int *jtb,*itb;
int  trans;
 double * qrmat; int * qrmatj; int * qrmati;
int *iw;double * w;
int  aOne; int  firstColumn;int  lastColumn;
int  nr;int  nc;int nonZeroNow;int nzmax;
int * jcn;
double * cntl;
int * icntl;
int * ip ;
int * np;
int * jfirst;
int * lenr;
int * lastr;
int * nextr;
int * ifirst;
int * lenc;
int * lastc;
int * nextc;
int * info;
double * rinfo;
int *lfact;
double * fact;
int *irnf;
int * iptrl;
int * iptru;
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
submat_(&qrows,&aOne,&aOne,&qrows,
&firstColumn,&lastColumn,
qmat,qmatj,qmati,&nr,&nc,
qrmat,qrmatj,qrmati);

nonZeroNow=qrmati[qrows]-qrmati[0];

ma50id_(cntl,icntl);
nzmax=*maxNumberHElements;
ma50ad_(&qrows,&qrows,&nonZeroNow,
&nzmax,qrmat,qrmatj,jcn,qrmati,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
/* restore odd since ad is destructive*/
submat_(&qrows,&aOne,&aOne,&qrows,
&firstColumn,&lastColumn,
qmat,qmatj,qmati,&nr,&nc,
qrmat,qrmatj,jcn);

@}


@d bmat factorize matrix
@{
ma50bd_(&qrows,&qrows,&nonZeroNow,&aOne,
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

submat_(&qrows,&aOne,
&aOne,&qrows,&firstColumn,&lastColumn,
qmat,qmatj,qmati,
&nr,&nc,
b,jb,ib);

csrdns_(&qrows,&aOne,b,jb,ib,
nsSumC,&qrows,&ierr);

ma50cd_(&qrows,&qrows,icntl,qrmati,np,&trans,
lfact,fact,irnf,iptrl,iptru,
nsSumC,x,
w,info);

nzmaxLeft= nzmax-cmatsExtent-1;
dnscsr_(&aOne,&qrows,&nzmaxLeft,x,
&aOne,tb+(itb[i]-1),jtb+(itb[i]-1),itb+i,
&ierr);
itb[i+1]=itb[i+1]+cmatsExtent;
itb[i]=itb[i]+cmatsExtent;
cmatsExtent=itb[i+1]-1;
}
aSmallDouble=DBL_EPSILON;
filter_(&cColumns,&aOne,&aSmallDouble,tb,jtb,itb,tb,jtb,itb,&nzmax,&ierr);
csrcsc_(&cColumns,&aOne,&aOne,tb,jtb,itb,bmat,bmatj,bmati);
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

OTHERSOURCE =  /msu/res2/m1gsa00/aim/frbus/miscLatexPkg.tex

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

pureStackView:   testStack
      purify -windows=no -log-file=stackViewPurify f77 -o pureStack testStack.o stackC.o testStackModel.o ma50ad.o $(SPARSELIB) $(LINKFLAGS)  $(AIMLIB)
      pureStack

@}

\section{Numerical Recipes Modifications of FPnewt}


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
  {printf("pathNewtAssert: processid=%ld\n",getpid());\
   printf ("%s:%u: failed assertion\n", file, lineno);\
   	kill(getpid(),SIGUSR2);}


@}

@o myNewt.c -d 
@{
@<FPnewt defines@>
void FPnewt(int * numberOfEquations,int * lags, int * leads,
void (* func)(),void (* dfunc)(),double * params,
double x[],
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,
int *check)
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
#include "sparseAMA.h"
#define NRANSI
#include "/msu/res2/m1gsa00/aim/frbus/nrutil.h"
#define MAXITS 200
#define TOLF 1.0e-12
#define TOLMIN 1.0e-6
#define TOLX 1.0e-10
#define STPMX 100.0
#define GAMMA 1.0

int nn;
double *fvec;
#define FREERETURN {cfree(fvec);cfree(xold);cfree(shockVec);\
	cfree(p);cfree(g);cfree(aOne);cfree(ierr);\
	cfree(indx);return;}
#define PFREERETURN {cfree(fvec);cfree(xold);cfree(xoldls);cfree(xdel);\
cfree(deviations);cfree(fullfvec);\
	cfree(p);cfree(g);cfree(aOne);cfree(ierr);\
	cfree(indx);cfree(rowDim);cfree(qColumns);return;}
@}

@d FPnewt declarations
@{
    int n;double * xdel;
/*	double fmin(double x[]);*/
	void lnsrch(int n,int np,int reps,
    double xold[], double   * fold, double g[], double p[],double * params,
		 double * shockVec,double * f, double stpmax, int *check, 
         void (*func)(double*,double*,double*,double*,int*,int*),double *x);
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int i,its,j,*indx,*aOne,*ndns,*ierr;
	double d,den,f,fold,stpmax,sum,temp,test,*g,*p,*xold,normSum;
    double * shockVec;
    n=*numberOfEquations*(*lags+*leads+1);
	aOne=(int *)calloc(1,sizeof(int));
    *aOne=1;
	ierr=(int *)calloc(1,sizeof(int));
	indx=(int *)calloc(n+1,sizeof(int));
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
        csrdns_(numberOfEquations,
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
        csrdns_(numberOfEquations,
        aOne,fmats[0],fmatsj[0],fmatsi[0],fvec+(*numberOfEquations**lags),
        numberOfEquations,ierr);

        amux_(numberOfEquations,fvec,
        g+(*numberOfEquations * *lags),
        smats[0],smatsj[0],smatsi[0]);
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
void pathNewt(int * numberOfEquations,int * lags, int * leads,int * pathLength,
void (* vecfunc)(),void (* fdjac)(),double * params,double * shockVec,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,double * qMat,int * qMatj,int * qMati,
double * fixedPoint,
double x[],
int *check)
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

copmat_(rowDim,qMat,qMatj,qMati,
smats[*pathLength],smatsj[*pathLength],smatsi[*pathLength],
aOne,aOne);


/*xxxxxxxxx add code for deviations using gmat*/
for(i=0;i<*numberOfEquations* (*lags+ *leads);i++){
deviations[i]=x[*numberOfEquations* *pathLength+i]-fixedPoint[i];}
amux_(rowDim,deviations,fullfvec,smats[*pathLength],smatsj[*pathLength],smatsi[*pathLength]);
dnscsr_(rowDim,aOne,rowDim,fullfvec,
aOne,
fmats[*pathLength],fmatsj[*pathLength],fmatsi[*pathLength],ierr);
for(i=0;i<*rowDim;i++){
  fvec[(*pathLength * *numberOfEquations)+i]=fullfvec[i];
            normSum += SQR(fullfvec[i]);}

        f= 0.5 * normSum;
for(tNow=0;tNow<*pathLength;tNow++) {
        csrdns_(numberOfEquations,
        aOne,fmats[tNow],fmatsj[tNow],fmatsi[tNow],
        fvec+(tNow * *numberOfEquations),numberOfEquations,ierr);
        }
/* modification end */
@}



@d Anderson-Moore algorithm  variable deallocations
@{
cfree(asymptoticLinearization);
cfree(cond);
cfree(epsi);
cfree(inform);
cfree(iq);
cfree(itsbad);
cfree(nbig);
cfree(nexa);
cfree(nnum);
cfree(nroot);
cfree(rooti);
cfree(rootr);
cfree(uprbnd);
cfree(qColumns);
cfree(hColumns);
@}


@d Anderson-Moore algorithm  variable allocations
@{
qColumns=(int *)calloc(1,sizeof(int));
*qColumns=(*numberOfEquations)*((*lags)+(*leads));
hColumns=(int *)calloc(1,sizeof(int));
*hColumns=(*numberOfEquations)*((*lags)+(*leads)+1);
asymptoticLinearization=(double *)
   calloc( (*numberOfEquations)*((*numberOfEquations)*((*leads)+(*lags)+1)),
   sizeof(double));

cond = (double *)calloc(1,sizeof(double));
epsi = (double *)calloc(1,sizeof(double));
inform = (int *)calloc(1,sizeof(int));
iq = (int *)calloc(1,sizeof(int));
itsbad = (int *)calloc(1,sizeof(int));
nbig = (int *)calloc(1,sizeof(int));
nexa = (int *)calloc(1,sizeof(int));
nnum = (int *)calloc(1,sizeof(int));
nroot = (int *)calloc(1,sizeof(int));
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
void * aPointerToVoid;/*adding since all the sparseAMA.h files have this arg*/
printf("from stackC.w line 2194\n");
cPrintSparse(*numberOfEquations,smats[0],smatsj[0],smatsi[0]);
sparseAMA(&maxHElementsForSparseAMA,DISCRETE_TIME,*numberOfEquations,
*numberOfEquations*(*lags+1+*leads),*leads,
smats[0],smatsj[0],smatsi[0],
newH,newHj,newHi,
&auxInit,&qRows,
qMat,qMatj,qMati,
&essential,rootr,rooti,&returnCode,
aPointerToVoid);
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


@o myNewt.c -d
@{

computeAsymptoticQMatrix(
int * numberOfEquations,int * lags, int * leads,
void (* func)(),void (* dfunc)(),double * params,
double canadaFP[],int * pthLngth,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,
double * AMqMatrix,
int * ierr
)
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
int * inform;
int * iq;
int * itsbad;
int * nbig;
int * nexa;
int * nnum;
int * nroot;
double * zeroVector;
int * qColumns;
double * rootr;
double * rooti;
double * uprbnd;
double * wts;
double * err;
int * hColumns;
int i;
@}

@o myNewt.c -d
@{

altComputeAsymptoticQMatrix(
int * numberOfEquations,int * lags, int * leads,
void (* func)(),void (* dfunc)(),double * params,
double canadaFP[],int * pthLngth,
double ** fmats, int ** fmatsj, int ** fmatsi,
double ** smats, int ** smatsj, int ** smatsi,
int * maxNumberElements,
double * qMat,int * qMatj,int * qMati,
int * ierr
)
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
cfree(newH);
cfree(newHj);
cfree(newHi);
cfree(rootr);
cfree(rooti);
@}
@d altComputeAsymptoticQMatrix variable allocations
@{
newH=(double *)calloc(maxHElements,sizeof(double));
newHj=(int *)calloc(maxHElements,sizeof(int));
newHi=(int *)calloc(maxHElements,sizeof(int));
rootr=(double *) calloc((*numberOfEquations)*((*lags)+(*leads)),
sizeof(double));
rooti=(double *) calloc((*numberOfEquations)*((*lags)+(*leads)),
sizeof(double));

@}
@d altComputeAsymptoticQMatrix variable declarations
@{
int auxInit;
int qRows;
int maxHElements=50000;
int maxHElementsForSparseAMA;
int essential;
int returnCode=0;
double * zeroVector;
int * qColumns;
double * newH;
int * newHj;
int * newHi;
double * rootr;
double * rooti;
double * uprbnd;
double * wts;
double * err;
int * hColumns;
int i;
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
	rowDim=(int *)calloc(1,sizeof(int));
	*rowDim=*numberOfEquations**leads;
    aOne=(int *)calloc(1,sizeof(int));
    *aOne=1;
    qColumns=(int *)calloc(1,sizeof(int));
    *qColumns=*numberOfEquations*(*leads+*lags);
    deviations=(double *)calloc(*qColumns,sizeof(double));
    fullfvec=(double *)calloc(*rowDim,sizeof(double));
	ierr=(int *)calloc(1,sizeof(int));
	indx=(int *)calloc(n+1,sizeof(int));
	g=(double *)calloc(n+1,sizeof(double));
	p=(double *)calloc(n+1,sizeof(double));
	xdel=(double *)calloc(n+1,sizeof(double));
	xold=(double *)calloc(n+1,sizeof(double));
	xoldls=(double *)calloc(n+1,sizeof(double));
	fvec=(double *)calloc(n+1,sizeof(double));


	nn=n;

@}
@d pathNewt declarations
@{    int  n; 
    int tNow;int * rowDim;int * qColumns;
    double * deviations;double * fullfvec;
/*	double fmin(double x[]);*/
	void lnsrch(int  n, int np,int reps,double xold[], double  *fold, double g[], double p[],double * params,double * shockVec,
		 double * f, double stpmax, int *check, void (*func)(),double * x);


	int i,its,j,*indx,*aOne,*ndns,*ierr;
	double d,den,f,fold,stpmax,sum,temp,test,*g,*p,*xold,*xoldls,*xdel,normSum;
@}


\subsection{lnsrch}
\label{sec:lnsrch}

@o myNewt.c -d
@{
#include <math.h>
#define NRANSI
#include "/msu/res2/m1gsa00/aim/frbus/nrutil.h"
#define ALF 1.0e-4
#define TOLX 1.0e-10

extern void dgemm_(char * noTransp1,char * noTransp2,
int * neq,int * aOne,int *BMatrixColumns1,
double *aFloatOne,
double *__asymptoticBmatrices,int *BMatrixRows1,
double *deviationsFromPeriodicPath,int *BMatrixColumns2,
double *aZero,
double *pathNow,
int *BMatrixRows2);


extern double  dnrm2_();
#include <string.h>

void lnsrch(int  n,int np,int reps,
double xold[], double * fold, double g[], double p[], 
		 double * params,double * shockVec,double * f,double stpmax, int *check,
			void (*func)(double*,double*,double*,double*,int*,int*),double * x)
{
@< lnsrch declarations@>
@< lnsrch preloop@>
@< lnsrch loop@>
@< lnsrch postloop@>
cfree(aOne);cfree(aZero);cfree(aTwo);cfree(fvec);cfree(fvecj);cfree(fveci);
}
#undef ALF
#undef TOLX
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software '>9m_L31.. */

@}

@d lnsrch declarations
@{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;
	int * aOne;int * aTwo,*ierr,*aZero;
	double * fnow,*fvec,*fnorm, *xorig,*aDoubleZero;
    int * fvecj;
    double *dir,*aDoubleOne;
    char * transp, *noTransp;
	int * fveci;
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
	ierr=(int *)calloc(1,sizeof(int));
	fveci=(int *)calloc(n+1,sizeof(int));
	aOne=(int *)calloc(1,sizeof(int));
	*aOne=1;
	aZero=(int *)calloc(1,sizeof(int));
	*aOne=0;
	aTwo=(int *)calloc(1,sizeof(int));
	*aTwo=2;
	fvec=(double *)calloc(n,sizeof(double));
	xorig=(double *)calloc(n,sizeof(double));
	fvecj=(int *)calloc(n,sizeof(double));

@}
@d lnsrch preloop
@{


/*
    fold = dnrm2_(&np,foldp,aOne);
*/

*fold=0;
		for(i=0;i<reps;i++){
((*func)(xold+(i*np),params,shockVec,fvec,fvecj,fveci));
        cnrms_(&np,aZero,fvec,fvecj,fveci,xorig);
        *fold += xorig[0];
        }
      *fold *= 0.5;

/*
csrdns_(&n,aOne,fvec,fvecj,fveci,xorig,&n,ierr);

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

			cfree(transp);
            cfree(ierr);
			cfree(noTransp);
			cfree(dir);
			cfree(aDoubleZero);
			cfree(aDoubleOne);
			cfree(fnow);
			cfree(xorig);
			cfree(fvec);
			cfree(fvecj);
			cfree(fveci);
			cfree(aOne);
			cfree(aZero);
			cfree(aTwo);
@}

@d lnsrch loop
@{
	for (;;) {
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
		for(i=0;i<1;i++)((*func)(x+(i*np),params,shockVec,fvec,fvecj,fveci));
        cnrms_(&np,aTwo,fvec,fvecj,fveci,xorig);*f = xorig[0] *  0.5;
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



@o notStack.c -d
@{
@<backSub definition@>
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



@d cPrintMatrix definition
@{
void cPrintMatrix(nrows,ncols,matrix)
int  nrows;
int  ncols;
double * matrix;
{
int i,j;
for(i=0;i<nrows;i++)
for(j=0;j<ncols;j++)printf("[%d] [%d] %f\n",i,j,matrix[i+(j*nrows)]);
}
void cPrintMatrixNonZero(nrows,ncols,matrix,zerotol)
int  nrows;
int  ncols;
double * matrix;
double zerotol;
{
int i,j;
double fabs();
for(i=0;i<nrows;i++)
for(j=0;j<ncols;j++)
    if(fabs(matrix[i+(j*nrows)]) > zerotol)
    printf("[%d] [%d] %f\n",i,j,matrix[i+(j*nrows)]);
}

void cPrintSparse(int rows,double * a,int * aj,int * ai)
{
int i,j,numEls;
numEls=ai[rows]-ai[0];
printf("matrix has %d non zero element/s\n",numEls);
for(i=0;i<rows;i++){
for(j=ai[i];j<ai[i+1];j++){
printf("row=%d,col=%d,val=%f\n",i+1,aj[j-1],a[j-1]);}}
}

@}
\subsection{SPARSEKIT2}
\label{sec:sparsekit2}

\subsubsection{dump}
@o stack.h -d
@{
/*void dump_();*/
@}


\begin{verbatim}
c----------------------------------------------------------------------- 
      subroutine dump (i1,i2,values,a,ja,ia,iout)
      integer i1, i2, ia(*), ja(*), iout
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
c i1    = first row (column) to print out
c i2    = last row (column) to print out 
c values= logical. indicates whether or not to print real values.
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
@o stack.h -d
@{
void amub_();
@}

\begin{verbatim}

/*
c----------------------------------------------------------------------c
       subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                  c,jc,ic,nzmax,iw,ierr) 
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix by matrix product C = A B 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A = row dimension of C
c ncol  = integer. The column dimension of B = column dimension of C
c job   = integer. Job indicator. When job = 0, only the structure
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
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic    = resulting matrix C in compressed sparse row sparse format.
c           
c ierr  = integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = integer work array of length equal to the number of
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

@o stack.h -d
@{
void amux_();
@}


\begin{verbatim}

/*
c----------------------------------------------------------------------c
      subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
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
@o stack.h -d
@{
void aplb_();
@}


\begin{verbatim}

/*
c-----------------------------------------------------------------------
      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix sum  C = A+B. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
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
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic    = resulting matrix C in compressed sparse row sparse format.
c       
c ierr  = integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = integer work array of length equal to the number of
c         columns in A.
c
c-----------------------------------------------------------------------
*/
\end{verbatim}


\subsubsection{csrcsc}

@o stack.h -d
@{
void csrcsc2_();
@}

\begin{verbatim}
c-----------------------------------------------------------------------
      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
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
c job   = integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
c     for any other normal usage, enter ipos=1.
c a = real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja    = integer array of length nnz containing the column positions
c     of the corresponding elements in a.
c ia    = integer of size n+1. ia(k) contains the position in a, ja of
c     the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao    = real array of size nzz containing the "a" part of the transpose
c jao   = integer array of size nnz containing the column indices.
c iao   = integer array of size n+1 containing the "ia" index array of
c     the transpose. 
c
c----------------------------------------------------------------------- 
\end{verbatim}
@o stack.h -d
@{
void csrcsc_();
@}


\subsubsection{csrcsc}
\begin{verbatim}
c----------------------------------------------------------------------- 
      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
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
c job	= integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 

\end{verbatim}

\subsubsection{csrdns}
@o stack.h -d
@{
void csrdns_();
@}


\begin{verbatim}

/*c----------------------------------------------------------------------c
      subroutine csrdns(nrow,ncol,a,ja,ia,dns,ndns,ierr) 
      real*8 dns(ndns,*),a(*)
      integer ja(*),ia(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row    to    Dense 
c-----------------------------------------------------------------------
c
c converts a row-stored sparse matrix into a densely stored one
c
c On entry:
c---------- 
c
c nrow  = row-dimension of a
c ncol  = column dimension of a
c a, 
c ja, 
c ia    = input matrix in compressed sparse row format. 
c         (a=value array, ja=column array, ia=pointer array)
c dns   = array where to store dense matrix
c ndns  = first dimension of array dns 
c
c on return: 
c----------- 
c dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
c 
c ierr  = integer error indicator. 
c         ierr .eq. 0  means normal return
c         ierr .eq. i  means that the code has stopped when processing
c         row number i, because it found a column number .gt. ncol.
c 
c----------------------------------------------------------------------- 
*/
\end{verbatim}


\subsubsection{dnscsr}
@o stack.h -d
@{
void dnscsr_();
@}


\begin{verbatim}

/*
c----------------------------------------------------------------------- 
      subroutine dnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)
      real*8 dns(ndns,*),a(*)
      integer ia(*),ja(*)
c-----------------------------------------------------------------------
c Dense     to    Compressed Row Sparse 
c----------------------------------------------------------------------- 
c
c converts a densely stored matrix into a row orientied
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
c a, ja, ia = value, column, pointer  arrays for output matrix 
c
c ierr  = integer error indicator: 
c         ierr .eq. 0 means normal retur
c         ierr .eq. i means that the the code stopped while
c         processing row number i, because there was no space left in
c         a, and ja (as defined by parameter nzmax).
c----------------------------------------------------------------------- 

*/
\end{verbatim}


\subsubsection{submat}
@o stack.h -d
@{
void submat_();
@}

\begin{verbatim}

/*
      subroutine submat (n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)
      integer n,job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
      real*8 a(*),ao(*) 
c-----------------------------------------------------------------------
c extracts the submatrix A(i1:i2,j1:j2) and puts the result in 
c matrix ao,iao,jao
c---- In place: ao,jao,iao may be the same as a,ja,ia.
c-------------- 
c on input
c---------
c n = row dimension of the matrix 
c i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
c          extracted. 
c j1,j2 = two integers with j2 .ge. j1 indicating the range of columns 
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
c   the column indices,and iao being the pointer to the beginning 
c   of the row,in arrays a,ja.
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
*/
\end{verbatim}



\subsection{HARWELL Routines}
\label{sec:harwell}


\subsubsection{MA50AD}

@o stack.h -d
@{
void ma50ad_();
@}

\begin{verbatim}
      SUBROUTINE MA50AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP,JFIRST,
     +                  LENR,LASTR,NEXTR,IW,IFIRST,LENC,LASTC,NEXTC,
     +                  INFO,RINFO)

C MA50A/AD chooses a pivot sequence using a Markowitz criterion with
C     threshold pivoting.

C If  the user requires a more convenient data interface then the MA48
C     package should be used. The MA48 subroutines call the MA50
C     subroutines after checking the user's input data and optionally
C     permute the matrix to block triangular form.

      INTEGER M,N,NE,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),JCN(LA),IQ(N)
      DOUBLE PRECISION CNTL(4)
      INTEGER ICNTL(7),IP(M),NP,JFIRST(M),LENR(M),LASTR(M),NEXTR(M),
     +        IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(7)
      DOUBLE PRECISION RINFO

C M is an integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries
C      in the input matrix. It is not altered by the subroutine.
C LA is an integer variable that must be set to the size of A, IRN, and
C      JCN. It is not altered by the subroutine.
C A is an array that holds the input matrix on entry and is used as
C      workspace.
C IRN  is an integer array.  Entries 1 to NE must be set to the
C      row indices of the corresponding entries in A.  IRN is used
C      as workspace and holds the row indices of the reduced matrix.
C JCN  is an integer array that need not be set by the user. It is
C      used to hold the column indices of entries in the reduced
C      matrix.
C IQ is an integer array of length N. On entry, it holds pointers
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
C        point of view of rank.
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
C IP is an integer array of length M that need not be set on entry
C      unless ICNTL(7)=2 (see ICNTL(7) for details of this case).
C      During execution, IP(i) holds the position of the start of row i
C      of the reduced matrix or -IP(i) holds the row index in the
C      permuted matrix of row i. Before exit, IP(i) is made positive.
C NP is an integer variable. It need not be set on entry. On exit,
C     it will be set to the number of columns to be processed in
C     packed storage.
C JFIRST is an integer workarray of length M. JFIRST(i) is the
C      first column of the reduced matrix to have i entries or is
C      zero if no column has i entries.
C LENR is an integer workarray of length M that is used to hold the
C      numbers of entries in the rows of the reduced matrix.
C LASTR is an integer workarray of length M, used only if ICNTL(4) = 0.
C      For rows in the reduced matrix, LASTR(i) indicates the previous
C      row to i with the same number of entries. LASTR(i) is zero if
C      no such row exists.
C NEXTR is an integer workarray of length M, used only if ICNTL(4) = 0
C      or ICNTL(7)=2. If ICNTL(4)=0, for rows in the reduced matrix,
C      NEXTR(i) indicates the next row to i with the same number of
C      entries; and if row i is the last in the chain, NEXTR is
C      equal to zero. If ICNTL(7)=2, NEXTR is a copy of the value of
C      IP on entry.
C IW is an integer array of length M used as workspace and is used to
C     assist the detection of duplicate entries and the sparse SAXPY
C     operations. It is reset to zero each time round the main loop.
C IFIRST is an integer array of length N, used only if ICNTL(4) = 0
C      or ICNTL(7)=2. If ICNTL(4) = 0, it is a workarray; IFIRST(i)
C      points to the first row of the reduced matrix to have i entries
C      or is zero if no row has i entries. If ICNTL(7)=2, IFIRST
C      must be set on entry (see ICNTL(7) for details of this case).
C LENC is an integer workarray of length N that is used to hold
C      the numbers of entries in the columns of the reduced matrix.
C LASTC is an integer workarray of length N.  For columns in the reduced
C      matrix, LASTC(j) indicates the previous column to j with the same
C      number of entries.  If column j is the first in the chain,
C      LASTC(j) is equal to zero.
C NEXTC is an integer workarray of length N.  For columns in the reduced
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
C    floating-point operations needed for the factorization.


\end{verbatim}


\subsubsection{CNRMS}
\begin{verbatim}
c----------------------------------------------------------------------- 
      subroutine cnrms   (nrow, nrm, a, ja, ia, diag) 
      real*8 a(*), diag(nrow) 
      integer ja(*), ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each column of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
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

@o stack.h -d
@{
void ma50bd_();
@}

\begin{verbatim}
      
SUBROUTINE MA50BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP,IQ,NP,
     +                  LFACT,FACT,IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)
C MA50B/BD factorizes the matrix in AA/IRNA/IPTRA as P L U Q where
C     P and Q are permutations, L is lower triangular, and U is unit
C     upper triangular. The prior information that it uses depends on
C     the value of the parameter JOB.
C
      INTEGER M,N,NE,JOB
      DOUBLE PRECISION AA(NE)
      INTEGER IRNA(NE),IPTRA(N)
      DOUBLE PRECISION CNTL(4)
      INTEGER ICNTL(7),IP(M),IQ(N),NP,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION W(M)
      INTEGER IW(M+2*N),INFO(7)
      DOUBLE PRECISION RINFO
C
C M is an integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries
C      in the input matrix.  It is not altered by the subroutine.
C JOB is an integer variable that must be set to the value 1, 2, or 3.
C     If JOB is equal to 1 and any of the first NP recommended pivots
C      fails to satisfy the threshold pivot tolerance, the row is
C      interchanged with the earliest row in the recommended sequence
C      that does satisfy the tolerance. Normal row interchanges are
C      performed in the last N-NP columns.
C     If JOB is equal to 2, then M, N, NE, IRNA, IPTRA, IP, IQ,
C      LFACT, NP, IRNF, IPTRL, and IPTRU must be unchanged since a
C      JOB=1 entry for the same matrix pattern and no interchanges are
C      performed among the first NP pivots; if ICNTL(6) > 0, the first
C      N-ICNTL(6) columns of AA must also be unchanged.
C     If JOB is equal to 3, ICNTL(6) must be in the range 1 to N-1.
C      The effect is as for JOB=2 except that interchanges are
C      performed.
C     JOB is not altered by the subroutine.
C AA is an array that holds the entries of the matrix and
C      is not altered.
C IRNA is an integer array of length NE that must be set to hold the
C      row indices of the corresponding entries in AA. It is not
C      altered.
C IPTRA is an integer array that holds the positions of the starts of
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
C        point of view of rank.
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
C IP is an integer array. If JOB=1, it must be set so that IP(I) < IP(J)
C      if row I is recommended to precede row J in the pivot sequence.
C      If JOB>1, it need not be set. If JOB=1 or JOB=3, IP(I) is set
C      to -K when row I is chosen for pivot K and IP is eventually
C      reset to recommend the chosen pivot sequence to a subsequent
C      JOB=1 entry. If JOB=2, IP is not be referenced.
C IQ is an integer array that must be set so that either IQ(J) is the
C      column in position J in the pivot sequence, J=1,2,...,N,
C      or IQ(1)=0 and the columns are taken in natural order.
C      It is not altered by the subroutine.
C NP is an integer variable that holds the number of columns to be
C      processed in packed storage. It is not altered by the subroutine.
C LFACT is an integer variable set to the size of FACT and IRNF.
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
C IRNF is an integer array of length LFACT that need not be set on
C      a JOB=1 entry and must be unchanged since the previous entry
C      if JOB>1. On exit, IRNF(1) holds the number of dropped entries,
C      IRNF(2) holds the number of rows MF in full storage,
C      IRNF(3:IPTRL(N)) holds the row numbers of the packed part
C      of L/U, IRNF(IPTRL(N)+1:IPTRL(N)+MF) holds the row indices
C      of the full part of L/U, and IRNF(IPTRL(N)+MF+I), I=1,2,..,N-NP
C      holds the vector IPIV output by MA50GD.
C      If JOB=2, IRNF will be unaltered.
C IPTRL is an integer array that need not be set on a JOB=1 entry and
C     must be unchanged since the previous entry if JOB>1.
C     For J = 1,..., NP, IPTRL(J) holds the position in
C     FACT and IRNF of the end of column J of L.
C     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
C IPTRU is an integer array that need not be set on a JOB=1 entry and
C     must be unchanged since the previous entry if JOB>1.
C     For J = 1,..., N, IPTRU(J) holds the position in
C     FACT and IRNF of the end of the packed part of column J of U.
C W is an array of length M used as workspace for holding
C      the expanded form of a sparse vector.
C IW is an integer array of length M+2*N used as workspace.
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
C    floating-point operations performed.
\end{verbatim}

\subsubsection{MA50CD}

@o stack.h -d
@{
void ma50cd_();
@}

\begin{verbatim}

      SUBROUTINE MA50CD(M,N,ICNTL,IQ,NP,TRANS,LFACT,FACT,IRNF,IPTRL,
     +                  IPTRU,B,X,W,INFO)
C MA50C/CD uses the factorization produced by
C     MA50B/BD to solve A x = b or (A trans) x = b.
C
      INTEGER M,N,ICNTL(7),IQ(N),NP
      LOGICAL TRANS
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION B(*),X(*),W(*)
      INTEGER INFO(7)
C
C M  is an integer variable set to the number of rows.
C     It is not altered by the subroutine.
C N  is an integer variable set to the number of columns.
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
C IQ is an integer array holding the permutation Q.
C     It is not altered by the subroutine.
C NP is an integer variable that must be unchanged since calling
C     MA50B/BD. It holds the number of rows and columns in packed
C     storage. It is not altered by the subroutine.
C TRANS a logical variable thatmust be set to .TRUE. if (A trans)x = b
C     is to be solved and to .FALSE. if A x = b is to be solved.
C     TRANS is not altered by the subroutine.
C LFACT is an integer variable set to the size of FACT and IRNF.
C     It is not altered by the subroutine.
C FACT is an array that must be unchanged since calling MA50B/BD. It
C     holds the packed part of L/U by columns, and the full part of L/U
C     by columns. U has unit diagonal entries, which are not stored, and
C     the signs of the off-diagonal entries are inverted.  In the packed
C     part, the entries of U precede the entries of L; also the diagonal
C     entries of L head each column of L and are reciprocated.
C     FACT is not altered by the subroutine.
C IRNF is an integer array that must be unchanged since calling
C     MA50B/BD. It holds the row numbers of the packed part of L/U, and
C     the row numbers of the full part of L/U.
C     It is not altered by the subroutine.
C IPTRL is an integer array that must be unchanged since calling
C     MA50B/BD. For J = 1,..., NP, IPTRL(J) holds the position in
C     FACT and IRNF of the end of column J of L.
C     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
C     It is not altered by the subroutine.
C IPTRU is an integer array that must be unchanged since calling
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
\bibliographystyle{plainnat}
\bibliography{files,anderson}
\end{document}

