
/*declare prototypes for functions*/
float  dtime_(float * userSystemTime);
double atof();
#define rbcNEQS 4
#define rbcNLAGS 1
#define rbcNLEADS 1
#define PATHLENGTH 1000


FILE * outFile;
 /* char outFileName[250];*/
static char flnm[50] = "stochOut.m";
/*a counter*/
unsigned int i;/*int j;*/
/*modelDimensions call determines these*/
unsigned int  numberOfEquations[1]={rbcNEQS};
unsigned int  lags[1]={rbcNLAGS};
unsigned int  leads[1]={rbcNLEADS};
unsigned int pathLength[1]={PATHLENGTH};
unsigned int stochasticPathLength[1]={PATHLENGTH};
unsigned int t0[1]={0};
unsigned int tf[1]={0};
unsigned int replications[1]={1};
double totalTime[1];
double userSystemTime[2];
/*int shockIndex[1];*/
unsigned int  numberOfParameters;
/*int  numberOfDataValues;*/
unsigned int  numberOfShocks;
/*int  numberExog;*/


//int intControlParameters[widthIntControlInfo];
//double doubleControlParameters[widthDoubleControlInfo];
/* int intOutputInfo[REPLICATIONS*widthIntOutputInfo];
//double doubleOutputInfo[REPLICATIONS*widthDoubleOutputInfo];*/


/*int qRows=0;
 int auxInit=0;
 int aZero=0;*/

/*processCommandLine() determines defaults for these*/


/*timing routine variables*/


/*workspace*/
double **fmats;int  **fmatsj;int  **fmatsi;
double **smats;int  **smatsj;int  **smatsi;
/*success indicators for stochSims*/

/*csr q matrix*/
double * AMqMatrix;
int * AMqMatrixj;
int * AMqMatrixi;
/*csr b matrix*/
/*double * AMbMatrix;
int * AMbMatrixj;
int * AMbMatrixi;
 double * rootr;
 double * rooti;*/
/* double * brootr;
 double * brooti;
double*upsilonMatrix;int*upsilonMatrixj;int*upsilonMatrixi;
double*hMat;int*hMatj;int*hMati;
double*hzMat;int*hzMatj;int*hzMati;
double*cstar;int*cstarj;int*cstari;
double*phiInvMat;int*phiInvMatj;int*phiInvMati;
double*fmat;int*fmatj;int*fmati;
double*impact;int*impactj;int*impacti;*/
/*double*selectZmat;int*selectZmatj;int*selectZmati;*/
/*double*varthetaC;int*varthetaCj;int*varthetaCi;*/
/*double*varthetaZstar;int*varthetaZstarj;int*varthetaZstari;*/

 /*int * ma50bdJob;
 int * ma50bdIq;
 double * ma50bdFact;
 int * ma50bdIrnf;
 int * ma50bdIptrl;
 int * ma50bdIptru;*/
 /*int * cmpma50bdJob;
 int * cmpma50bdIq;
 double * cmpma50bdFact;
 int * cmpma50bdIrnf;
 int * cmpma50bdIptrl;
 int * cmpma50bdIptru;*/
/* int sysDim;*/



