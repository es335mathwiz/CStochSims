
#line 378 "stochProto.w"

void allocStochSim(unsigned int stochasticPathLength,unsigned int replications,unsigned int ** failedQ);
void freeStochSim(unsigned int ** failedQ);

#line 737 "stochProto.w"

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









#line 991 "stochProto.w"


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





#line 1063 "stochProto.w"



//double FMAX(double a,double b);
//double FMIN(double a,double b);
double FABS(double a);
double doRightSmaller(double a,double b);
double doSign(double a);





#line 1079 "stochProto.w"


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



void stochSim(
#line 361 "stochProto.w"
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
              unsigned int *failedQ
              
#line 1098 "stochProto.w"
);
void generateDraws(unsigned int t0Index,unsigned int tfIndex,unsigned int replications,unsigned int shocksAvailable,
unsigned int * iarray);

void processCommandLine(unsigned int argc, char * argv[],char ** namesArray,unsigned int modelNEQS,char ** paramNamesArray,unsigned int numberOfParameters,double * parameters,
double * dataValues,unsigned int numberDataValues,unsigned int numShockValues,
unsigned int * pathLength,unsigned int * replications,unsigned int * t0,unsigned int * stochasticPathLength,
unsigned int * intControlParameters,double* doubleControlParameters,char * flnm);

