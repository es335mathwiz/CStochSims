/*julliard example from paper describing stack*/
#include <math.h>
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

#define ey(t)     (stateVector[(t-(-1))*5+0])
#define pdot(t)     (stateVector[(t-(-1))*5+1])
#define rr(t)     (stateVector[(t-(-1))*5+2])
#define rs(t)     (stateVector[(t-(-1))*5+3])
#define y(t)     (stateVector[(t-(-1))*5+4])


#define eps (parameters[0])
#define g (parameters[1])





void julModPeriodicPointGuesser
(double * parameters,int period,
	double guessVector[7][5])
{
int i,j;
double svalue;
int timeOffset;
for(timeOffset=0;
	timeOffset<period+ 7 - 1;
			timeOffset++)
	{
guessVector[timeOffset][0]=0;
guessVector[timeOffset][1]=0;
guessVector[timeOffset][2]=0;
guessVector[timeOffset][3]=0;
guessVector[timeOffset][4]=0;
}
}



void julMod(double *stateVector,double *parameters,
double * shockVec,
double * aMat,int * jaMat,int *iaMat
)
{
int i;
double okay[2500];
okay[1]=pdot(0);
okay[2]=pdot(tMaOne);
okay[3]=pdot(tPaOne);
okay[4]=-g;
okay[5]=pow(g,2.);
okay[6]=y(0);
okay[7]=-okay[6];
okay[8]=y(tMaOne);
okay[9]=rr(0);
okay[10]=rs(0);
aMat[0]=ey(0);
aMat[1]=okay[1]-0.586*okay[2]-0.414*okay[3]-0.196*(okay[4]+okay[5]/(g+o\
kay[7]))-0.276*(okay[4]+okay[5]/(g-okay[8]))+eps*pdot(tPaFive);
aMat[2]=0.586*okay[2]+0.414*okay[3]+okay[9]-okay[10];
aMat[3]=-3.*okay[1]+okay[7]+okay[10];
aMat[4]=ey(tMaOne)+okay[6]-0.304*okay[8]+0.98*okay[9]+0.315*rr(tMaOne);
for(i=0;i<5;i++){aMat[i]=aMat[i]-shockVec[i];};
iaMat[0]=1.;
iaMat[1]=2.;
iaMat[2]=3.;
iaMat[3]=4.;
iaMat[4]=5.;
iaMat[5]=6.;
jaMat[0]=1.;
jaMat[1]=1.;
jaMat[2]=1.;
jaMat[3]=1.;
jaMat[4]=1.;
}

void julModDerivative(double *stateVector,double *parameters,
double * aMat,int * jaMat,int *iaMat
)
{
double okay[2500];
okay[1]=pow(g,2.);
aMat[0]=1.;
aMat[1]=-0.586;
aMat[2]=-0.276*okay[1]*pow(g-y(tMaOne),-2.);
aMat[3]=1.;
aMat[4]=-0.196*okay[1]*pow(g-y(0),-2.);
aMat[5]=-0.414;
aMat[6]=eps;
aMat[7]=0.586;
aMat[8]=1.;
aMat[9]=-1.;
aMat[10]=0.414;
aMat[11]=-3.;
aMat[12]=1.;
aMat[13]=-1.;
aMat[14]=1.;
aMat[15]=0.315;
aMat[16]=-0.304;
aMat[17]=0.98;
aMat[18]=1.;
iaMat[0]=1.;
iaMat[1]=2.;
iaMat[2]=8.;
iaMat[3]=12.;
iaMat[4]=15.;
iaMat[5]=20.;
jaMat[0]=6.;
jaMat[1]=2.;
jaMat[2]=5.;
jaMat[3]=7.;
jaMat[4]=10.;
jaMat[5]=12.;
jaMat[6]=32.;
jaMat[7]=2.;
jaMat[8]=8.;
jaMat[9]=9.;
jaMat[10]=12.;
jaMat[11]=7.;
jaMat[12]=9.;
jaMat[13]=10.;
jaMat[14]=1.;
jaMat[15]=3.;
jaMat[16]=5.;
jaMat[17]=8.;
jaMat[18]=10.;
}
/*made up data from normal dist*/
void julModData(int t,double * vectorOfVals)
{
int i;
static double theData[50][5]=
{-0.422401, 0.910338, -0.531898, -0.755108, -0.0329909, 0.148738, -0.182084, -0.88668, 2.69952, -0.120583, -0.21818, 0.0699103, -0.437539, -0.691171, 1.20786, 0.30338, -0.520588, 1.11671, 1.22693, -0.816738, 0.360885, -0.858474, 1.66715, -1.71189, 0.190519, 2.21293, 0.193504, 0.83057, -0.479637, -0.405844, -0.895405, 0.194228, -0.705833, -0.767932, -1.83828, 0.2259, -0.730635, -0.215615, 0.336557, -0.363111, -1.58305, -0.274715, -0.168059, 0.122736, -0.357303, -0.356944, 0.311806, -0.0803014, 0.0634804, 0.843273, 1.32363, -0.648111, -2.76966, 0.671527, -1.83073, -0.827917, -1.38832, 1.27912, 1.10909, -1.32497, -0.233533, -0.140782, 1.64871, -0.0797257, -0.241159, -1.1172, -0.884532, 1.08328, -0.984211, -0.730269, 3.12159, -0.905032, 0.724856, 0.238532, -0.583846, -0.919136, 0.7983, -0.760735, 0.653634, 1.68352, 0.295189, -0.895331, 0.168534, 0.00206746, -1.5253, -0.769304, 0.593886, -0.176777, -0.218989, 0.408153, -0.115861, -0.276384, 1.03232, -1.82229, 0.816341, -0.898487, 0.318052, 0.353191, -1.10885, -1.30401, 0.534568, -0.00540664, -1.27414, 0.417203, -1.05741, 0.658341, 0.233608, 0.49623, -0.877044, 0.187255, 0.572616, 0.739403, -0.343473, 0.54418, 1.08486, 1.45904, -1.15012, 1.45695, 2.2939, -0.650277, 1.01336, 1.63173, 1.39341, -0.0654548, -0.388803, 0.944093, -0.385784, 0.0721783, -0.0921586, 0.0932976, -0.356145, -1.40987, -0.498208, 0.245969, 0.00644888, -0.860894, -0.182544, -0.576631, -0.374438, 0.398192, -0.86794, 0.407687, 0.936861, 0.667791, -0.122578, 0.216126, 0.341977, 1.89751, -1.63772, -1.30438, 0.381605, 0.761368, 1.33408, -0.283563, -0.446898, -0.747549, 0.332385, 1.08601, -1.8959, -0.343493, -1.01972, 0.714939, 0.0602815, 1.53819, -0.700014, 1.19701, 0.692441, 0.205182, -1.37425, -0.742087, -0.711536, -0.500304, -0.130913, -0.767792, 0.307287, -0.679097, 2.66511, -0.443184, 0.695541, 1.69415, 1.53609, -0.635037, -0.135144, -0.450537, -0.105337, -0.844073, -0.794148, -1.26993, -0.029989, 1.26344, 0.614086, -0.410467, 0.250929, 0.911403, -0.0311083, 0.853893, -0.486828, -1.30948, 1.04435, 1.30154, -0.0989484, 0.0285836, -0.118457, 1.13655, 0.827322, -1.84675, -2.42383, -0.370874, -0.718247, -0.203151, 0.0320798, -0.283783, -0.127031, -0.00701016, 3.39039, -1.17658, -0.759631, -0.469243, -0.581371, -0.208948, 0.285828, 0.680108, -0.784031, 1.01646, 1.0188, -0.337822, -1.25625, -0.264415, 0.449394, 0.671978, -0.390034, -0.173835, -0.6541, 0.0684648, -0.424449, 1.01155, -0.813964, 1.49071, -1.23271, -0.161876, -0.395272, -1.30981, 2.26533, -1.31113, -0.229239, 0.260292, 0.970752, -0.251326, 0.304083, 0.13871};
for(i=0;i<5;i++)vectorOfVals[i]=theData[t][i];
}
/*made up data from normal dist*/
void julModShocks(int t,double * vectorOfVals)
{
int i;
static double theShocks[30][5]=
{0.795354, -2.85793, -0.422771, -3.88047, -0.150895, 0.571316, -0.893172, -0.225891, 0.236909, -1.94969, -0.535032, -0.942932, -1.3336, -1.9105, 0.73469, -1.1672, 0.422671, -0.458395, -0.395806, 0.123075, 0.268679, -0.569483, -0.917477, 0.0239555, -0.758296, -1.24116, -1.83314, 0.275472, 0.667468, 0.136316, -1.73999, -0.218102, -0.369508, 0.166112, -1.19774, 1.24544, -0.0936375, 0.885556, -0.395165, 0.0395947, 0.499787, 0.429799, -0.168071, 1.95667, -2.08641, -1.84801, -1.49703, -1.01412, 1.1298, -2.27293, 1.22116, -1.7592, 1.64448, 1.04116, -0.290173, 0.688534, 0.512071, -1.88905, 2.09408, -0.264069, -0.041646, 1.75545, 0.675526, 0.363968, -0.63124, -0.853355, 0.16116, -0.322536, 0.559482, 0.934535, -0.823338, 0.0691906, 0.777358, -0.965679, 0.900888, 0.0275707, -0.719353, 1.04527, -2.2276, -0.430615, -0.550617, -0.242605, -1.01772, 1.71713, 0.55643, 0.834238, -0.94795, -0.460843, -0.0974618, 0.518653, 1.76149, 0.52998, 1.30878, 0.511408, 1.31918, -0.331072, 0.953485, 2.24237, 0.576425, 0.623208, -1.68135, 0.277209, -0.978947, -0.347117, -0.882467, -0.74744, -1.86631, 0.854616, 0.431551, 0.217047, 1.2475, 0.221028, -0.092043, -0.387996, 0.73221, -1.26857, 0.807966, 1.13433, 1.17779, -0.99727, -1.51307, -2.19654, 2.73203, -0.590451, 1.12581, -0.0817024, 0.481503, 0.522134, -0.0928781, -1.02454, 1.23563, -0.452997, 0.232239, -0.868058, 0.357871, -1.07495, 0.98599, -0.368629, -0.494486, 0.498145, 0.670809, -0.025101, 0.639217, -0.493355, 1.47446, -0.442414, -0.408367, 1.28708, -0.198385, 1.26833};
for(i=0;i<5;i++)vectorOfVals[i]=theShocks[t][i];
}
