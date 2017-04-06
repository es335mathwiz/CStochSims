
/*dvalsInfo*/
void rbcExampleShocks(int t,double * vectorOfVals)
{
int i;
#include "rbcTryCShocksForInclude.h"/*dstr;*/
for(i=0;i<4;i++)vectorOfVals[i]=0;
for(i=0;i<4-0;i++)vectorOfVals[i]=theShocks[t][i];
}


