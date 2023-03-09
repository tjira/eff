#include <time.h>
#include "eff_timing.h"

double time_elapsed;
clock_t t1, t2;

void StartTimer()
{
  t1 = clock();
}

void StopTimer()
{
  t2 = clock();
}

double TimeElapsed()
{
  return (double) (t2 - t1) / CLOCKS_PER_SEC;
}



