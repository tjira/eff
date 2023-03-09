#include "eff_access.h"
#include "eff_cores.h"
#include <stdlib.h>
#include <math.h>

int *array_isCore;

int TestCore(int i)
{
  // If the electron is close to a heavy nucleus and small enough,
  // mark it as a core electron.

  int j;
  double ex, ey, ez, er;
  GetElectronPosition(i, &ex, &ey, &ez, &er);

  int numnuc = NumNuclei();
  for (j = 0; j < numnuc; j++)
  {
    double Z = GetNuclearCharge(j);
    if (Z > 1)
    {
      double r_core = calc_rcore(Z);
      double factor1 = 1.5;
      double factor2 = 0.5;

      // First condition: electron needs to be smaller than factor1 * rcore for that nucleus
      if (er < factor1 * r_core)
      {
        double nx, ny, nz, dx, dy, dz;
        GetNuclearPosition(j, &nx, &ny, &nz);
        dx = nx - ex; dy = ny - ey; dz = nz - ez;

        // Second condition: electron needs to be within factor2 * rcore for that nucleus
        if (dx * dx + dy * dy + dz * dz < factor2 * factor2 * r_core * r_core)
          return 1;
      }
    }
  }
  return 0;
}

double calc_rcore(double Z)
{
  return 3 * sqrt(3.14159265359) / (2 * sqrt(2.0) * Z - 1.0);
}

void InitializeCores()
{
  int numelecs = NumElectrons();

  array_isCore = (int *) malloc(sizeof(int) * numelecs);

  int i;
  int numCores = 0;
  for (i = 0; i < numelecs; i++)
  {
    array_isCore[i] = TestCore(i);
    numCores += array_isCore[i];
  }

//  if (numCores != NumHeavyNuclei() * 2)
//    error("Number of core electrons (%i) doesn't match number of nuclei (%i).\n", numCores, NumHeavyNuclei());
}

int isCore(int i)
{
  return array_isCore[i];
}

