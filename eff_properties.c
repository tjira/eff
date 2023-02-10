#include "eff_properties.h"
#include "eff_update.h"
#include "eff_access.h"
#include <math.h>

void CalcDipole(double *dipole_x, double *dipole_y, double *dipole_z)
{
  int i;
  double x, y, z, r, q;
  double ux, uy, uz;
  ux = uy = uz = 0;

  for (i = 0; i < NumNuclei(); i++)
  {
    q = GetNuclearCharge(i);
    GetNuclearPosition(i, &x, &y, &z);
    ux += q * x;
    uy += q * y;
    uz += q * z;
  }

  for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronPosition(i, &x, &y, &z, &r);
    ux += -x;
    uy += -y;
    uz += -z;
  }

  *dipole_x = ux;
  *dipole_y = uy;
  *dipole_z = uz;
}

double ElectronDensity(double x, double y, double z)
{
  double cx, cy, cz, re;
  double dx, dy, dz;

  double density = 0;
  int i;
  for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronPosition(i, &cx, &cy, &cz, &re);

    dx = x - cx; dy = y - cy; dz = z - cz;
    double s = sqrt(dx * dx + dy * dy + dz * dz) / re;
if (i == 2 || i == 3)
    density += pow(2.0 / (3.1415926535 * re * re), 1.5) * exp(-2 * s * s);
//    density += (0.743401 / (re * re * re)) * exp(-4.44605 * s * s / (1 + 0.745491 * s));
  }
  return density;
}

