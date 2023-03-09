#include "eff_bounds.h"
#include "eff_global.h"
#include "eff_access.h"

double GetVolume()
{
  double lx, ly, lz;
  lx = (params.x_bound[1] - params.x_bound[0]);
  ly = (params.y_bound[1] - params.y_bound[0]);
  lz = (params.z_bound[1] - params.z_bound[0]);
  return lx * ly * lz;
}

/* Change atomic positions at boundaries */

double BoundX(double x)
{
  if (!params.periodic) return x;
  if (x < params.x_bound[0])
      return x + (params.x_bound[1] - params.x_bound[0]);
  if (x > params.x_bound[1])
      return x - (params.x_bound[1] - params.x_bound[0]);
  return x;
}

double BoundY(double y)
{
  if (!params.periodic) return y;
  if (y < params.y_bound[0])
      return y + (params.y_bound[1] - params.y_bound[0]);
  if (y > params.y_bound[1])
      return y - (params.y_bound[1] - params.y_bound[0]);
  return y;
}

double BoundZ(double z)
{
  if (!params.periodic) return z;
  if (z < params.z_bound[0])
      return z + (params.z_bound[1] - params.z_bound[0]);
  if (z > params.z_bound[1])
      return z - (params.z_bound[1] - params.z_bound[0]);
  return z;
}

/* Change interatomic distances at boundaries */

double BoundDX(double dx)
{
  if (!params.periodic) return dx;
  double length = params.x_bound[1] - params.x_bound[0];
  if (dx < -0.5 * length) return dx + length;
  if (dx > 0.5 * length) return dx - length;
  return dx;
}

double BoundDY(double dy)
{
  if (!params.periodic) return dy;
  double length = params.y_bound[1] - params.y_bound[0];
  if (dy < -0.5 * length) return dy + length;
  if (dy > 0.5 * length) return dy - length;
  return dy;
}

double BoundDZ(double dz)
{
  if (!params.periodic) return dz;
  double length = params.z_bound[1] - params.z_bound[0];
  if (dz < -0.5 * length) return dz + length;
  if (dz > 0.5 * length) return dz - length;
  return dz;
}

/* Applying conditions all at once */

void ApplyPositionBounds()
{
  int i;
  double x, y, z, r;
  for (i = 0; i < NumNuclei(); i++)
  {
    GetNuclearPosition(i, &x, &y, &z);
    x = BoundX(x); y = BoundY(y); z = BoundZ(z);
    SetNuclearPosition(i, x, y, z);
  }
  for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronPosition(i, &x, &y, &z, &r);
    x = BoundX(x); y = BoundY(y); z = BoundZ(z);
    SetElectronPosition(i, x, y, z, r);
  }
}
