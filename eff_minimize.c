#include "eff_access.h"
#include "eff_util.h"
#include "eff_update.h"
#include "eff_minimize.h"
#include "eff_shanno.h"
#include "eff_global.h"
#include "eff_output.h"
#include "eff_bounds.h"
#include "eff_efield.h"
#include <math.h>
#include <stdio.h>

int numvars;
double *min_x, *min_g, *min_temp, *min_ref_x;

void AllocateMinimize(enum MinimizeType s_minmethod)
{
  numvars = NumNuclei() * 3 + NumElectrons() * 4;
  min_x = (double *) malloc(sizeof(double) * numvars);
  min_g = (double *) malloc(sizeof(double) * numvars);
  min_temp = (double *) malloc(sizeof(double) * MinStorageRequirements(s_minmethod));
  min_ref_x = (double *) malloc(sizeof(double) * numvars);
}

int MinStorageRequirements(enum MinimizeType s_minmethod)
{
  if (s_minmethod == CONJUGATE_GRADIENT)
    return 5 * numvars + 2;
  else if (s_minmethod == NEWTON)
    return numvars * (numvars + 7) / 2;
  else
    error("Minimization method not recognized.\n");
	return 0;
}

void CalcContractFG(double *x, double *f, double *g)
{
  SetMinimizePositions(x);
  ApplyPositionBounds();
  UpdateEnergyForces();
  ApplyExternalField(0);
  *f = GetTotalPE();
  GetMinimizeForces(g);
  if (fabs(*f + 510.773833) < 0.0001)
  {
    printf("%f\n", *f);
    double t;
    int i, n;
    for (n = 0; n < 100; n++)
    {
      for (i = 0; i < numvars; i++)
      {
        x[i] -= 0.0001 * g[i];
      } 
      SetMinimizePositions(x);
      //ApplyPositionBounds();
      UpdateEnergyForces();
      ApplyExternalField(0);
      printf("%f\n", GetTotalPE());
    }
    exit(1);
  }
}

enum MinimizeResult Minimize(enum MinimizeType s_minmethod, double eps, double acc, int max_numsteps, int print_every, void (*MinCallback)(int, int, double, double))
{
  // Interface with Shanno's CG/BFGS minimizer
  double min_f;
  int num_fevals, num_iterations;
  int result;
  
  GetInitialPosition(min_x);

  // Callback function called with 
  //   # iters, # func evals, f, grad2

  result = conmin(numvars, min_x, &min_f, min_g, &num_fevals, &num_iterations,
         eps, max_numsteps, min_temp, print_every, MinStorageRequirements(s_minmethod), acc, s_minmethod,
         &CalcContractFG, MinCallback);

  if (result == 0)
			 return NORMAL;
  else if (result == 1)
			 return FAIL_EXCEED;
  else if (result == 2)
			 return FAIL_LINSEARCH;
  else if (result == 3)
			 return FAIL_NODESCENT;
  return NORMAL;

}

void GetMinimizeForces(double *forces)
{
  int i, idx;
	double fx, fy, fz, fr;
  double x, y, z, r;

  idx = 0;
  for (i = 0; i < NumNuclei(); i++)
	{
    GetNuclearForce(i, &fx, &fy, &fz);
  //fx = fy = fz = 0;
    forces[idx] = -fx; idx++;  /* x */
    forces[idx] = -fy; idx++;  /* y */
    forces[idx] = -fz; idx++;  /* z */
	}
	for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronForce(i, &fx, &fy, &fz, &fr);
  //fx = fy = fz = fr = 0;
    GetElectronPosition(i, &x, &y, &z, &r);
    forces[idx] = -fx; idx++;   /* x */
    forces[idx] = -fy; idx++;   /* y */
    forces[idx] = -fz; idx++;   /* z */
    forces[idx] = -fr * r; idx++;   /* r */
    /* transform radius variable: r = exp(t) => df/dt = df/dr dr/dt = df/dr * r */
  }
}

void GetMinimizeReferencePositions(double *positions)
{
	int i, idx;
	double x, y, z, r;

  idx = 0;
  for (i = 0; i < NumNuclei(); i++)
	{
    GetNuclearPosition(i, &x, &y, &z);
    positions[idx] = x; idx++;  /* x */
    positions[idx] = y; idx++;  /* y */
    positions[idx] = z; idx++;  /* z */
	}
	for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronPosition(i, &x, &y, &z, &r);
    positions[idx] = x; idx++;   /* x */
    positions[idx] = y; idx++;   /* y */
    positions[idx] = z; idx++;   /* z */
    positions[idx] = r; idx++;   /* r */
  }
}

void GetInitialPosition(double *positions)
{
  int i, idx;
	double x, y, z, r;

  idx = 0;
  for (i = 0; i < NumNuclei(); i++)
	{
    GetNuclearPosition(i, &x, &y, &z);
    positions[idx] = x; idx++;  /* x */
    positions[idx] = y; idx++;  /* y */
    positions[idx] = z; idx++;  /* z */
	}
	for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronPosition(i, &x, &y, &z, &r);
    positions[idx] = x; idx++;   /* x */
    positions[idx] = y; idx++;   /* y */
    positions[idx] = z; idx++;   /* z */
    positions[idx] = log(r); idx++;   /* r */
  }
}

void SetMinimizePositions(double *positions)
{
	int i, idx;
  idx = 0;
 	for (i = 0; i < NumNuclei(); i++)
  {
		SetNuclearPosition(i, positions[idx], positions[idx+1], positions[idx+2]);
    idx += 3;
  }
	for (i = 0; i < NumElectrons(); i++)
  {
    /* transform radius variable: r = exp(t) => df/dt = df/dr dr/dt = df/dr * r */
		SetElectronPosition(i, positions[idx], positions[idx+1], positions[idx+2], exp(positions[idx+3]));
    idx += 4;
  }
}
