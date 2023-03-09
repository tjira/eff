#include <math.h>
#include <stdio.h>
#include "eff_access.h"
#include "eff_util.h"
#include "eff_update.h"
#include "eff_dynamics.h"
#include "eff_bounds.h"
#include "eff_constants.h"
#include "eff_pressure.h"
#include "eff_efield.h"
#include "eff_collapse.h"

double dt = 0.005;
FILE *logfile;

double temp = 0;
int g_useTemp = 0;
int numvars;

double *dyn_x, *dyn_v, *dyn_f, *dyn_prev_f, *dyn_mass;
double dynamics_time;

int GetNumDOF()
{
  return numvars;
}

void AllocateDynamics()
{
  numvars = NumNuclei() * 3 + NumElectrons() * 4;
  dyn_x      = (double *) malloc(sizeof(double) * numvars);
  dyn_v      = (double *) malloc(sizeof(double) * numvars);
  dyn_f      = (double *) malloc(sizeof(double) * numvars);
  dyn_prev_f = (double *) malloc(sizeof(double) * numvars);
  dyn_mass   = (double *) malloc(sizeof(double) * numvars);
	if (dyn_x == 0 || dyn_f == 0 || dyn_mass == 0) error("Could not allocate space for dynamics objects.\n");
}

void InitializeDynamics()
{
	int i, idx;
  double mass;
  
  /* initialize positions */
	GetDynamicsPositions();

  /* initialize velocities */
  for (i = 0; i < numvars; i++)
    dyn_v[i] = 0;

  /* initialize masses */
  idx = 0;
  for (i = 0; i < NumNuclei(); i++)
	{
    mass = GetMass((int) GetNuclearCharge(i));
    dyn_mass[idx] = mass; idx++;  /* x */
    dyn_mass[idx] = mass; idx++;  /* y */
    dyn_mass[idx] = mass; idx++;  /* z */
	}
	for (i = 0; i < NumElectrons(); i++)
  {
    dyn_mass[idx] = params.electron_mass; idx++;   /* x */
    dyn_mass[idx] = params.electron_mass; idx++;   /* y */
    dyn_mass[idx] = params.electron_mass; idx++;   /* z */
    dyn_mass[idx] = 0.75 * params.electron_mass; idx++;   /* r */
  }

  // initialize forces
	UpdateEnergyForces();
	GetDynamicsForces();

  // initialize time
  dynamics_time = 0;
}

void InitializeRandomVelocities(double temp)
{
  int i;
  double m, sigma, randn;

  /* Set all velocities to zero */
  for (i = 0; i < numvars; i++)
    dyn_v[i] = 0;

  /* Initialize nuclear velocities */
  /*
  for (i = 0; i < NumNuclei() * 3; i++)
  {
    sigma = sqrt(KB * temp / dyn_mass[i]);
    dyn_v[i] = rand_gauss() * sigma;  
  }
  */

  /* Initialize electron velocities */
  for (i = NumNuclei() * 3; i < NumNuclei() * 3 + NumElectrons() * 4; i++)
  {
    sigma = sqrt(KB * temp / dyn_mass[i]);
    dyn_v[i] = rand_gauss() * sigma;  
  }

  /* Initialize all velocities */
  /*
  for (i = 0; i < NumNuclei() * 3 + NumElectrons() * 4; i++)
  {
    sigma = sqrt(KB * temp / dyn_mass[i]);
    dyn_v[i] = rand_gauss() * sigma;  
  }
  */
}

void Dynamics(enum ThermostatType s_thermostat, double temperature)
{
  // Try to collpase electrons that are expanding too fast.
  int i;
  if (params.collapse_move > 0)
  {
    for (i = 0; i < NumElectrons(); i++)
      CollapseElectron(i);
    GetDynamicsPositions();  // real --> dyn
  }

	UpdateDynamicsPositions();
  
  SetDynamicsPositions();  // dyn --> real
  ApplyPositionBounds();
  UpdateEnergyForces();
  GetDynamicsPositions();  // real --> dyn

  ApplyExternalField(dynamics_time);
  UpdateDynamicsVirial();

	GetDynamicsForces();
	UpdateDynamicsVelocities(s_thermostat, temperature);
  dynamics_time += dt;
}

void UpdateDynamicsPositions()
{
	int i;
	double m;

  for (i = 0; i < numvars; i++)
  {
    m = 0.5 / dyn_mass[i];
    dyn_x[i] += dt * (dyn_v[i] + m * dt * dyn_f[i]);
    dyn_prev_f[i] = dyn_f[i];
  }
}

void UpdateDynamicsVelocities(enum ThermostatType s_thermostat, double temperature)
{
	int i;
	double m;
	double sigma, randn;

  for (i = 0; i < numvars; i++)
  {
    m = 0.5 / dyn_mass[i];
    if (s_thermostat == NOSE_HOOVER)
    {
      sigma = sqrt(KB * temperature / dyn_mass[i]);
      randn = rand_uni();
      if (randn < 0.1 * dt)  /* larger parameter --> more contact with heat bath */
        dyn_v[i] = rand_gauss() * sigma;  /* random temperature bath */
      else
        dyn_v[i] += m * dt * (dyn_f[i] + dyn_prev_f[i]); /* regular forces */
    }
    else
    {
      dyn_v[i] += m * dt * (dyn_f[i] + dyn_prev_f[i]); /* regular forces */
    }
	}
}

void GetDynamicsForces()
{
	int i, idx;
	double fx, fy, fz, fr;

  idx = 0;
  for (i = 0; i < NumNuclei(); i++)
	{
    GetNuclearForce(i, &fx, &fy, &fz);
    dyn_f[idx] = fx; idx++;  /* x */
    dyn_f[idx] = fy; idx++;  /* y */
    dyn_f[idx] = fz; idx++;  /* z */
	}
	for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronForce(i, &fx, &fy, &fz, &fr);
    dyn_f[idx] = fx; idx++;   /* x */
    dyn_f[idx] = fy; idx++;   /* y */
    dyn_f[idx] = fz; idx++;   /* z */
    dyn_f[idx] = fr; idx++;   /* r */
  }
}

void GetDynamicsPositions()
{
	int i, idx;
	double x, y, z, r;

  idx = 0;
  for (i = 0; i < NumNuclei(); i++)
	{
    GetNuclearPosition(i, &x, &y, &z);
    dyn_x[idx] = x; idx++;  /* x */
    dyn_x[idx] = y; idx++;  /* y */
    dyn_x[idx] = z; idx++;  /* z */
	}
	for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronPosition(i, &x, &y, &z, &r);
    dyn_x[idx] = x; idx++;   /* x */
    dyn_x[idx] = y; idx++;   /* y */
    dyn_x[idx] = z; idx++;   /* z */
    dyn_x[idx] = r; idx++;   /* r */
  }
}

void SetDynamicsPositions()
{
	int i, idx;
  idx = 0;
 	for (i = 0; i < NumNuclei(); i++)
  {
		SetNuclearPosition(i, dyn_x[idx], dyn_x[idx+1], dyn_x[idx+2]);
    idx += 3;
  }
	for (i = 0; i < NumElectrons(); i++)
  {
		SetElectronPosition(i, dyn_x[idx], dyn_x[idx+1], dyn_x[idx+2], dyn_x[idx+3]);
    idx += 4;
  }
}

/* dynamics accessor functions */
void SetTimeStep(double s_dt)
{
  /* input in fs, but our program works internally in modified atomic time units
     such that one atu = 1.0327499 fs 
  */
	dt = s_dt / T0_IN_FS;
}

double GetTemperature()
{
  int i;
  double meas_temp = 0;
  for (i = 0; i < numvars; i++)
    meas_temp += dyn_mass[i] * dyn_v[i] * dyn_v[i]; 

  /* T = sum((1/2)mv^2 / (1/2)k N) -- depracated */
  /* Assume 3/2 k T per nucleus */
  return meas_temp / (KB * 3 * NumNuclei());
}

// Accessor functions for velocities

void SetNuclearVelocity(int i, double vx, double vy, double vz)
{
  dyn_v[3*i  ] = vx;
  dyn_v[3*i+1] = vy;
  dyn_v[3*i+2] = vz;
}

void SetElectronVelocity(int i, double vx, double vy, double vz, double vr)
{
  int i0 = NumNuclei() * 3;
  dyn_v[4*i  +i0] = vx;
  dyn_v[4*i+1+i0] = vy;
  dyn_v[4*i+2+i0] = vz;
  dyn_v[4*i+3+i0] = vr;
}

void GetNuclearVelocity(int i, double *vx, double *vy, double *vz)
{
  *vx = dyn_v[3*i  ];
  *vy = dyn_v[3*i+1];
  *vz = dyn_v[3*i+2];
}

void GetElectronVelocity(int i, double *vx, double *vy, double *vz, double *vr)
{
  int i0 = NumNuclei() * 3;
  *vx = dyn_v[4*i  +i0];
  *vy = dyn_v[4*i+1+i0];
  *vz = dyn_v[4*i+2+i0];
  *vr = dyn_v[4*i+3+i0];
}

// Accessor functions for masses

void SetNuclearMass(int i, double mx, double my, double mz)
{
  dyn_mass[3*i  ] = mx;
  dyn_mass[3*i+1] = my;
  dyn_mass[3*i+2] = mz;
}

void SetElectronMass(int i, double mx, double my, double mz, double mr)
{
  int i0 = NumNuclei() * 3;
  dyn_mass[4*i  +i0] = mx;
  dyn_mass[4*i+1+i0] = my;
  dyn_mass[4*i+2+i0] = mz;
  dyn_mass[4*i+3+i0] = mr;
}

void GetNuclearMass(int i, double *mx, double *my, double *mz)
{
  *mx = dyn_mass[3*i  ];
  *my = dyn_mass[3*i+1];
  *mz = dyn_mass[3*i+2];
}

void GetElectronMass(int i, double *mx, double *my, double *mz, double *mr)
{
  int i0 = NumNuclei() * 3;
  *mx = dyn_mass[4*i  +i0];
  *my = dyn_mass[4*i+1+i0];
  *mz = dyn_mass[4*i+2+i0];
  *mr = dyn_mass[4*i+3+i0];
}

double GetTotalKE()
{
  int i;
  double total_ke = 0;
  for (i = 0; i < NumNuclei(); i++)
    total_ke += GetNucleusKE(i);
  for (i = 0; i < NumElectrons(); i++)
    total_ke += GetElectronKE(i);
  
  return total_ke;
}

double GetNucleusKE(int i)
{
  double ke = 0;
  ke += dyn_mass[3*i  ] * dyn_v[3*i  ] * dyn_v[3*i  ];
  ke += dyn_mass[3*i+1] * dyn_v[3*i+1] * dyn_v[3*i+1];
  ke += dyn_mass[3*i+2] * dyn_v[3*i+2] * dyn_v[3*i+2];
  return 0.5 * ke;
}

double GetElectronKE(int i)
{
  int i0 = NumNuclei() * 3;
  double ke = 0;
  ke += dyn_mass[4*i  +i0] * dyn_v[4*i  +i0] * dyn_v[4*i  +i0];
  ke += dyn_mass[4*i+1+i0] * dyn_v[4*i+1+i0] * dyn_v[4*i+1+i0];
  ke += dyn_mass[4*i+2+i0] * dyn_v[4*i+2+i0] * dyn_v[4*i+2+i0];
  ke += dyn_mass[4*i+3+i0] * dyn_v[4*i+3+i0] * dyn_v[4*i+3+i0];
  return 0.5 * ke;
}

double GetElectronTranslationKE(int i)
{
  /* Gets kinetic energy of only the translational DOF of electron */
  int i0 = NumNuclei() * 3;
  double ke = 0;
  ke += dyn_v[4*i  +i0] * dyn_v[4*i  +i0];
  ke += dyn_v[4*i+1+i0] * dyn_v[4*i+1+i0];
  ke += dyn_v[4*i+2+i0] * dyn_v[4*i+2+i0];
  return 0.5 * dyn_mass[4*i+i0] * ke;
}

void UpdateDynamicsVirial()
{
  int i;
  for (i = 0; i < NumNuclei(); i++)
    AddKineticEnergyVirial(GetNucleusKE(i), BOTH);

  for (i = 0; i < NumElectrons(); i++)
  {
    AddKineticEnergyVirial(GetElectronTranslationKE(i), RIGID);
    AddKineticEnergyVirial(GetElectronKE(i), FLEXIBLE);
  }
}
