#include "eff_global.h"
#include "eff_driver.h"
#include "eff_update.h"
#include "eff_create.h"
#include "eff_access.h"
#include "eff_dynamics.h"
#include "eff_util.h"
#include "eff_output.h"
#include "eff_input.h"
#include "eff_initial.h"
#include "eff_ewald.h"
#include "eff_timing.h"
#include "eff_pressure.h"
#include "eff_bounds.h"
#include "eff_efield.h"
#include "eff_properties.h"
#include "eff_constants.h"
#include "eff_cores.h"
#include "eff_eigen.h"
#include "eff_minimize.h"
#include "eff_params.h"

#include <stdio.h>
#include <math.h>

/* Output and restart file pointers */
FILE *out_fp;
char restart_filname[2048];

int main(int argc, char **argv)
{
	int i;
	double energy, fx, fy, fz, fr;
	double x, y, z;

  if (argc != 2)
  {
    printf("eff_p 1-28-06\n");
    printf("Usage: eff <input.cfg>\n");
    printf("Output files produced: .out, .cfg.restart\n\n");
    return 1;
  }


  /* Allocate enough space for all the nuclei, electrons, initial velocities, and initial masses */
  /* Peek into the input file first */
  FILE *fp = fopen(argv[1], "r");
  int a_numnuclei = CountLines(fp, "@nuclei");
  int a_numelectrons = CountLines(fp, "@electrons");
  int a_numnucv = CountLines(fp, "@nuc_velocities");
  int a_numelecv = CountLines(fp, "@elec_velocities");
  int a_numnucm = CountLines(fp, "@nuc_masses");
  int a_numelecm = CountLines(fp, "@elec_masses");

  Initialize(a_numnuclei, a_numelectrons);
  AllocateInitialVelocities(a_numnucv, a_numelecv);
  AllocateInitialMasses(a_numnucm, a_numelecm);

  /* Set default parameters and load in input file */
  InitializeParameters();
  SetDefaultParameters();
  ProcessInputFile(fp);
  fclose(fp);

  /* Set up output files */
  char basename[2048], extension[255];
  GetExtension(argv[1], extension);
  GetBasename(argv[1], basename);
  
  char out_filname[2048];
  sprintf(out_filname, "%s.eff", basename);
  sprintf(restart_filname, "%s.restart", argv[1]);
  out_fp = fopen(out_filname, "w");

  /* Output options and fixed header */
  OutputParams("[params]", out_fp);
  OutputFixedHeader(out_fp);  

  /* Set up energy updater */
  IntitializeEnergyUpdater(params.s_cutoff, params.periodic);

  /* Set up core electrons */
  InitializeCores();

  /* Set up external electric field */
  SetupExternalField(params.Ex, params.Ey, params.Ez, params.Efreq, params.Epacket_duration);

  /* Set up ewald solver if periodic boundary conditions are present */
  if (params.periodic)
  {
    AllocateEwald();
    if (params.ewald_autoset)
    {
      params.ewald_r_cutoff = EwaldRCutoff(params.ewald_log_precision, 2.0 / (params.ewald_re_cutoff * params.ewald_re_cutoff), 2.0 / (params.ewald_max_re * params.ewald_max_re));
      params.ewald_k_cutoff = EwaldKCutoff(params.ewald_log_precision, 2.0 / (params.ewald_re_cutoff * params.ewald_re_cutoff));
    }
    InitializeEwald(params.x_bound[1] - params.x_bound[0], 
                    params.y_bound[1] - params.y_bound[0],
                    params.z_bound[1] - params.z_bound[0],
                    2.0 / (params.ewald_re_cutoff * params.ewald_re_cutoff), params.ewald_r_cutoff, params.ewald_k_cutoff, 2.0 / (params.ewald_nuc_r * params.ewald_nuc_r));
  }

// Run calculation 

  StartTimer();
  if (params.calc == SINGLE_PT)
    RunSinglePoint(out_fp);
  else if (params.calc == MINIMIZE)
    RunMinimize(out_fp);
  else if (params.calc == DYNAMICS)
    RunDynamics(out_fp);
  StopTimer();
  OutputTimeElapsed(out_fp);
  
  //PrintAnalyticForces();
  //PrintNumericalForces();

  /* Output restart file */
  OutputRestartFile(restart_filname);

  fclose(out_fp);
}

void RunSinglePoint()
{
  UpdateEnergyForces();
  ApplyExternalField(0);  // external electric field

  //printf("rigid pressure = %f\n", GetRigidPEPressure(GetVolume()));
  //printf("flexible pressure = %f\n", GetFlexiblePEPressure(GetVolume()));
  //printf("energy = %f\n", GetTotalPE());

double r;
for (r = 0; r >= -11.5; r -= 0.1)
{
// hollow
// ShiftMonomer(r, 0.007252167, -0.008538281, -0.99993725);

// atop
ShiftMonomer(r, 0.006484466, -0.614436608, -0.788939545);

// bridge
//ShiftMonomer(r, 0.31742366, -0.18847568, -0.92936491);

//  UpdateEnergyForces();
//    OutputPositions(out_fp);
  RunMinimize(out_fp);
  printf("%f\t%f\n", r, GetTotalPE());

// hollow
// ShiftMonomer(-r, 0.007252167, -0.008538281, -0.99993725);

// atop
ShiftMonomer(-r, 0.006484466, -0.614436608, -0.788939545);

// bridge
//ShiftMonomer(-r, 0.31742366, -0.18847568, -0.92936491);
}

/*
double r;
for (r = 0.5; r < 2.0; r += 0.1)
{
  ScaleElectrons(r);
  UpdateEnergyForces();
  printf("%f\t%f\n", r, GetTotalPE());
  ScaleElectrons(1.0/r);
}
*/

  if (params.output_energy_forces == ALL || params.output_energy_forces == END) 
    OutputEnergyForces(out_fp);
  if (params.output_position == ALL || params.output_position == END)
    OutputPositions(out_fp);
  OutputEnergy(out_fp);
  if (params.output_restart == ALL || params.output_restart == END)
    OutputRestartFile(restart_filname);
}

void MinimizeCallback(int num_iterations, int num_fevals, double f_val, double g_sq)
{
  fprintf(out_fp, "[Minimize_iter]\t%i\t%i\t%e\t%e\n", num_iterations, num_fevals, f_val, g_sq);
  if (params.output_energy_forces == ALL) OutputEnergyForces(out_fp);
  if (params.output_position == ALL)      OutputPositions(out_fp);
  if (params.output_restart == ALL)       OutputRestartFile(restart_filname);
}

void RunMinimize()
{
  AllocateMinimize(params.min);
  if (params.min == CONJUGATE_GRADIENT)
    fprintf(out_fp, "[Minimize_type]\tConjugate gradient.\n");
  else if (params.min == NEWTON)
    fprintf(out_fp, "[Minimize_type]\tNewton.\n");

  /* Minimization parameters */  
  double eps = .00001;
  double acc = 1e-20;

  /* Minimize */
  fprintf(out_fp, "[Minimize_iter_header]\t# iters\t# func evals\tf\tgrad2\n");
  enum MinimizeResult min_result = Minimize(params.min, eps, acc, params.num_steps, params.print_every, &MinimizeCallback);

  /* Results of minimization */
  if (min_result == NORMAL)
    fprintf(out_fp, "[Minimize_result] Converged.\n");
  else if (min_result == FAIL_EXCEED)
    fprintf(out_fp, "[Minimize_result] Failed, maximum number of function evaluations used.\n");
  else if (min_result == FAIL_LINSEARCH)
    fprintf(out_fp, "[Minimize_result] Failed, linear search fails to improve function value, is gradient correct?\n");
  else if (min_result == FAIL_NODESCENT)
    fprintf(out_fp, "[Minimize_result] Failed, search vector not a descent direction, convergence critereon too strict?\n");

  /* End energy */
  if (params.output_energy_forces == ALL || params.output_energy_forces == END)
    OutputEnergyForces(out_fp);
  if (params.output_position == ALL || params.output_position == END)
    OutputPositions(out_fp);
  if (params.output_restart == ALL || params.output_restart == END)
    OutputRestartFile(restart_filname);
  OutputEnergy(out_fp);

/*
  double r;
  for (r = -3.0; r <= 3.1; r += 0.1)
    printf("%e\n", ElectronDensity(r/sqrt(3), r/sqrt(3), -r/sqrt(3)));
    //printf("%e\n", ElectronDensity(r, 0, 0));
  return;
*/
}

void RunDynamics()
{
  AllocateDynamics();
  InitializeDynamics();
  rand_seed(params.rand_seed);

  /* Initialize masses with those given in starting file */
  double mx, my, mz, mr;
  int i, idx;
  for (i = 0; i < GetNumNucM(); i++)
  {
    GetNucM(i, &idx, &mx, &my, &mz);
    SetNuclearMass(idx, mx, my, mz);
  }
  for (i = 0; i < GetNumElecM(); i++)
  {
    GetElecM(i, &idx, &mx, &my, &mz, &mr);
    SetElectronMass(idx, mx, my, mz, mr);
  }

  InitializeRandomVelocities(params.start_temperature);
  OutputDynamicsMasses(out_fp);

  /* Initialize velocites with those given in starting file */
  double vx, vy, vz, vr;
  for (i = 0; i < GetNumNucV(); i++)
  {
    GetNucV(i, &idx, &vx, &vy, &vz);
    SetNuclearVelocity(idx, vx, vy, vz);
  }
  for (i = 0; i < GetNumElecV(); i++)
  {
    GetElecV(i, &idx, &vx, &vy, &vz, &vr);
    SetElectronVelocity(idx, vx, vy, vz, vr);
  }

  /* Dynamics parameters */
  SetTimeStep(params.dt);

  /* Dynamics */
  fprintf(out_fp, "[dynamics_header]\t(fs)\tkinetic_e\tpotential_e\ttotal_e\ttarget_temp\tmeasured_temp\tp_ke_rigid\tp_ke_flexible\tp_pe_rigid\tp_pe_flexible\n");
	for (i = 0; i < params.num_steps; i++)
	{
		if (i % params.print_every == 0) 
    {
      double ke, pe;
      ke = GetTotalKE();
      pe = GetTotalPE();
      fprintf(out_fp, "[dynamics_iter]\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%e\t%e\t%e\n", i * params.dt, 
        ke, pe, ke + pe, params.start_temperature, GetTemperature(), 
        GetRigidKEPressure(GetVolume()), 
        GetFlexibleKEPressure(GetVolume()),
        GetRigidPEPressure(GetVolume()), 
        GetFlexiblePEPressure(GetVolume())
        );
    }
		Dynamics(params.thermostat, params.start_temperature);
		if (i % params.print_every == 0) 
    {
      if (params.output_energy_forces == ALL) OutputEnergyForces(out_fp);
      if (params.output_position == ALL)      OutputPositions(out_fp);
      if (params.output_velocity == ALL)      OutputVelocities(out_fp);
      if (params.output_restart == ALL)       OutputRestartFile(restart_filname);
    }
	}

  /* End energy */
  if (params.output_energy_forces == ALL || params.output_energy_forces == END)
    OutputEnergyForces(out_fp);
  if (params.output_position == ALL || params.output_position == END)
    OutputPositions(out_fp);
  if (params.output_velocity == ALL || params.output_velocity == END)
    OutputVelocities(out_fp);
  if (params.output_restart == ALL || params.output_restart == END)
    OutputRestartFile(restart_filname);
  OutputEnergy(out_fp);
}

