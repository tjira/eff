#include "eff_global.h"
#include "eff_create.h"
#include "eff_dynamics.h"
#include "eff_output.h"
#include "eff_access.h"
#include "eff_util.h"
#include "eff_timing.h"
#include "eff_constants.h"
#include "eff_efield.h"
#include "eff_update.h"
#include "eff_initial.h"
#include "eff_properties.h"
#include <stdarg.h>
#include <stdio.h>
#include <math.h>

void OutputEnergyForces(FILE *out_fp)
{
  fprintf(out_fp, "[energy_force_nuc_header]\t#\tE\tfx\tfy\tfz\n");
  int i;
  double fx, fy, fz, fr, E;
  for (i = 0; i < NumNuclei(); i++)
  {
    E = GetNucleusPE(i);
    GetNuclearForce(i, &fx, &fy, &fz);
    fprintf(out_fp, "[energy_force_nuc]\t%i\t%f\t%f\t%f\t%f\n", i+1, E, fx, fy, fz);
  }
  fprintf(out_fp, "[energy_force_elec_header]\t#\tE\tfx\tfy\tfz\tfr\n");
  for (i = 0; i < NumElectrons(); i++)
  {
    E = GetElectronPE(i);
    GetElectronForce(i, &fx, &fy, &fz, &fr);
    fprintf(out_fp, "[energy_force_elec]\t%i\t%f\t%f\t%f\t%f\t%f\n", i+1, E, fx, fy, fz, fr);
  }
}

void OutputEnergy(FILE *out_fp)
{
  double ux, uy, uz;
  CalcDipole(&ux, &uy, &uz);
  fprintf(out_fp, "[total_dipole]\t%f\n", DEBYE_PER_AU * sqrt(ux * ux + uy * uy + uz * uz));

  fprintf(out_fp, "[total_energy]\t%f\n", GetTotalPE());
}

void OutputTimeElapsed(FILE *out_fp)
{
  fprintf(out_fp, "[time_elapsed]\t%f\n", TimeElapsed());
}

void OutputPositions(FILE *out_fp)
{
  fprintf(out_fp, "[position_nuc_header]\t#\tx\ty\tz\n");
  int i;
  double x, y, z, r;
  for (i = 0; i < NumNuclei(); i++)
  {
    GetNuclearPosition(i, &x, &y, &z);
    fprintf(out_fp, "[position_nuc]\t%i\t%f\t%f\t%f\n", i+1, x, y, z);
  }
  fprintf(out_fp, "[position_elec_header]\t#\tE\tx\ty\tz\tr\n");
  for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronPosition(i, &x, &y, &z, &r);
    fprintf(out_fp, "[position_elec]\t%i\t%f\t%f\t%f\t%f\n", i+1, x, y, z, r);
  }
}

void OutputVelocities(FILE *out_fp)
{
  fprintf(out_fp, "[velocity_nuc_header]\t#\tvx\tvy\tvz\n");
  int i;
  double vx, vy, vz, vr;
  for (i = 0; i < NumNuclei(); i++)
  {
    GetNuclearVelocity(i, &vx, &vy, &vz);
    fprintf(out_fp, "[velocity_nuc]\t%i\t%f\t%f\t%f\n", i+1, vx, vy, vz);
  }
  fprintf(out_fp, "[velocity_elec_header]\t#\tE\tvx\tvy\tvz\tvr\n");
  for (i = 0; i < NumElectrons(); i++)
  {
    GetElectronVelocity(i, &vx, &vy, &vz, &vr);
    fprintf(out_fp, "[velocity_elec]\t%i\t%f\t%f\t%f\t%f\n", i+1, vx, vy, vz, vr);
  }
}

void OutputRestartFile(char *restart_filname)
{
  int i, spin;
  double x, y, z, r, q;

  FILE *restart_fp = fopen(restart_filname, "w");
  fprintf(restart_fp, "@params\n");
  OutputParams("", restart_fp);

  fprintf(restart_fp, "@nuclei\n");
  for (i = 0; i < NumNuclei(); i++)
  {
    q = GetNuclearCharge(i);
    GetNuclearPosition(i, &x, &y, &z);
    fprintf(restart_fp, "%f %f %f %f\n", x, y, z, q);
  }

  fprintf(restart_fp, "@electrons\n");
  for (i = 0; i < NumElectrons(); i++)
  {
    spin = GetElectronSpin(i);
    GetElectronPosition(i, &x, &y, &z, &r);
    //fprintf(fp, "%i 1 1 1\n", spin);
    //fprintf(fp, "1.0 000 %f %f %f %f\n", x, y, z, r);
    fprintf(restart_fp, "%f %f %f %i %f\n", x, y, z, spin, r);
  }

  if (params.calc == DYNAMICS)
  {
    fprintf(restart_fp, "@nuc_velocities\n");
    double vx, vy, vz, vr;
    for (i = 0; i < NumNuclei(); i++)
    {
      GetNuclearVelocity(i, &vx, &vy, &vz);
      fprintf(restart_fp, "%i\t%f\t%f\t%f\n", i+1, vx, vy, vz);
    }

    fprintf(restart_fp, "@elec_velocities\n");
    for (i = 0; i < NumElectrons(); i++)
    {
      GetElectronVelocity(i, &vx, &vy, &vz, &vr);
      fprintf(restart_fp, "%i\t%f\t%f\t%f\t%f\n", i+1, vx, vy, vz, vr);
    }

    int i, idx;
    double mx, my, mz, mr;
  
    /* List explicitly specified masses if any */
    if (GetNumNucM() > 0)
    {
      fprintf(restart_fp, "@nuc_masses\n");
      for (i = 0; i < GetNumNucM(); i++)
      {
        GetNucM(i, &idx, &mx, &my, &mz);
        fprintf(restart_fp, "%i\t%e\t%e\t%e\n", idx+1, mx, my, mz);
      }
    }
    if (GetNumElecM() > 0)
    {
      fprintf(restart_fp, "@elec_masses\n");
      for (i = 0; i < GetNumElecM(); i++)
      {
        GetElecM(i, &idx, &mx, &my, &mz, &mr);
        fprintf(restart_fp, "%i\t%e\t%e\t%e\t%e\n", idx+1, mx, my, mz, mr);
      }
    }
  }
  fclose(restart_fp);
}

void OutputFixedHeader(FILE *out_fp)
{
  /* Number of nuclei and electrons */
  fprintf(out_fp, "[num_nuclei]\t%i\n", NumNuclei());
  fprintf(out_fp, "[num_electrons]\t%i\n", NumElectrons());

  /* Nucleus charge and label and mass */
  fprintf(out_fp, "[nuc_header]\tindex\tZ\tlabel\n");
  int i;
  double Z;
  for (i = 0; i < NumNuclei(); i++)
  {
    Z = GetNuclearCharge(i);
    fprintf(out_fp, "[nuc]\t%i\t%f\t%s\n", i+1, Z, GetType((int) Z));
  }

  /* Electron spin */
  fprintf(out_fp, "[elec_header]\tindex\tspin\n");
  for (i = 0; i < NumElectrons(); i++)
  {
    fprintf(out_fp, "[elec]\t%i\t%i\n", i+1, GetElectronSpin(i));
  }
}

void OutputDynamicsMasses(FILE *out_fp)
{
  int i, idx;
  double mx, my, mz, mr;
  
  /* List explicitly specified masses if any */
  if (GetNumNucM() > 0)
  {
    fprintf(out_fp, "[nuc_mass_header]\tmx\tmy\tmz\n");
    for (i = 0; i < GetNumNucM(); i++)
    {
      GetNucM(i, &idx, &mx, &my, &mz);
      fprintf(out_fp, "[nuc_mass]\t%i\t%f\t%f\t%f\n", idx+1, mx, my, mz);
    }
  }
  if (GetNumElecM() > 0)
  {
    fprintf(out_fp, "[elec_mass_header]\tmx\tmy\tmz\tmr\n");
    for (i = 0; i < GetNumElecM(); i++)
    {
      GetElecM(i, &idx, &mx, &my, &mz, &mr);
      fprintf(out_fp, "[elec_mass]\t%i\t%f\t%f\t%f\t%f\n", idx+1, mx, my, mz, mr);
    }
  }
}

void OutputParams(char *tag, FILE *out_fp)
{
  /* Cutoffs */
  fprintf(out_fp, "%s\ts_cutoff = %f\n", tag, params.s_cutoff);

  /* Boundaries */
  char bool_string[128];
  fprintf(out_fp, "%s\tx_bound = %f %f\n", tag, params.x_bound[0], params.x_bound[1]);
  fprintf(out_fp, "%s\ty_bound = %f %f\n", tag, params.y_bound[0], params.y_bound[1]);
  fprintf(out_fp, "%s\tz_bound = %f %f\n", tag, params.z_bound[0], params.z_bound[1]);
  BoolToString(params.periodic, bool_string);
  fprintf(out_fp, "%s\tperiodic = %s\n", tag, bool_string);

  /* Ewald */
  fprintf(out_fp, "%s\tewald_log_precision = %f\n", tag, params.ewald_log_precision);
  fprintf(out_fp, "%s\tewald_re_cutoff = %f\n", tag, params.ewald_re_cutoff);
  fprintf(out_fp, "%s\tewald_max_re = %f\n", tag, params.ewald_max_re);
  BoolToString(params.ewald_autoset, bool_string);
  fprintf(out_fp, "%s\tewald_autoset = %s\n", tag, bool_string);
  fprintf(out_fp, "%s\tewald_r_cutoff = %f\n", tag, params.ewald_r_cutoff);
  fprintf(out_fp, "%s\tewald_k_cutoff = %f\n", tag, params.ewald_k_cutoff);
  fprintf(out_fp, "%s\tewald_nuc_r = %e\n", tag, params.ewald_nuc_r);

  /* Minimization */
  char min_str[128];
  MinimizeTypeToString(params.min, min_str);
  fprintf(out_fp, "%s\tmin = %s\n", tag, min_str);

  /* Dynamics */
  char thermo_str[128];
  ThermostatTypeToString(params.thermostat, thermo_str);
  fprintf(out_fp, "%s\tthermostat = %s\n", tag, thermo_str);

  fprintf(out_fp, "%s\tstart_temperature = %f\n", tag, params.start_temperature);
  fprintf(out_fp, "%s\tdt = %f\n", tag, params.dt);

  /* Updating */
  fprintf(out_fp, "%s\tprint_every = %i\n", tag, params.print_every);
  fprintf(out_fp, "%s\tnum_steps = %i\n", tag, params.num_steps);

  /* Action */
  char action_str[128];
  ActionTypeToString(params.calc, action_str);
  fprintf(out_fp, "%s\tcalc = %s\n", tag, action_str);

  /* Output */
  char output_str[128];
  OutputTypeToString(params.output_position, output_str);
  fprintf(out_fp, "%s\toutput_position = %s\n", tag, output_str);

  OutputTypeToString(params.output_velocity, output_str);
  fprintf(out_fp, "%s\toutput_velocity = %s\n", tag, output_str);

  OutputTypeToString(params.output_energy_forces, output_str);
  fprintf(out_fp, "%s\toutput_energy_forces = %s\n", tag, output_str);

  OutputTypeToString(params.output_restart, output_str);
  fprintf(out_fp, "%s\toutput_restart = %s\n", tag, output_str);

  fprintf(out_fp, "%s\telectron_mass = %f\n", tag, params.electron_mass);

  /* Random numbers */
  fprintf(out_fp, "%s\trand_seed = %i\n", tag, params.rand_seed);

  /* External electric field */
  fprintf(out_fp, "%s\te_field = %f %f %f\n", tag, params.Ex, params.Ey, params.Ez);
  fprintf(out_fp, "%s\te_field_freq = %f\n", tag, params.Efreq);
  fprintf(out_fp, "%s\te_field_duration = %f\n", tag, params.Epacket_duration);

  /* Electron collapse */
  fprintf(out_fp, "%s\tcollapse_move = %f\n", tag, params.collapse_move);
}

void MinimizeTypeToString(enum MinimizeType min_type, char *str)
{
  if (min_type == CONJUGATE_GRADIENT)
    sprintf(str, "conjugate_gradient");
  else if (min_type == NEWTON)
    sprintf(str, "newton");
  else
    error("Minimize type %s not recognized.\n", str);
}

void ThermostatTypeToString(enum ThermostatType thermo_type, char *str)
{
  if (thermo_type == NONE)
    sprintf(str, "none");
  else if (thermo_type == NOSE_HOOVER)
    sprintf(str, "nose_hoover");
  else
    error("Thermostat type not recognized.\n");
}

void ActionTypeToString(enum ActionType action_type, char *str)
{
  if (action_type == SINGLE_PT)
    sprintf(str, "single_pt");
  else if (action_type == MINIMIZE)
    sprintf(str, "minimize");
  else if (action_type == DYNAMICS)
    sprintf(str, "dynamics");
  else
    error("Calculation type not recognized.\n");
}

void OutputTypeToString(enum OutputType output_type, char *str)
{
  if (output_type == NONE2)
    sprintf(str, "none");
  else if (output_type == ALL)
    sprintf(str, "all");
  else if (output_type == END)
    sprintf(str, "end");
  else
    error("Output type not recognized.\n");
}

void BoolToString(int val, char *str)
{
  if (val == 1)
    sprintf(str, "true");
  else
    sprintf(str, "false");
}

void PrintAnalyticForces()
{
	int i;
	double fx, fy, fz, fr;
	for (i = 0; i < NumNuclei(); i++)
	{
		GetNuclearForce(i, &fx, &fy, &fz);
		printf("n%i %e %e %e\n", i, fx, fy, fz);
	}
	for (i = 0; i < NumElectrons(); i++)
	{
		GetElectronForce(i, &fx, &fy, &fz, &fr);
		printf("e%i %e %e %e %e\n", i, fx, fy, fz, fr);
	}
}

void PrintNumericalForces()
{
	double en1, en2, dt = 0.00001;
	double fx, fy, fz, fr;
	int i;
  double energy;
  UpdateEnergyForces(); energy = GetTotalPE();

  // Calculate left handed, right handed, and centered derivative.

  double fx1, fy1, fz1, fr1, fx2, fy2, fz2, fr2;

	for (i = 0; i < NumNuclei(); i++)
	{
		nuc[i].x += dt; UpdateEnergyForces(); ApplyExternalField(0); en1 = GetTotalPE();
		nuc[i].x -= 2 * dt; UpdateEnergyForces(); ApplyExternalField(0); en2 = GetTotalPE();
		nuc[i].x += dt;
    fx1 = (en2 - energy) / dt;
    fx2 = (energy - en1) / dt;
		fx = -(en1 - en2) / (2 * dt);

		nuc[i].y += dt; UpdateEnergyForces(); ApplyExternalField(0); en1 = GetTotalPE();
		nuc[i].y -= 2 * dt; UpdateEnergyForces(); ApplyExternalField(0); en2 = GetTotalPE();
		nuc[i].y += dt;
    fy1 = (en2 - energy) / dt;
    fy2 = (energy - en1) / dt;
		fy = -(en1 - en2) / (2 * dt);

		nuc[i].z += dt; UpdateEnergyForces(); ApplyExternalField(0); en1 = GetTotalPE();
		nuc[i].z -= 2 * dt; UpdateEnergyForces(); ApplyExternalField(0); en2 = GetTotalPE();
		nuc[i].z += dt;
    fz1 = (en2 - energy) / dt;
    fz2 = (energy - en1) / dt;
		fz = -(en1 - en2) / (2 * dt);
		printf("n%i %e %e %e\n", i, fx, fy, fz);
    //printf(":n%i %e %e %e\n", i, fx1, fy1, fz1);
    //printf(":n%i %e %e %e\n", i, fx2, fy2, fz2);
	}

	for (i = 0; i < NumElectrons(); i++)
	{
		elec[i].x += dt; UpdateEnergyForces(); ApplyExternalField(0); en1 = GetTotalPE();
		elec[i].x -= 2 * dt; UpdateEnergyForces(); ApplyExternalField(0); en2 = GetTotalPE();
		elec[i].x += dt;
    fx1 = (en2 - energy) / dt;
    fx2 = (energy - en1) / dt;
		fx = -(en1 - en2) / (2 * dt);

		elec[i].y += dt; UpdateEnergyForces(); ApplyExternalField(0); en1 = GetTotalPE();
		elec[i].y -= 2 * dt; UpdateEnergyForces(); ApplyExternalField(0); en2 = GetTotalPE();
		elec[i].y += dt;
    fy1 = (en2 - energy) / dt;
    fy2 = (energy - en1) / dt;
		fy = -(en1 - en2) / (2 * dt);

		elec[i].z += dt; UpdateEnergyForces(); ApplyExternalField(0); en1 = GetTotalPE();
		elec[i].z -= 2 * dt; UpdateEnergyForces(); ApplyExternalField(0); en2 = GetTotalPE();
		elec[i].z += dt;
    fz1 = (en2 - energy) / dt;
    fz2 = (energy - en1) / dt;
		fz = -(en1 - en2) / (2 * dt);

		elec[i].r += dt; UpdateEnergyForces(); ApplyExternalField(0); en1 = GetTotalPE();
		elec[i].r -= 2 * dt; UpdateEnergyForces(); ApplyExternalField(0); en2 = GetTotalPE();
		elec[i].r += dt;
    fr1 = (en2 - energy) / dt;
    fr2 = (energy - en1) / dt;
		fr = -(en1 - en2) / (2 * dt);
		printf("e%i %e %e %e %e\n", i, fx, fy, fz, fr);
    //printf(":e%i %e %e %e %e\n", i, fx1, fy1, fz1, fr1);
    //printf(":e%i %e %e %e %e\n", i, fx2, fy2, fz2, fr2);
	}
}
