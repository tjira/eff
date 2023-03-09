#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "eff_erf.h"
#include "eff_ewald.h"
#include "eff_global.h"
#include "eff_update.h"
#include "eff_constants.h"
#include "eff_pressure.h"

int numcharges;

double *ewald_x,  *ewald_y,  *ewald_z, *ewald_alpha, *ewald_q;  // charge position, exponent, and magnitude
double *ewald_fx, *ewald_fy, *ewald_fz, *ewald_falpha, *ewald_energy;  // forces and energies of charges
double *ewald_S_real, *ewald_S_imag, *ewald_S_cap_real, *ewald_S_cap_imag;   // structure factor

double *lattice_kx, *lattice_ky, *lattice_kz, *lattice_k2;
double *lattice_rx, *lattice_ry, *lattice_rz;

double ewald_L[3];
double ewald_max_alpha, ewald_max_r, ewald_max_k, ewald_nuc_alpha;
double ewald_total_charge;

int n_total_klattice, n_total_rlattice;

void AllocateEwald()
{
  numcharges = NumNuclei() + NumElectrons();

  ewald_x      = (double *) malloc(sizeof(double) * numcharges);
  ewald_y      = (double *) malloc(sizeof(double) * numcharges);
  ewald_z      = (double *) malloc(sizeof(double) * numcharges);
  ewald_alpha  = (double *) malloc(sizeof(double) * numcharges);
  ewald_q      = (double *) malloc(sizeof(double) * numcharges);

  ewald_fx     = (double *) malloc(sizeof(double) * numcharges);
  ewald_fy     = (double *) malloc(sizeof(double) * numcharges);
  ewald_fz     = (double *) malloc(sizeof(double) * numcharges);
  ewald_falpha = (double *) malloc(sizeof(double) * numcharges);
  ewald_energy = (double *) malloc(sizeof(double) * numcharges);

  ewald_S_real = (double *) malloc(sizeof(double) * numcharges);
  ewald_S_imag = (double *) malloc(sizeof(double) * numcharges);
  ewald_S_cap_real = (double *) malloc(sizeof(double) * numcharges);
  ewald_S_cap_imag = (double *) malloc(sizeof(double) * numcharges);

	if (ewald_x == 0 || ewald_y == 0 || ewald_z == 0 || ewald_alpha == 0 || ewald_q == 0 ||
      ewald_fx == 0 || ewald_fy == 0 || ewald_falpha == 0 ||
      ewald_S_real == 0 || ewald_S_imag == 0 || ewald_S_cap_real == 0 || ewald_S_cap_imag == 0)
  error("Could not allocate space for ewald objects.\n");
}

void InitializeEwald(double Lx, double Ly, double Lz, double max_alpha, double max_r, double max_k, double nuc_alpha)
{
  ewald_L[0] = Lx;
  ewald_L[1] = Ly;
  ewald_L[2] = Lz;

  ewald_max_alpha = max_alpha;
  ewald_max_r = max_r;
  ewald_max_k = max_k;

  ewald_nuc_alpha = nuc_alpha;

  /* Check neutrality */
  SetEwaldCharges(nuc_alpha);
  ewald_total_charge = TotalCharge();

  /* Initialize k space and real space lattices */
  InitializeKSpace();
  InitializeRSpace();
}

void SetEwaldCharges(double nuc_alpha)
{
  int i, idx = 0;
  for (i = 0; i < NumNuclei(); i++)
  {
    ewald_x[idx] = nuc[i].x;
    ewald_y[idx] = nuc[i].y;
    ewald_z[idx] = nuc[i].z;
    ewald_alpha[idx] = 2 * nuc_alpha;                 // density = |wavefunction|^2
    ewald_q[idx] = nuc[i].q;
    idx++;
  }

  for (i = 0; i < NumElectrons(); i++)
  {
    ewald_x[idx] = elec[i].x;
    ewald_y[idx] = elec[i].y;
    ewald_z[idx] = elec[i].z;
    ewald_alpha[idx] = 2.0 / (elec[i].r * elec[i].r);   // density = |wavefunction|^2
    ewald_q[idx] = -1;
    idx++;
  }
}

double EwaldRCutoff(double precision, double alpha_cutoff, double min_alpha)
{
  return sqrt((-log(10.0) * precision + 3) / (alpha_cutoff / (1 + alpha_cutoff / min_alpha)));
}

double EwaldKCutoff(double precision, double alpha_cutoff)
{
  return sqrt(4 * alpha_cutoff * (-log(10.0) * precision + 5));
}

double TotalCharge()
{
  int i;
  double total_q = 0;
  for (i = 0; i < numcharges; i++)
  {
    total_q += ewald_q[i];
  }
  return total_q;
}

double TotalVolume()
{
  return ewald_L[0] * ewald_L[1] * ewald_L[2];
}

void AddEwaldEnergyForce(double *fx, double *fy, double *fz, double *falpha, double *energy)
{
  int i, idx = 0;
  for (i = 0; i < NumNuclei(); i++)
  {
    nuc[i].fx += fx[idx];
    nuc[i].fy += fy[idx];
    nuc[i].fz += fz[idx];
    nuc[i].energy += energy[idx];
    idx++;
  }

  for (i = 0; i < NumElectrons(); i++)
  {
    elec[i].fx += fx[idx];
    elec[i].fy += fy[idx];
    elec[i].fz += fz[idx];
    elec[i].fr += falpha[idx] * -4 / (elec[i].r * elec[i].r * elec[i].r);
    elec[i].energy += energy[idx];
    idx++;
  }
}

void EwaldVirial(double *ewald_energy, double *ewald_falpha)
{
  /* Update virial contribution */
  /* For coulombic systems, virial = -U */
  int i;
  double total_elec_energy = 0;
  for (i = 0; i < numcharges; i++)
    total_elec_energy += ewald_energy[i];
  AddPotentialEnergyVirial(total_elec_energy, BOTH);

  /* Subtract away contribution from changing electron size */
  int idx = NumNuclei();
  for (i = 0; i < NumElectrons(); i++)
  {
    AddSizeForceVirial(elec[i].r, -ewald_falpha[idx] * -4 / (elec[i].r * elec[i].r * elec[i].r), RIGID);
    idx++;
  }
}

void UpdateEwaldEnergy()
{
  SetEwaldCharges(ewald_nuc_alpha);

  int i;
  for (i = 0; i < numcharges; i++)
    ewald_fx[i] = ewald_fy[i] = ewald_fz[i] = ewald_falpha[i] = ewald_energy[i] = 0;
  
  SelfEnergy(ewald_falpha, ewald_energy);
  AddEwaldEnergyForce(ewald_fx, ewald_fy, ewald_fz, ewald_falpha, ewald_energy);
  EwaldVirial(ewald_energy, ewald_falpha);

  if (ewald_total_charge != 0)
  {
    UniformChargeEnergy(ewald_falpha, ewald_energy);
    AddEwaldEnergyForce(ewald_fx, ewald_fy, ewald_fz, ewald_falpha, ewald_energy);
    EwaldVirial(ewald_energy, ewald_falpha);
  }

  KSpaceEnergy(ewald_fx, ewald_fy, ewald_fz, ewald_falpha, ewald_energy);
  AddEwaldEnergyForce(ewald_fx, ewald_fy, ewald_fz, ewald_falpha, ewald_energy);
  EwaldVirial(ewald_energy, ewald_falpha);

  RSpaceEnergy(ewald_fx, ewald_fy, ewald_fz, ewald_falpha, ewald_energy);
  AddEwaldEnergyForce(ewald_fx, ewald_fy, ewald_fz, ewald_falpha, ewald_energy);
  EwaldVirial(ewald_energy, ewald_falpha);
 
  return;
}

void InitializeKSpace()
{
  double dk[3];
  int n_klattice[3];

  int i;
  for (i = 0; i < 3; i++)
  {
    dk[i] = 2 * PI / ewald_L[i];
    n_klattice[i] = (int) (ewald_max_k / dk[i]);
  }
 
  double max_k_squared = ewald_max_k * ewald_max_k;

  /* Count points in k-lattice */
  int nx, ny, nz;
  double ksquared;
  int idx = 0;
  for (nx = -n_klattice[0]; nx <= n_klattice[0]; nx++)
    for (ny = -n_klattice[1]; ny <= n_klattice[1]; ny++)
      for (nz = -n_klattice[2]; nz <= n_klattice[2]; nz++)
      {
        ksquared = pow(nx * dk[0], 2) + pow(ny * dk[1], 2) + pow(nz * dk[2], 2);
        if (ksquared < max_k_squared) idx++;
      }

  n_total_klattice = idx - 1;  /* exclude nx,ny,nz = 0 */
  lattice_kx = (double *) malloc(sizeof(double) * n_total_klattice);
  lattice_ky = (double *) malloc(sizeof(double) * n_total_klattice);
  lattice_kz = (double *) malloc(sizeof(double) * n_total_klattice);
  lattice_k2 = (double *) malloc(sizeof(double) * n_total_klattice);

  /* Create k-lattice */
  idx = 0;
  for (nx = -n_klattice[0]; nx <= n_klattice[0]; nx++)
    for (ny = -n_klattice[1]; ny <= n_klattice[1]; ny++)
      for (nz = -n_klattice[2]; nz <= n_klattice[2]; nz++)
      {
        if (nx != 0 || ny != 0 || nz != 0)
        {
          ksquared = pow(nx * dk[0], 2) + pow(ny * dk[1], 2) + pow(nz * dk[2], 2);
          if (ksquared < max_k_squared)
          {
            lattice_kx[idx] = nx * dk[0];
            lattice_ky[idx] = ny * dk[1];
            lattice_kz[idx] = nz * dk[2];
            lattice_k2[idx] = ksquared;
            idx++;
          }
        }
      }
}

void InitializeRSpace()
{
  /* Compute n_rlattice from maximum alpha allowed */
  double dr[3];
  int n_rlattice[3];

  int i;
  for (i = 0; i < 3; i++)
  {
    dr[i] = ewald_L[i];
    n_rlattice[i] = (int) (ewald_max_r / dr[i] + 1);
  }

  double max_r_squared = ewald_max_r * ewald_max_r;

  /* Count number of lattice points */
  int nx, ny, nz;
  n_total_rlattice = (2 * n_rlattice[0] + 1) * (2 * n_rlattice[1] + 1) * (2 * n_rlattice[2] + 1);
  lattice_rx = (double *) malloc(sizeof(double) * n_total_rlattice);
  lattice_ry = (double *) malloc(sizeof(double) * n_total_rlattice);
  lattice_rz = (double *) malloc(sizeof(double) * n_total_rlattice);

  /* Make the zeroth element the zero lattice vector */
  int idx = 0;
  double rsquared;
  lattice_rx[0] = lattice_ry[0] = lattice_rz[0] = 0; 
  idx++;

  for (nx = -n_rlattice[0]; nx <= n_rlattice[0]; nx++)
    for (ny = -n_rlattice[1]; ny <= n_rlattice[1]; ny++)
      for (nz = -n_rlattice[2]; nz <= n_rlattice[2]; nz++)
      {
        if (nx != 0 || ny != 0 || nz != 0)
        {
          rsquared = pow(nx * dr[0], 2) + pow(ny * dr[1], 2) + pow(nz * dr[2], 2);
          lattice_rx[idx] = nx * dr[0];
          lattice_ry[idx] = ny * dr[1];
          lattice_rz[idx] = nz * dr[2];
          idx++;
        }
      }
}

void KSpaceEnergy(double *fx, double *fy, double *fz, double *falpha, double *energy)
{
  int i;
  for (i = 0; i < numcharges; i++)
    fx[i] = fy[i] = fz[i] = falpha[i] = energy[i] = 0;

  /* Calculate Ewald sum */
  int j;
  double ewald_energy = 0;
  for (i = 0; i < n_total_klattice; i++)
  {
    double rho_k_real = 0, rho_k_imag = 0;
    double rho_k_cap_real = 0, rho_k_cap_imag = 0;
    double angle, mag, mag_cap;
   
    for (j = 0; j < numcharges; j++)
    {
      angle   = -(lattice_kx[i] * ewald_x[j] + lattice_ky[i] * ewald_y[j] + lattice_kz[i] * ewald_z[j]);
      mag     = ewald_q[j] * exp(-lattice_k2[i] / (4 * ewald_alpha[j]));
      ewald_S_real[j] = mag * cos(angle);
      ewald_S_imag[j] = mag * sin(angle);

      if (ewald_alpha[j] > ewald_max_alpha)
      {
        mag_cap = ewald_q[j] * exp(-lattice_k2[i] / (4 * ewald_max_alpha));
        ewald_S_cap_real[j] = mag_cap * cos(angle);
        ewald_S_cap_imag[j] = mag_cap * sin(angle);
      }
      else
      {
        mag_cap = mag;
        ewald_S_cap_real[j] = ewald_S_real[j];
        ewald_S_cap_imag[j] = ewald_S_imag[j];
      }
 
      rho_k_real += ewald_S_real[j]; rho_k_cap_real += ewald_S_cap_real[j];
      rho_k_imag += ewald_S_imag[j]; rho_k_cap_imag += ewald_S_cap_imag[j];
    }

    /* Calculate forces and energy */
    double fxyz_mag;
    for (j = 0; j < numcharges; j++)
    {
      fxyz_mag = (ewald_S_real[j] * rho_k_cap_imag - ewald_S_imag[j] * rho_k_cap_real + ewald_S_cap_real[j] * rho_k_imag - ewald_S_cap_imag[j] * rho_k_real) / lattice_k2[i];
      fx[j] += fxyz_mag * lattice_kx[i];
      fy[j] += fxyz_mag * lattice_ky[i];
      fz[j] += fxyz_mag * lattice_kz[i];

      if (ewald_alpha[j] > ewald_max_alpha)
      {
        falpha[j] += (ewald_S_real[j] * rho_k_cap_real + ewald_S_imag[j] * rho_k_cap_imag) / (ewald_alpha[j] * ewald_alpha[j]);
      }
      else
      {
        falpha[j] += (ewald_S_real[j] * rho_k_cap_real + ewald_S_imag[j] * rho_k_cap_imag + 
                           ewald_S_cap_real[j] * rho_k_real + ewald_S_cap_imag[j] * rho_k_imag) / 
                          (ewald_alpha[j] * ewald_alpha[j]);
      }
    }
    ewald_energy += (1.0 / lattice_k2[i]) * (rho_k_real * rho_k_cap_real + rho_k_imag * rho_k_cap_imag);
  }

  /* Smear long range contribution uniformly over all charges */
  double V = TotalVolume();
  ewald_energy *= (2 * PI / V) / (double) numcharges;
  for (i = 0; i < numcharges; i++)
    energy[i] = ewald_energy;

  /* Multiply forces and energies by appropriate factors */
  double fxyz_factor, falpha_factor;
  fxyz_factor   = 2 * PI / V;
  falpha_factor = -0.5 * PI / V; 
  for (i = 0; i < numcharges; i++)
  {
    fx[i] *= fxyz_factor;
    fy[i] *= fxyz_factor;
    fz[i] *= fxyz_factor;
    falpha[i] *= falpha_factor;
  }
  return;
}

void SelfEnergy(double *falpha, double *energy)
{
  int i;
  for (i = 0; i < numcharges; i++)
  {
    double a1, a2;
    if (ewald_alpha[i] >= ewald_max_alpha)
    {
      a2 = sqrt(ewald_alpha[i] * ewald_max_alpha / (ewald_alpha[i] + ewald_max_alpha)); 
      energy[i] = -ewald_q[i] * ewald_q[i] * a2 / sqrt(PI);
      falpha[i] = ewald_q[i] * ewald_q[i] * (0.5 * a2 * a2 * a2 / (ewald_alpha[i] * ewald_alpha[i])) / sqrt(PI);
    }
    else
    {
//      a1 = sqrt(ewald_alpha[i] / 2);
      a1 = sqrt(ewald_alpha[i] * ewald_alpha[i] / (ewald_alpha[i] + ewald_alpha[i]));
      energy[i] = -ewald_q[i] * ewald_q[i] * a1 / sqrt(PI);
      falpha[i] = ewald_q[i] * ewald_q[i] / (4 * a1) / sqrt(PI);
    }
  }
}

void UniformChargeEnergy(double *falpha, double *energy)
{
  /* E = 1/4 Q (Pi / V) * sum((r_cap^2 - re^2) * q_i) */
  int i;
  double factor = -0.5 * ewald_total_charge * PI / TotalVolume();

  for (i = 0; i < numcharges; i++)
  {
    if (ewald_alpha[i] >= ewald_max_alpha)
    {
      energy[i] = factor * ewald_q[i] * (1.0 / ewald_max_alpha - 1.0 / ewald_alpha[i]);
      falpha[i] = -factor * ewald_q[i] / (ewald_alpha[i] * ewald_alpha[i]);
    }
    else
    {
      energy[i] = 0;
      falpha[i] = 0;
    }
  }
}

void RSpaceEnergy(double *fx, double *fy, double *fz, double *falpha, double *energy)
{
  /* Calculate interaction of unit cell with other charges on real lattice */
  int i;
  for (i = 0; i < numcharges; i++)
    fx[i] = fy[i] = fz[i] = falpha[i] = energy[i] = 0;

  double dx, dy, dz, r, a1, a2;
  int j, k;
  double factor;
  int idx = 0;
  for (i = 0; i < numcharges; i++)
  {
    for (j = 0; j < numcharges; j++)
    {
      if (ewald_alpha[j] > ewald_max_alpha)
      {
        a1  = sqrt(ewald_alpha[i] * ewald_alpha[j] / (ewald_alpha[i] + ewald_alpha[j])); 
        a2  = sqrt(ewald_alpha[i] * ewald_max_alpha / (ewald_alpha[i] + ewald_max_alpha)); 

        for (k = 0; k < n_total_rlattice; k++)
        { 
          if (i != j || k != 0)
          {
            dx = ewald_x[i] - ewald_x[j] - lattice_rx[k];
            dy = ewald_y[i] - ewald_y[j] - lattice_ry[k];
            dz = ewald_z[i] - ewald_z[j] - lattice_rz[k];
            r  = sqrt(dx * dx + dy * dy + dz * dz);

            if (r < ewald_max_r) 
            {
              double f1, df1, f2, df2;
              f1 = erfoverx1(r * a1, &df1);
              f2 = erfoverx1(r * a2, &df2);

              energy[i] += ewald_q[i] * ewald_q[j] * (a1 * f1 - a2 * f2);

              double mag_fxyz = ewald_q[i] * ewald_q[j] * (a1 * a1 * df1 - a2 * a2 * df2);
              double forcex, forcey, forcez;
              SmallRForce(dx, dy, dz, r, mag_fxyz, &forcex, &forcey, &forcez);

              fx[i] += forcex;
              fy[i] += forcey;
              fz[i] += forcez;
              fx[j] -= forcex;
              fy[j] -= forcey;
              fz[j] -= forcez;
          
              falpha[i] += ewald_q[i] * ewald_q[j] * (0.5 / (ewald_alpha[i] * ewald_alpha[i])) * (a1 * a1 * a1 * (r * a1 * df1 + f1) - a2 * a2 * a2 * (r * a2 * df2 + f2));
              falpha[j] += ewald_q[i] * ewald_q[j] * (0.5 / (ewald_alpha[j] * ewald_alpha[j])) * (a1 * a1 * a1 * (r * a1 * df1 + f1));
              idx++;
            }
          }
        }
      }
    }
  }

  for (i = 0; i < numcharges; i++)
  {
    fx[i] *= -0.5;
    fy[i] *= -0.5;
    fz[i] *= -0.5;
    falpha[i] *= -0.5;
    energy[i] *= 0.5;
  }

  return;
}
