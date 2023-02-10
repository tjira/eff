#include "eff_access.h"
#include "eff_update.h"
#include "eff_util.h"
#include "eff_global.h"

double GetTotalPE()
{
  int i;
	double total_energy = 0;
	for (i = 0; i < numnuclei; i++)
		total_energy += GetNucleusPE(i);
	for (i = 0; i < numelectrons; i++)
		total_energy += GetElectronPE(i);
  return total_energy;
}

double GetNucleusPE(int i)
{
	if (i < 0 || i >= numnuclei) error("Nucleus %i does not exist.\n", i);
	return nuc[i].energy;
}

double GetElectronPE(int i)
{
	if (i < 0 || i >= numelectrons) error("Electron %i does not exist.\n", i);
	return elec[i].energy;
}

void GetNuclearForce(int i, double *fx, double *fy, double *fz)
{
	if (i < 0 || i >= numnuclei) error("Nucleus %i does not exist.\n", i);
	*fx = nuc[i].fx;
	*fy = nuc[i].fy;
	*fz = nuc[i].fz;
}

void GetElectronForce(int i, double *fx, double *fy, double *fz, double *fr)
{
	if (i < 0 || i >= numelectrons) error("Electron %i does not exist.\n", i);
	*fx = elec[i].fx;
	*fy = elec[i].fy;
	*fz = elec[i].fz;
	*fr = elec[i].fr;
}

void GetElectronPosition(int i, double *x, double *y, double *z, double *r)
{
	*x = elec[i].x; *y = elec[i].y; *z = elec[i].z; *r = elec[i].r;
}

void SetElectronPosition(int i, double x, double y, double z, double r)
{
	elec[i].x = x; elec[i].y = y; elec[i].z = z; elec[i].r = r;
}

void SetElectronRadius(int i, double r)
{
	elec[i].r = r;
}

double GetElectronRadiusForce(int i)
{
	if (i < 0 || i >= numelectrons) error("Electron %i does not exist.\n", i);
        return elec[i].fr;
}

void GetNuclearPosition(int i, double *x, double *y, double *z)
{
	*x = nuc[i].x; *y = nuc[i].y; *z = nuc[i].z;
}

void SetNuclearPosition(int i, double x, double y, double z)
{
	nuc[i].x = x; nuc[i].y = y; nuc[i].z = z;
}

double GetNuclearCharge(int i)
{
	return nuc[i].q;
}

int GetElectronSpin(int i)
{
  return elec[i].spin;
}

int NumNuclei()
{
	return numnuclei;
}

int NumElectrons()
{
	return numelectrons;
}

void AddNuclearForce(int i, double fx, double fy, double fz)
{
	nuc[i].fx += fx;
	nuc[i].fy += fy;
	nuc[i].fz += fz;
}

void AddElectronForce(int i, double fx, double fy, double fz, double fr)
{
	elec[i].fx += fx;
	elec[i].fy += fy;
	elec[i].fz += fz;
	elec[i].fr += fr;
}

void AddNucleusPE(int i, double energy)
{
	nuc[i].energy += energy;
}

void AddElectronPE(int i, double energy)
{
	elec[i].energy += energy;
}

void ShiftMonomer(double r, double dx, double dy, double dz)
{
  int i;
  for (i = 0; i < 1; i++)
  {
    nuc[i].x += r * dx;
    nuc[i].y += r * dy;
    nuc[i].z += r * dz;
  }
/*
  for (i = 0; i < 2; i++)
  {
    elec[i].x += r * dx;
    elec[i].y += r * dy;
    elec[i].z += r * dz;
  }
*/
}

void ScaleElectrons(double r)
{
  int i;
  for (i = 2; i < numelectrons; i++)
  {
    elec[i].x *= r;
    elec[i].y *= r;
    elec[i].z *= r;
  }
}

int NumHeavyNuclei()
{
  int i, count = 0;
  for (i = 0; i < numnuclei; i++)
    if (nuc[i].q > 1) count++;
  return count;
}

