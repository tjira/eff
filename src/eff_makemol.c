#include <stdio.h>
#include <math.h>
#include "eff_create.h"
#include "eff_global.h"
#include "eff_makemol.h"
#include "eff_constants.h"
#include "eff_util.h"

MMXYZ *mmxyz;
int numxyz;

MMLIB mmlib[100];
int numlib = 0;

void AllocateMakemol(int numatoms)
{
  mmxyz = (MMXYZ *) malloc(sizeof(MMXYZ) * numatoms);
}

void XYZMakemol(FILE *fp)
{
	char str[2048], s1[2048];
	double d1, d2, d3, d4;

	if (fp == NULL) error("File not found.\n");
  int lines = 0;
  /* Read in number of atoms */
  fgets(str, 2048, fp);
  sscanf(str, "%i", &numxyz);
  AllocateMakemol(numxyz);
  /* Read in comment line */
  fgets(str, 2048, fp);
  /* Read in XYZ coordinates */
  int i = 0;
	while(fgets(str, 2048, fp))
	{
    sscanf(str, "%s %lf %lf %lf", s1, &mmxyz[i].x, &mmxyz[i].y, &mmxyz[i].z);
    mmxyz[i].q = GetZ(s1);
    i++;
    if (i > numxyz) error("Number of atoms in XYZ file doesn't match number declared.\n");
  }
  /* Intitialize library */
  InitializeMakemol();
  /* Add electrons and nuclei appropriately */
  ExecuteMakemol(fp);
}

void InitializeMakemol()
{
  /* Z1, Z2, fractional position, radius of electron, bond length in A */
  numlib = 0;
	Store(1, 6, 1, 0.293, 1.543, 1.2);   // C-H bonds
	Store(6, 6, 1, 0.5, 1.258, 1.6);     // C-C, C=C, CtrC bonds.
	Store(6, 6, 2, 0.5, 1.348, 1.45);
	Store(6, 6, 3, 0.5, 1.619, 1.24);
}

void Store(int c1, int c2, int bo, double x, double r, double length)
{
	mmlib[numlib].c1 = c1;
	mmlib[numlib].c2 = c2;
	mmlib[numlib].bo = bo;
	mmlib[numlib].x = x;
	mmlib[numlib].r = r;
  mmlib[numlib].length = length;
	numlib++;
  if (c1 < c2)
    Store(c2, c1, bo, 1 - x, r, length); 
}

double CoreRadius(double Z)
{
  return 3 * sqrt(PI) / (2 * sqrt(2.0) * Z - 1.0);
}

void ExecuteMakemol(FILE *fp)
{
	// Output nuclei.
  int i, j, k;
	for (i = 0; i < numxyz; i++)
    AddNucleus(mmxyz[i].x / A0, mmxyz[i].y / A0, mmxyz[i].z / A0, mmxyz[i].q);

  // Output electrons
  int numelecs = 0;

  // Core electrons on non-H nuclei
	for (i = 0; i < numxyz; i++)
	{
		if (mmxyz[i].q > 1)
		{
      double r = CoreRadius(mmxyz[i].q);
      AddElectron(mmxyz[i].x / A0, mmxyz[i].y / A0, mmxyz[i].z / A0, 1, r);
      AddElectron(mmxyz[i].x / A0, mmxyz[i].y / A0, mmxyz[i].z / A0, -1, r);
			numelecs += 2;
		}
	}

  // Assign an electron pair to each chemical bond
  MMLIB *lib;
  double dx, dy, dz, r;
  for (i = 0; i < numxyz; i++)
  {
    for (j = 0; j < i; j++)
    {
      dx = mmxyz[i].x - mmxyz[j].x;
      dy = mmxyz[i].y - mmxyz[j].y;
      dz = mmxyz[i].z - mmxyz[j].z;
      r = sqrt(dx * dx + dy * dy + dz * dz);

      if (r < 2)
      {
        lib = GetBond(mmxyz[i].q, mmxyz[j].q, r);
        if (lib) /* found a bond */
        {
          int bo = lib->bo;
          double frac = lib->x;
          double re = lib->r;
          double xe, ye, ze;
          xe = frac * mmxyz[i].x + (1 - frac) * mmxyz[j].x;
          ye = frac * mmxyz[i].y + (1 - frac) * mmxyz[j].y;
          ze = frac * mmxyz[i].z + (1 - frac) * mmxyz[j].z;
          for (k = 0; k < bo; k++)
          {
            /* one electron pair per bond order */
            if (k > 0)
            {
              /* shift electrons around */
               xe += 0.2 * rand_uni() - 0.1;
			         ye += 0.2 * rand_uni() - 0.1;
			         ze += 0.2 * rand_uni() - 0.1;
            }
            AddElectron(xe / A0, ye / A0, ze / A0, 1, re);
            AddElectron(xe / A0, ye / A0, ze / A0, -1, re);
            numelecs += 2;
          }
        }
      }
    }
  }

  int sumZ = 0;
  for (i = 0; i < numxyz; i++)
    sumZ += mmxyz[i].q;
           
	if (numelecs != sumZ)
		printf("Warning: sum of nuclear charges (%i) != number of electrons (%i)\n", sumZ, numelecs);
}

MMLIB *GetBond(int q1, int q2, double r)
{
  int i;
  double min_r = 10000;
  MMLIB *match_lib = 0;

  /* Search to find type of bond connecting two charges */
  /* We want to find the smallest bond that is bigger than ours */
  for (i = 0; i < numlib; i++)
  {
    if (mmlib[i].c1 == q1 && mmlib[i].c2 == q2)
    {
      if (mmlib[i].length > r && mmlib[i].length < min_r)
      {
        min_r = mmlib[i].length;
        match_lib = &(mmlib[i]);
      }
    }
  }
  return match_lib;
}
