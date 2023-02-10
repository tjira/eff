#include <stdio.h>

typedef struct
{
  double x, y, z;
  int q;
} MMXYZ;

typedef struct
{
	int c1, c2, bo;
	double x, r, length;
} MMLIB;

void AllocateMakemol(int numatoms);
void XYZMakemol(FILE *fp);
void InitializeMakemol();
void Store(int c1, int c2, int bo, double x, double r, double length);
double CoreRadius(double Z);
void ExecuteMakemol(FILE *fp);
MMLIB *GetBond(int q1, int q2, double r);

