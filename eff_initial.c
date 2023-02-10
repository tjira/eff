#include "eff_global.h"
#include "eff_initial.h"
#include "eff_util.h"

/* Variables and functions to store initial velocities */
int num_nuc_v = 0, num_elec_v = 0;
int max_nuc_v, max_elec_v;

/* Variables and functions to store initial masses */
int num_nuc_m = 0, num_elec_m = 0;
int max_nuc_m, max_elec_m;

typedef struct
{
  int i;
  double vx, vy, vz;
} NUC_V;

typedef struct
{
  int i;
  double vx, vy, vz, vr;
} ELEC_V;

typedef struct
{
  int i;
  double mx, my, mz;
} NUC_M;

typedef struct
{
  int i;
  double mx, my, mz, mr;
} ELEC_M;

NUC_V *nuc_v_array;
ELEC_V *elec_v_array;
NUC_M *nuc_m_array;
ELEC_M *elec_m_array;

void AllocateInitialVelocities(int s_max_nuc_v, int s_max_elec_v)
{
  max_nuc_v = s_max_nuc_v;
  max_elec_v = s_max_elec_v;

  nuc_v_array = (NUC_V *) malloc(sizeof(NUC_V) * max_nuc_v);
  elec_v_array = (ELEC_V *) malloc(sizeof(ELEC_V) * max_elec_v);
}

void AllocateInitialMasses(int s_max_nuc_m, int s_max_elec_m)
{
  max_nuc_m = s_max_nuc_m;
  max_elec_m = s_max_elec_m;

  nuc_m_array = (NUC_M *) malloc(sizeof(NUC_M) * max_nuc_m);
  elec_m_array = (ELEC_M *) malloc(sizeof(ELEC_M) * max_elec_m);
}

int GetNumNucV()
{
  return num_nuc_v;
}

int GetNumElecV()
{
  return num_elec_v;
}

int GetNumNucM()
{
  return num_nuc_m;
}

int GetNumElecM()
{
  return num_elec_m;
}

void AddNucV(int i, double vx, double vy, double vz)
{
  nuc_v_array[num_nuc_v].i = i;
  nuc_v_array[num_nuc_v].vx = vx;
  nuc_v_array[num_nuc_v].vy = vy;
  nuc_v_array[num_nuc_v].vz = vz;
  num_nuc_v++;
  if (num_nuc_v > max_nuc_v) error("More initial nuc velocities added than space allocated.\n");
}

void AddElecV(int i, double vx, double vy, double vz, double vr)
{
  elec_v_array[num_elec_v].i = i;
  elec_v_array[num_elec_v].vx = vx;
  elec_v_array[num_elec_v].vy = vy;
  elec_v_array[num_elec_v].vz = vz;
  elec_v_array[num_elec_v].vr = vr;
  num_elec_v++;
  if (num_elec_v > max_elec_v) error("More initial elec velocities added than space allocated.\n");
}

void GetNucV(int idx, int *i, double *vx, double *vy, double *vz)
{
  if (idx < 0 || idx >= num_nuc_v) error("Invalid nuc v index.\n");
  *i =  nuc_v_array[idx].i;
  *vx = nuc_v_array[idx].vx;
  *vy = nuc_v_array[idx].vy;
  *vz = nuc_v_array[idx].vz;
}

void GetElecV(int idx, int *i, double *vx, double *vy, double *vz, double *vr)
{
  if (idx < 0 || idx >= num_elec_v) error("Invalid elec v index.\n");
  *i =  elec_v_array[idx].i;
  *vx = elec_v_array[idx].vx;
  *vy = elec_v_array[idx].vy;
  *vz = elec_v_array[idx].vz;
  *vr = elec_v_array[idx].vr;
}

void AddNucM(int i, double mx, double my, double mz)
{
  nuc_m_array[num_nuc_m].i = i;
  nuc_m_array[num_nuc_m].mx = mx;
  nuc_m_array[num_nuc_m].my = my;
  nuc_m_array[num_nuc_m].mz = mz;
  num_nuc_m++;
  if (num_nuc_m > max_nuc_m) error("More initial nuc masses added than space allocated.\n");
}

void AddElecM(int i, double mx, double my, double mz, double mr)
{
  elec_m_array[num_elec_m].i = i;
  elec_m_array[num_elec_m].mx = mx;
  elec_m_array[num_elec_m].my = my;
  elec_m_array[num_elec_m].mz = mz;
  elec_m_array[num_elec_m].mr = mr;
  num_elec_m++;
  if (num_elec_m > max_elec_m) error("More initial elec masses added than space allocated.\n");
}

void GetNucM(int idx, int *i, double *mx, double *my, double *mz)
{
  if (idx < 0 || idx >= num_nuc_m) error("Invalid nuc m index.\n");
  *i =  nuc_m_array[idx].i;
  *mx = nuc_m_array[idx].mx;
  *my = nuc_m_array[idx].my;
  *mz = nuc_m_array[idx].mz;
}

void GetElecM(int idx, int *i, double *mx, double *my, double *mz, double *mr)
{
  if (idx < 0 || idx >= num_elec_m) error("Invalid elec m index.\n");
  *i =  elec_m_array[idx].i;
  *mx = elec_m_array[idx].mx;
  *my = elec_m_array[idx].my;
  *mz = elec_m_array[idx].mz;
  *mr = elec_m_array[idx].mr;
}


