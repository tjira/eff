#include <stdlib.h>

void AllocateInitialVelocities(int s_max_nuc_v, int s_max_elec_v);
void AllocateInitialMasses(int s_max_nuc_m, int s_max_elec_m);
int GetNumNucV();
int GetNumElecV();
int GetNumNucM();
int GetNumElecM();
void AddNucV(int i, double vx, double vy, double vz);
void AddElecV(int i, double vx, double vy, double vz, double vr);
void GetNucV(int idx, int *i, double *vx, double *vy, double *vz);
void GetElecV(int idx, int *i, double *vx, double *vy, double *vz, double *vr);
void AddNucM(int i, double mx, double my, double mz);
void AddElecM(int i, double mx, double my, double mz, double mr);
void GetNucM(int idx, int *i, double *mx, double *my, double *mz);
void GetElecM(int idx, int *i, double *mx, double *my, double *mz, double *mr);
