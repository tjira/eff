#include "eff_global.h"

void AllocateDynamics();
void InitializeDynamics();
void InitializeRandomVelocities(double temp);
void Dynamics();
void UpdateDynamicsPositions();
void UpdateDynamicsVelocities();
void GetDynamicsForces();
void GetDynamicsPositions();
void SetDynamicsPositions();
void UpdateDynamicsVirial();

/* dynamics accessor functions */
void Dynamics(enum ThermostatType s_thermostat, double temperature);
void SetTimeStep(double s_dt);
double GetTemperature();

void SetNuclearVelocity(int i, double vx, double vy, double vz);
void SetElectronVelocity(int i, double vx, double vy, double vz, double vr);
void GetNuclearVelocity(int i, double *vx, double *vy, double *vz);
void GetElectronVelocity(int i, double *vx, double *vy, double *vz, double *vr);

void SetNuclearMass(int i, double mx, double my, double mz);
void SetElectronMass(int i, double mx, double my, double mz, double mr);
void GetNuclearMass(int i, double *mx, double *my, double *mz);
void GetElectronMass(int i, double *mx, double *my, double *mz, double *mr);

double GetTotalKE();
double GetNucleusKE(int i);
double GetElectronKE(int i);
double GetElectronTranslationKE(int i);
