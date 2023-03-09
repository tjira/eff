#include <stdlib.h>

void IntitializeEnergyUpdater(double s_cutoff, int periodic);
void UpdateEnergyForces();
void UpdateKineticEnergy();
void UpdateElecAndPauli();
void UpdatePauliPeriodic();
void UpdateNucNuc(int i, int j, double energy, double fx, double fy, double fz);
void UpdateNucElec(int i, int j, double energy, double fx, double fy, double fz, double fre1);
void UpdateElecElec(int i, int j, double energy, double fx, double fy, double fz, double fre1, double fre2);
void RForce(double dx, double dy, double dz, double rc, double force, double *fx, double *fy, double *fz);
void SmallRForce(double dx, double dy, double dz, double rc, double force, double *fx, double *fy, double *fz);
void ElecNucNuc(double q, double rc, double *energy, double *frc);
void ElecNucElec(double q, double rc, double re1, double *energy, double *frc, double *fre1, int i, int j);
void ElecElecElec(double rc, double re1, double re2, double *energy, double *frc, double *fre1, double *fre2, int i, int j);
void PauliElecElec(int samespin, double rc, double re1, double re2, double *energy, double *frc, double *fre1, double *fre2, double a, double b, double c, double d, int i, int j);
double unit_cubic_spline(double x);
double d_unit_cubic_spline(double x);
double cubic_spline(double x, double *fx, double r_slope, double r_cut);
double cubic_spline2(double x, double *fx, double r_cut);
double quintic_spline(double x, double *fx, double r_slope, double r_cut, double rd_1, double rd_2, double rd_3);
double unit_quintic_spline(double x, double *fx, double x1, double y1, double dy1, double d2y1, double x2, double y2, double dy2, double d2y2);
double quintic_spline2(double x, double *fx, double y0, double r_cut);
void HOverrides();
double quintic_spline3(double x, double *fx, double scale_slope, double deriv2, double r_cut, double deriv22);
void SlaterPotential(double Z, double q, double rc, double *energy, double *frc);
double OverlapSum();
void TotalPauli(int i, int j, double rc, double re1, double re2, double *pauli_energy, double *pauli_frc, double *pauli_fre1, double *pauli_fre2);
