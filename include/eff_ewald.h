void AllocateEwald();
void InitializeEwald(double Lx, double Ly, double Lz, double max_alpha, double max_r, double max_k, double nuc_alpha);
void SetEwaldCharges(double nuc_alpha);
double TotalCharge();
double TotalVolume();
void AddEwaldEnergyForce(double *fx, double *fy, double *fz, double *falpha, double *energy);
void EwaldVirial(double *ewald_energy, double *ewald_falpha);
void UpdateEwaldEnergy();
void InitializeKSpace();
void InitializeRSpace();
void KSpaceEnergy(double *fx, double *fy, double *fz, double *falpha, double *energy);
void SelfEnergy(double *falpha, double *energy);
void RSpaceEnergy(double *fx, double *fy, double *fz, double *falpha, double *energy);
void UniformChargeEnergy(double *falpha, double *energy);

double EwaldRCutoff(double precision, double alpha_cutoff, double min_alpha);
double EwaldKCutoff(double precision, double alpha_cutoff);

