// Make contribution to which virial (rigid electron, compressible electron, both)
enum ElectronType {RIGID, FLEXIBLE, BOTH};
void ClearVirial();
void AddForceVirial(double rx, double ry, double rz, double rr, double fx, double fy, double fz, double fr, enum ElectronType electron_type);
void AddRadialForceVirial(double r, double f, enum ElectronType electron_type);
void AddSizeForceVirial(double r, double f, enum ElectronType electron_type);
void AddPotentialEnergyVirial(double energy, enum ElectronType electron_type);
void AddKineticEnergyVirial(double energy, enum ElectronType electron_type);
double GetRigidPEPressure(double volume);
double GetFlexiblePEPressure(double volume);
double GetRigidKEPressure(double volume);
double GetFlexibleKEPressure(double volume);
