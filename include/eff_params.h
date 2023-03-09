#include <stdlib.h>
#include <string.h>

void InitializeParameters();
void AddParameter(char *name, void (*fp)(int, char **));
void ExecuteParameter(char *name, int argc, char **argv);

void SetDefaultParameters();

void ParamECutoff(int argc, char **argv);
void ParamSCutoff(int argc, char **argv);
void ParamXYZBound(int argc, char **argv);
void ReadValue(char *str, double *val);
void ParamPeriodic(int argc, char **argv);
void ParamEwald(int argc, char **argv);
void ParamMin(int argc, char **argv);
void ParamThermostat(int argc, char **argv);
void ParamStartTemperature(int argc, char **argv);
void ParamDt(int argc, char **argv);
void ParamPrintEvery(int argc, char **argv);
void ParamNumSteps(int argc, char **argv);
void ParamCalc(int argc, char **argv);
void ParamOutput(int argc, char **argv);
void ParamElectronMass(int argc, char **argv);
void ParamRandSeed(int argc, char **argv);
void ParamEField(int argc, char **argv);
void ParamEFieldFreq(int argc, char **argv);
void ParamEFieldPacketDuration(int argc, char **argv);
void ParamCollapseMove(int argc, char **argv);

