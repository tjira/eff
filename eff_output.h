#include <stdio.h>
#include "eff_global.h"

void OutputEnergyForces(FILE *out_fp);
void OutputEnergy(FILE *out_fp);
void OutputTimeElapsed(FILE *out_fp);
void OutputPositions(FILE *out_fp);
void OutputVelocities(FILE *out_fp);
void OutputRestartFile(char *restart_filname);
void OutputFixedHeader(FILE *out_fp);
void OutputDynamicsMasses(FILE *out_fp);
void OutputParams(char *tag, FILE *out_fp);

void MinimizeTypeToString(enum MinimizeType min_type, char *str);
void ThermostatTypeToString(enum ThermostatType thermo_type, char *str);
void ActionTypeToString(enum ActionType action_type, char *str);
void OutputTypeToString(enum OutputType output_type, char *str);
void BoolToString(int val, char *str);

void PrintAnalyticForces();
void PrintNumericalForces();
