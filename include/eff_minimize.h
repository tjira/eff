#include "eff_global.h"
#include <stdlib.h>

void AllocateMinimize(enum MinimizeType s_minmethod);
int MinStorageRequirements(enum MinimizeType s_minmethod);
void CalcContractFG(double *x, double *f, double *g);
enum MinimizeResult Minimize(enum MinimizeType s_minmethod, double eps, double acc, int max_numsteps, int print_every, void (*MinCallback)(int, int, double, double));
void GetMinimizeForces(double *forces);
void GetInitialPosition(double *positions);
void SetMinimizePositions(double *positions);

