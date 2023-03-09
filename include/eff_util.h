#include <stdlib.h>

void error(const char *format, ...);
double rand_uni();
double rand_gauss();
void lcase(char *str);
double GetMass(int Z);
char *GetType(int Z);
int GetZ(char *str);
void GetBasename(char *str, char *basename);
void GetExtension(char *str, char *extension);
void rand_seed(unsigned int seed);
