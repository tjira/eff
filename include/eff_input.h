#include <stdio.h>

int CountLines(FILE *fp, char *header);
void ProcessInputFile(FILE *fp);
void ProcessNucleus(char *str);
void ProcessElectron(char *str);
void ProcessNucVelocities(char *str);
void ProcessElecVelocities(char *str);
void ProcessNucMasses(char *str);
void ProcessElecMasses(char *str);
void ProcessParam(char *str);
int string_match(char *str1, char *str2);
void StripNewline(char *str);
void ParseString(char *buffer, int *argc, char **argv);
