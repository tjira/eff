#include <stdio.h>
#include "eff_global.h"
#include "eff_params.h"
#include "eff_input.h"
#include "eff_initial.h"
#include "eff_util.h"
#include "eff_access.h"
#include "eff_create.h"
#include "eff_initial.h"

int CountLines(FILE *fp, char *header)
{
  /* Count lines in an input file under a given header */
  char str[2048];
  if (fp == NULL) error("File not found.\n");

  rewind(fp);
  int in_section = 0, count = 0;
  while (fgets(str, 2048, fp))
  {
    StripNewline(str);
    if (strstr(str, "@") || strlen(str) == 0) in_section = 0;
    if (in_section) count++;
    if (string_match(str, header)) in_section = 1;
  }
  rewind(fp);

  return count;
}

void ProcessInputFile(FILE *fp)
{
	char str[2048];
  enum FileSection {NONE, PARAMS, NUC_VELOCITIES, ELEC_VELOCITIES, NUC_MASSES, ELEC_MASSES, NUCLEI, ELECTRONS};
	enum FileSection readmode = NONE;
	
	if (fp == NULL) error("File not found.\n");
	
	while(fgets(str, 2048, fp))
	{
    StripNewline(str);

		if (strstr(str, "@") || strlen(str) == 1 || strlen(str) == 0)
			readmode = NONE;
		
		if (readmode == NUCLEI)
      ProcessNucleus(str);
		else if (readmode == ELECTRONS)
      ProcessElectron(str);
    else if (readmode == PARAMS)
      ProcessParam(str);
    else if (readmode == NUC_VELOCITIES)
      ProcessNucVelocities(str);
    else if (readmode == ELEC_VELOCITIES)
      ProcessElecVelocities(str);
    else if (readmode == NUC_MASSES)
      ProcessNucMasses(str);
    else if (readmode == ELEC_MASSES)
      ProcessElecMasses(str);

    if (strstr(str, "@nuclei"))
			readmode = NUCLEI;
		else if (strstr(str, "@electrons"))
			readmode = ELECTRONS;
    else if (strstr(str, "@params"))
      readmode = PARAMS;
    else if (strstr(str, "@nuc_velocities"))
      readmode = NUC_VELOCITIES;
    else if (strstr(str, "@elec_velocities"))
      readmode = ELEC_VELOCITIES;
    else if (strstr(str, "@nuc_masses"))
      readmode = NUC_MASSES;
    else if (strstr(str, "@elec_masses"))
      readmode = ELEC_MASSES;
	}
}

void ProcessNucleus(char *str)
{
  /* nucleus: x y z charge */
  double x, y, z, charge;
  sscanf(str, "%lf %lf %lf %lf", &x, &y, &z, &charge);
	AddNucleus(x, y, z, charge);
}

void ProcessElectron(char *str)
{
	/* electron: x y z spin re */
  int spin;
  double x, y, z, re;
	sscanf(str, "%lf %lf %lf %i %lf", &x, &y, &z, &spin, &re);
	AddElectron(x, y, z, spin, re);
}

void ProcessNucVelocities(char *str)
{
  /* nuc velocity: # vx vy vz */
  int i;
  double vx, vy, vz;
  sscanf(str, "%i %lf %lf %lf", &i, &vx, &vy, &vz);
  AddNucV(i-1, vx, vy, vz);
}

void ProcessElecVelocities(char *str)
{
  /* elec velocity: # vx vy vz vr */
  int i;
  double vx, vy, vz, vr;
  sscanf(str, "%i %lf %lf %lf %lf", &i, &vx, &vy, &vz, &vr);
  AddElecV(i-1, vx, vy, vz, vr);
}

void ProcessNucMasses(char *str)
{
  /* nuc mass: # mx my mz */
  int i;
  double mx, my, mz;
  sscanf(str, "%i %lf %lf %lf", &i, &mx, &my, &mz);
  AddNucM(i-1, mx, my, mz);
}

void ProcessElecMasses(char *str)
{
  /* elec mass: # mx my mz mr */
  int i;
  double mx, my, mz, mr;
  sscanf(str, "%i %lf %lf %lf %lf", &i, &mx, &my, &mz, &mr);
  AddElecM(i-1, mx, my, mz, mr);
}

void ProcessParam(char *str)
{
  int argc;
  char *argv[128];

  /* parameter: param = value */
  ParseString(str, &argc, argv);
  if (argc < 3) error("Not enough parameters in: %s\n", str);
  if (!string_match(argv[1], "="))
    error("Expected equals in: %s\n", str);

  char *param_name = argv[0];
  ExecuteParameter(param_name, argc, argv);
}

int string_match(char *str1, char *str2)
{
  return !strcmp(str1, str2);
}

void StripNewline(char *str)
{
  int i = strlen(str);
  if (str[i - 1] == '\n') str[i - 1] = 0;
}

void ParseString(char *buffer, int *argc, char **argv)
{
  /* tokenizes input string and generates argc, argv */
  int nw = 0;
  char *token = strtok(buffer, " \t");
  while (token != NULL)
  {
    if (token[0] != '\n')
    {
      argv[nw] = token;
      nw++;
    }
    token = strtok(NULL, " \t"); 
  }
  *argc = nw;
}
