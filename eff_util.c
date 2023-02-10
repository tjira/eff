#include "eff_util.h"
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

void error(const char *format, ...)
{
	/* adapted from panic function in mpsim code */
  va_list ap;
	va_start(ap, format);
	vprintf(format, ap);
	va_end(ap);
	exit(1);
}

void rand_seed(unsigned int seed)
{
  srand(seed);
}

double rand_uni()
{
	return ((rand() % 1025) / 1024.0);
}

double rand_gauss()
{
	double r = 2.0, v1, v2;
	while (r > 1.0)
	{
		v1 = (rand() % 1025 - 512) / 512.0;
		v2 = (rand() % 1025 - 512) / 512.0;
		r = v1 * v1 + v2 * v2;
	}
	return v1 * sqrt(-2.0 * log(r) / r);
}

void lcase(char *str)
{
  while (*str != '\0')
  {
    *str = tolower(*str);
    str++;
  }
}

double GetMass(int Z)
{
	switch(Z)
        {
    case 0: return 1;
    case 1: return 1.008;
    case 2: return 4.003;
    case 3: return 6.941;
    case 4: return 9.012;
    case 5: return 10.81;
    case 6: return 12.01;
    case 7: return 14.01;
    case 8: return 16.00;
    case 9: return 18.998;
    case 10: return 20.180;
    default: return 1;

//error("Mass %i not found.\n", Z);
	}
}

char *GetType(int Z)
{
	switch(Z)
	{
case 0: return "Z";
case 1: return "H";
case 2: return "He";
case 3: return "Li";
case 4: return "Be";
case 5: return "B";
case 6: return "C";
case 7: return "N";
case 8: return "O";
case 9: return "F";
case 10: return "Ne";
    default: return "Z";

error("Nuclear charge %i not found.\n", Z);
	}
}

int GetZ(char *str)
{
  if (!strcmp("H", str))
    return 1;
  else if (!strcmp("He", str))
    return 2;
  else if (!strcmp("Li", str))
    return 3;
  else if (!strcmp("Be", str))
    return 4;
  else if (!strcmp("B", str))
    return 5;
  else if (!strcmp("C", str))
    return 6;
  else if (!strcmp("N", str))
    return 7;
  else if (!strcmp("O", str))
    return 8;
  else if (!strcmp("F", str))
    return 9;
  else if (!strcmp("Ne", str))
    return 10;
  else
    error("Element %s not found.\n", str);
	return -1;
}

void GetBasename(char *str, char *basename)
{
  /* find the last occurence of . in str */
  char *pos = strrchr(str, '.');
  if (pos == 0) error("File doesn't have an extension.");
  strncpy(basename, str, pos - str);
  basename[pos-str] = '\0';
}

void GetExtension(char *str, char *extension)
{
  /* find the last occrence of . in str */
  char *pos = strrchr(str, '.');
  strcpy(extension, pos + 1);
}

