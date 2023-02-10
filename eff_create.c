#include "eff_create.h"
#include "eff_global.h"
#include "eff_util.h"
#include <stdio.h>

int maxnuclei = 0;
int maxelectrons = 0;

void Initialize(int s_maxnuclei, int s_maxelectrons)
{
	maxnuclei = s_maxnuclei;
	maxelectrons = s_maxelectrons;

	nuc = (NUCLEUS *) malloc(sizeof(NUCLEUS) * maxnuclei);
	elec = (ELECTRON *) malloc(sizeof(ELECTRON) * maxelectrons);

	if (nuc == 0) error("Not enough memory for nuclei.");
	if (elec == 0) error("Not enough memory for electrons.");

	numnuclei = numelectrons = 0;
	return;
}

void AddNucleus(double x, double y, double z, double q)
{
	if (numnuclei > maxnuclei) error("More nuclei added than space allocated for them.");
	nuc[numnuclei].x = x;
	nuc[numnuclei].y = y;
	nuc[numnuclei].z = z;
	nuc[numnuclei].q = q;
	numnuclei++;
}

void AddElectron(double x, double y, double z, int spin, double re)
{
	if (numelectrons > maxelectrons) error("More electrons added than space allocated for them.");
	elec[numelectrons].x = x;
	elec[numelectrons].y = y;
	elec[numelectrons].z = z;
	elec[numelectrons].spin = spin;
	elec[numelectrons].r = re;
	numelectrons++;
}


