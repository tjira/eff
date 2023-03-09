#include <math.h>
#include <stdio.h>
#include "eff_access.h"
#include "eff_constants.h"
#include "eff_efield.h"

double Ex0, Ey0, Ez0, Efreq, Epacket_duration;

void SetupExternalField(double s_Ex0, double s_Ey0, double s_Ez0, double s_Efreq, double s_Epacket_duration)
{
  Ex0 = s_Ex0; //* V_PER_CM;
  Ey0 = s_Ey0; //* V_PER_CM;
  Ez0 = s_Ez0; //* V_PER_CM;
  Efreq = s_Efreq * 2 * PI / T0_IN_FS / 1e15;

  /* Relates FWHM of intensity to standard deviation of electric field */
  Epacket_duration = s_Epacket_duration / T0_IN_FS / (2.0 * sqrt(log(2.0)));
}

void GetExternalField(double t, double *Ex, double *Ey, double *Ez)
{
  double A;
  if (Efreq == 0)
  {
    // Constant electric field 
    A = 1;
  }
  else
  {
    if (Epacket_duration == 0)
    {
      // Continuous wave 
      A = sin(Efreq * t);
    }
    else
    {
      // Wave packet -- Gaussian shaped pulse, where Epacket_duration is the standard deviation
      // Begin three standard deviations away.
      double t0 = 3 * Epacket_duration;
      A = sin(Efreq * t) * exp(-(t - t0) * (t - t0) / (2 * Epacket_duration * Epacket_duration));
    }
  }

  *Ex = A * Ex0;
  *Ey = A * Ey0;
  *Ez = A * Ez0;
}


void ApplyExternalField(double t)
{
  double Ex, Ey, Ez;
  GetExternalField(t, &Ex, &Ey, &Ez);

//  static int count = 0;
//  if (count % 100 == 0) printf("%f %f\n", t * T0_IN_FS, Ex);
//  count++;

  double q, x, y, z;
  int i;
  for (i = 0; i < NumNuclei(); i++)
  {
    // Force contribution
    q = GetNuclearCharge(i);
    AddNuclearForce(i, q * Ex, q * Ey, q * Ez);

    // Energy contribution
    GetNuclearPosition(i, &x, &y, &z);
    AddNucleusPE(i, -(q * x * Ex + q * y * Ey + q * z * Ez));
  }

  double r;
  for (i = 0; i < NumElectrons(); i++)
  {
    // Force contribution
    q = -1;
    AddElectronForce(i, q * Ex, q * Ey, q * Ez, 0);

    // Energy contribution
    GetElectronPosition(i, &x, &y, &z, &r);
    AddElectronPE(i, -(q * x * Ex + q * y * Ey + q * z * Ez));
  }
}

