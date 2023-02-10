#include <math.h>
#include <stdio.h>
#include "eff_access.h"
#include "eff_pressure.h"
#include "eff_constants.h"

// rigid virial assumes electrons stay same size under compression
// flexible virial assumes electrons compress under compression

double rigid_ke_virial, flexible_ke_virial;
double rigid_pe_virial, flexible_pe_virial;

void ClearVirial()
{
  rigid_ke_virial = flexible_ke_virial = 0;
  rigid_pe_virial = flexible_pe_virial = 0;
}

void AddForceVirial(double rx, double ry, double rz, double rr, double fx, double fy, double fz, double fr, enum ElectronType electron_type)
{
  if (electron_type == RIGID    || electron_type == BOTH) rigid_pe_virial += (rx * fx + ry * fy + rz * fz + rr * fr);
  if (electron_type == FLEXIBLE || electron_type == BOTH) flexible_pe_virial += (rx * fx + ry * fy + rz * fz + rr * fr);
}

void AddRadialForceVirial(double r, double f, enum ElectronType electron_type)
{
  if (electron_type == RIGID    || electron_type == BOTH) rigid_pe_virial += r * f;
  if (electron_type == FLEXIBLE || electron_type == BOTH) flexible_pe_virial += r * f;
}

void AddSizeForceVirial(double r, double f, enum ElectronType electron_type)
{
  if (electron_type == RIGID    || electron_type == BOTH) rigid_pe_virial += r * f;
  if (electron_type == FLEXIBLE || electron_type == BOTH) flexible_pe_virial += r * f;
}

void AddPotentialEnergyVirial(double energy, enum ElectronType electron_type)
{
  if (electron_type == RIGID    || electron_type == BOTH) rigid_pe_virial += energy;
  if (electron_type == FLEXIBLE || electron_type == BOTH) flexible_pe_virial += energy;
}

void AddKineticEnergyVirial(double energy, enum ElectronType electron_type)
{
  if (electron_type == RIGID    || electron_type == BOTH) rigid_ke_virial += 2 * energy;
  if (electron_type == FLEXIBLE || electron_type == BOTH) flexible_ke_virial += 2 * energy;
}

double GetRigidPEPressure(double volume)
{
  /* P = n kB T / V + (1 / (3 V)) * virial 
       = (3 n Kb T + virial) / 3V
       = (2 kinetic energy + virial) / 3V 
   */
  double pressure = rigid_pe_virial / (3.0 * volume);
  return P0_IN_PA * pressure;  // return pressure in Pascals.
}

double GetFlexiblePEPressure(double volume)
{
  /* P = n kB T / V + (1 / (3 V)) * virial 
       = (3 n Kb T + virial) / 3V
       = (2 kinetic energy + virial) / 3V 
   */
  double pressure = flexible_pe_virial / (3.0 * volume);
  return P0_IN_PA * pressure;  // return pressure in Pascals.
}

double GetRigidKEPressure(double volume)
{
  /* P = n kB T / V + (1 / (3 V)) * virial 
       = (3 n Kb T + virial) / 3V
       = (2 kinetic energy + virial) / 3V 
   */
  double pressure = rigid_ke_virial / (3.0 * volume);
  return P0_IN_PA * pressure;  // return pressure in Pascals.
}

double GetFlexibleKEPressure(double volume)
{
  /* P = n kB T / V + (1 / (3 V)) * virial 
       = (3 n Kb T + virial) / 3V
       = (2 kinetic energy + virial) / 3V 
   */
  double pressure = flexible_ke_virial / (3.0 * volume);
  return P0_IN_PA * pressure;  // return pressure in Pascals.
}


