#include "eff_collapse.h"
#include "eff_dynamics.h"
#include "eff_access.h"
#include "eff_util.h"
#include "eff_global.h"
#include "eff_update.h"
#include <math.h>

void CollapseElectron(int i)
{
  /* Get position and velocity of electron */
  double x_old, y_old, z_old, r_old;
  GetElectronPosition(i, &x_old, &y_old, &z_old, &r_old);

  double vx_old, vy_old, vz_old, vr_old;
  GetElectronVelocity(i, &vx_old, &vy_old, &vz_old, &vr_old); 

  /* Find out new electron radius and size of position jump */
  double r_new, r_move;
  double mx, my, mz, mr;
  GetElectronMass(i, &mx, &my, &mz, &mr);
  r_new = GetNewR(r_old, vr_old, mr);
  if (r_new < 0) return;
  r_move = GetMoveR(r_old, r_new);

  /* If r_move is small, don't worry about collapsing electron */
  if (r_move < params.collapse_move) return;

  /* Otherwise, trial move the electron by a random Gaussian distance, where
     standard deviation = r_move / 2.
  */
  double x_new, y_new, z_new;
  double dx, dy, dz, r;
  dx = 0.5 * r_move * rand_gauss();
  dy = 0.5 * r_move * rand_gauss();
  dz = 0.5 * r_move * rand_gauss();
  dx = dy = dz = 0;
  x_new = x_old + dx;
  y_new = y_old + dy;
  z_new = z_old + dz;
 
  /* Calculate a revised r_new that conserves energy */
  double r_new2, r_move2;
  SetElectronPosition(i, x_new, y_new, z_new, r_new);
  r_new2 = SeekEnergy(GetTotalPE() + 0.5 * mr * vr_old * vr_old, i, r_new);
  SetElectronPosition(i, x_old, y_old, z_old, r_old);

  if (r_new2 < 0) return;
  r_move2 = GetMoveR(r_old, r_new2);
  if (r_move2 < params.collapse_move) return;

  /* Reject if the move does not fit within the new distribution */
  r = sqrt(dx * dx + dy * dy + dz * dz);
  //printf("r_move2 = %f, r_move = %f\n", r_move2, r_move);
  //printf("--- %f\n", exp(-(2 * r * r) * (1.0 / (r_move2 * r_move2) - 1.0 / (r_move * r_move))));
  if (rand_uni() > exp(-(2 * r * r) * (1.0 / (r_move2 * r_move2) - 1.0 / (r_move * r_move)))) return;

  /* Update the position and velocity of the electron.
     The translational velocity is unchanged, while the radial velocity is now zero.
  */
  SetElectronPosition(i, x_new, y_new, z_new, r_new2);
  SetElectronVelocity(i, vx_old, vy_old, vz_old, 0);
}

double GetNewR(double r_old, double vr_old, double m)
{
  /* Solution to the equation 
       (3/2) (1/r_old^2) + (1/2) m v_rold^2 = (3/2) (1/r_new^2);
   */

  double val = 1.0 / (r_old * r_old) + (m / 3.0) * vr_old * vr_old;
  if (val > 0)
    return 1.0 / sqrt(val);
  else
    return -1;
}

double GetMoveR(double r_old, double r_new)
{
  /* To conserve probability, we must have 
      r_old^2 = r_new^2 + r_move^2
  */
  if (r_new > r_old) return 0;
  return sqrt(r_old * r_old - r_new * r_new);
}

double SeekEnergy(double desiredE, int i, double r_guess)
{
  /* Tries to find radius of electron that will produce desired energy */
  /* Uses Newton's method */

  double r = r_guess;
  double f, df, delta_r;

  int iter = 1;
  do
  {
    SetElectronRadius(i, r);
    UpdateEnergyForces();
    f = GetTotalPE() - desiredE;
    df = -GetElectronRadiusForce(i);
    delta_r = -f / df;

    r += delta_r;
    if (r < 0 || r > 10) return -1;

    //printf("iterate: r = %f, f = %f, PE = %f, desired = %f\n", r, f, GetTotalPE(), desiredE);
    iter++;
    if (iter > 10) {return -1;}
  } while (fabs(delta_r) > 0.0001);
//  printf("converged.\n");

  return r;
}
