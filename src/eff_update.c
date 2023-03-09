#include "eff_global.h"
#include "eff_erf.h"
#include "eff_bounds.h"
#include "eff_ewald.h"
#include "eff_pressure.h"
#include "eff_update.h"
#include "eff_cores.h"
#include <math.h>

int is_periodic;
double s_cutoff;

double rd1 = 5, rd2 = -5, rd3 = 5;

// multibody nuclear forces
double *mb_nuc_fx, *mb_nuc_fy, *mb_nuc_fz;

void IntitializeEnergyUpdater(double set_s_cutoff, int periodic)
{
  s_cutoff = set_s_cutoff;
  is_periodic = periodic;

  /* Initialize temporary arrays for multibody nuclear forces */
  mb_nuc_fx = (double *) malloc(sizeof(double) * numnuclei); 
  mb_nuc_fy = (double *) malloc(sizeof(double) * numnuclei); 
  mb_nuc_fz = (double *) malloc(sizeof(double) * numnuclei); 
}

void UpdateEnergyForces()
{
	int i;
	for (i = 0; i < numnuclei; i++)
	{
		nuc[i].fx = nuc[i].fy = nuc[i].fz = 0;
		nuc[i].energy = 0;
	}
	for (i = 0; i < numelectrons; i++)
	{
		elec[i].fx = elec[i].fy = elec[i].fz = elec[i].fr = 0;
		elec[i].energy = 0;
	}
  ClearVirial();
  
  UpdateKineticEnergy();
  if (is_periodic)
  {
    UpdatePauliPeriodic();
    UpdateEwaldEnergy();
  }
  else
    UpdateElecAndPauli();
}

void UpdateKineticEnergy()
{
	double re;
	int i, k;
        double E, fr;

        double k5 = 1.0;
        double k6 = 1.2;
        double rho, r_cut;
        double r_slope = 0;

        double tot_ke = 0;
	for (i = 0; i < numelectrons; i++)
	{
          re = elec[i].r;
          E = 1.5 / (re * re);
          fr = 3.0 / (re * re * re);

         if (!isCore(i))
          {
          double repel, total_repel = 1.0;
          double total_override = 1.0;
          double ei_fx[2], ei_fy[2], ei_fz[2], ei_fr[2];
          for (k = 0; k < 2; k++)
            ei_fx[k] = ei_fy[k] = ei_fz[k] = ei_fr[k] = 0;

          for (k = 0; k < numnuclei; k++)
          {
            mb_nuc_fx[k] = mb_nuc_fy[k] = mb_nuc_fz[k] = 0;

            double dx1, dy1, dz1;
            double dRadial1_dS1;
            double s1, radial1, re1;
            double dS1_dX1, dS1_dY1, dS1_dZ1, dS1_dR1;
            double dF_dX1, dF_dY1, dF_dZ1, dF_dR1;
            double r1;

            dx1 = elec[i].x - nuc[k].x;
            dy1 = elec[i].y - nuc[k].y;
            dz1 = elec[i].z - nuc[k].z;
            r1 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
            
            if (nuc[k].q < 3) 
            {
              // Nearby protons make orbitals more s-like.
              rho = 1.0; r_cut = 2.5; 
              s1 = r1 / (rho * re); 
              
              radial1 = quintic_spline2(s1, &dRadial1_dS1, 0.5, r_cut);
//radial1 = 1; dRadial1_dS1 = 0;
              
              SmallRForce(dx1, dy1, dz1, r1, 1.0 / (rho * re), &dS1_dX1, &dS1_dY1, &dS1_dZ1);
              dS1_dR1 = -r1 / (rho * re * re);

              total_override *= radial1;
              if (radial1 == 0.0) continue;
             
              double c4 = -dRadial1_dS1 / radial1;
 
              dF_dX1 = dS1_dX1 * c4;
              dF_dY1 = dS1_dY1 * c4;
              dF_dZ1 = dS1_dZ1 * c4;
              dF_dR1 = dS1_dR1 * c4;

              ei_fx[0] += dF_dX1;
              ei_fy[0] += dF_dY1;
              ei_fz[0] += dF_dZ1;
              ei_fr[0] += dF_dR1;

              mb_nuc_fx[k] = -dF_dX1;
              mb_nuc_fy[k] = -dF_dY1;
              mb_nuc_fz[k] = -dF_dZ1;
            }
            else
            {
              // Nearby cores make orbitals more p-like
              rho = 0.6; r_cut = 2.5;
              s1 = r1 / (rho * re);
              radial1 = quintic_spline(s1, &dRadial1_dS1, r_slope, r_cut, rd1, rd2, rd3);

//radial1 = 1; dRadial1_dS1 = 0;
         
              SmallRForce(dx1, dy1, dz1, r1, 1.0 / (rho * re), &dS1_dX1, &dS1_dY1, &dS1_dZ1);
              dS1_dR1 = -r1 / (rho * re * re);
            
              repel = 1.0 - radial1;
              total_repel *= repel;
              if (repel == 0.0) continue;
 
              double c4 = dRadial1_dS1 / repel;

              dF_dX1 = dS1_dX1 * c4;
              dF_dY1 = dS1_dY1 * c4;
              dF_dZ1 = dS1_dZ1 * c4;
              dF_dR1 = dS1_dR1 * c4;

              ei_fx[1] += dF_dX1;
              ei_fy[1] += dF_dY1;
              ei_fz[1] += dF_dZ1;
              ei_fr[1] += dF_dR1;

              mb_nuc_fx[k] = -dF_dX1;
              mb_nuc_fy[k] = -dF_dY1;
              mb_nuc_fz[k] = -dF_dZ1;
            }
         }
        
        // E = (k6 - k5) total_override * (1 - total_repel) + k5
         
        // We calculate derivatives using x dy + y dx = x y (dx/x + dy/y)
        double c5 = E * (k6 - k5) * total_override; 
        double c6 = c5 * (1 - total_repel);
        double c7 = c5 * -total_repel;

        elec[i].fx += c6 * ei_fx[0] + c7 * ei_fx[1];
        elec[i].fy += c6 * ei_fy[0] + c7 * ei_fy[1];
        elec[i].fz += c6 * ei_fz[0] + c7 * ei_fz[1];
        elec[i].fr += c6 * ei_fr[0] + c7 * ei_fr[1];
        
        for (k = 0; k < numnuclei; k++)
        {
          if (nuc[k].q == 1)
          {
            nuc[k].fx += mb_nuc_fx[k] * c6;
            nuc[k].fy += mb_nuc_fy[k] * c6;
            nuc[k].fz += mb_nuc_fz[k] * c6;
          }
          else
          {
            nuc[k].fx += mb_nuc_fx[k] * c7;
            nuc[k].fy += mb_nuc_fy[k] * c7;
            nuc[k].fz += mb_nuc_fz[k] * c7;
          }
        }

        // Modify original derivatives for repulsion as well

        //printf("total_override = %f, total_repel = %f\n", total_override, total_repel);
        double factor = (k6 - k5) * total_override * (1 - total_repel) + k5;
        E *= factor;
        fr *= factor;
        }
        
        AddSizeForceVirial(re, fr, FLEXIBLE); 

        elec[i].energy += E;
        elec[i].fr += fr;
tot_ke += E;

      }
}

void UpdateElecAndPauli()
{
	int i, j, k;
	double energy, frc, fre1, fre2;
	double q1, q2, r1, r2;
	double x1, y1, z1, dx, dy, dz;
	double q, re1, re2, rc;
        double fx, fy, fz;
	int spin1;

	/* Nuclear-nuclear interaction */
  double total_elec_energy = 0;
	for (i = 0; i < numnuclei; i++)
	{
		q1 = nuc[i].q; x1 = nuc[i].x; y1 = nuc[i].y; z1 = nuc[i].z;
		for (j = 0; j < i; j++)
		{
			dx = nuc[j].x - x1;
			dy = nuc[j].y - y1;
			dz = nuc[j].z - z1;
			rc = sqrt(dx * dx + dy * dy + dz * dz);

			energy = frc = fx = fy = fz = 0;
			q = q1 * nuc[j].q;

			ElecNucNuc(q, rc, &energy, &frc);

      RForce(dx, dy, dz, rc, frc, &fx, &fy, &fz);
      UpdateNucNuc(i, j, energy, fx, fy, fz);
      total_elec_energy += energy;
		}
	}

double tot_nuc_elec = 0;
	/* Nuclear-electron interaction */
	for (i = 0; i < numnuclei; i++)
	{
		q1 = nuc[i].q; x1 = nuc[i].x; y1 = nuc[i].y; z1 = nuc[i].z;
		for (j = 0; j < numelectrons; j++)
		{
			dx = elec[j].x - x1;
			dy = elec[j].y - y1;
			dz = elec[j].z - z1;
			rc = sqrt(dx * dx + dy * dy + dz * dz);
			re1 = elec[j].r;
			q = -q1; 

			energy = frc = fre1 = fx = fy = fz = 0;
      ElecNucElec(q, rc, re1, &energy, &frc, &fre1, i, j);

      SmallRForce(dx, dy, dz, rc, frc, &fx, &fy, &fz);
      UpdateNucElec(i, j, energy, fx, fy, fz, fre1);

      // remove size contribution to virial
      AddSizeForceVirial(re1, -fre1, RIGID);
      total_elec_energy += energy;
tot_nuc_elec += energy;
		}
	}

double tot_elec_elec = 0;
double tot_pauli = 0;

	/* Electron-electron repulsion */
	for (i = 0; i < numelectrons; i++)
	{
		re1 = elec[i].r; x1 = elec[i].x; y1 = elec[i].y; z1 = elec[i].z; spin1 = elec[i].spin;
		for (j = 0; j < i; j++)
		{
			dx = elec[j].x - x1;
			dy = elec[j].y - y1;
			dz = elec[j].z - z1;
			rc = sqrt(dx * dx + dy * dy + dz * dz);
			re2 = elec[j].r;

			energy = frc = fre1 = fre2 = 0;

      double pauli_energy, pauli_frc, pauli_fre1, pauli_fre2;

      ElecElecElec(rc, re1, re2, &energy, &frc, &fre1, &fre2, i, j);
      TotalPauli(i, j, rc, re1, re2, &pauli_energy, &pauli_frc, &pauli_fre1, &pauli_fre2);
      SmallRForce(dx, dy, dz, rc, frc + pauli_frc, &fx, &fy, &fz);
      UpdateElecElec(i, j, energy + pauli_energy, fx, fy, fz, fre1 + pauli_fre1, fre2 + pauli_fre2);

      total_elec_energy += energy;

      // remove size contribution from virial
      AddSizeForceVirial(re1, -fre1, RIGID);
      AddSizeForceVirial(re2, -fre2, RIGID);

      // Add virial contribution from Pauli
      AddRadialForceVirial(rc, pauli_frc, BOTH);
      AddSizeForceVirial(re1, pauli_fre1, FLEXIBLE);    
      AddSizeForceVirial(re2, pauli_fre2, FLEXIBLE);
		}
	}

//printf("%f %f %f\n", tot_nuc_elec, tot_elec_elec, tot_pauli);
//printf("%f %f\n", test_s_energy, test_s2_energy);

  /* Add virial contribution from electrostatics */
  AddPotentialEnergyVirial(total_elec_energy, BOTH);
}

void UpdatePauliPeriodic()
{
  int i, j;
  int spin1;
  double x1, y1, z1, dx, dy, dz, rc, rc2, re1, re2;

  for (i = 0; i < numelectrons; i++)
        {
                re1 = elec[i].r; x1 = elec[i].x; y1 = elec[i].y; z1 = elec[i].z; spin1 = elec[i].spin;
                for (j = 0; j < i; j++)
                {
      dx = BoundDX(elec[j].x - x1);
                        dy = BoundDY(elec[j].y - y1);
                        dz = BoundDZ(elec[j].z - z1);

			rc2 = dx * dx + dy * dy + dz * dz;
			if (rc2 < s_cutoff * s_cutoff)
			{
        		rc = sqrt(rc2);
            re2 = elec[j].r;

	      	double energy, frc, fre1, fre2;
		      energy = frc = fre1 = fre2 = 0;
  	   	 TotalPauli(i, j, rc, re1, re2, &energy, &frc, &fre1, &fre2);

			AddRadialForceVirial(rc, frc, BOTH);
			AddSizeForceVirial(re1, fre1, FLEXIBLE);    
			AddSizeForceVirial(re2, fre2, FLEXIBLE);

      	double fx, fy, fz;
	      SmallRForce(dx, dy, dz, rc, frc, &fx, &fy, &fz);
			UpdateElecElec(i, j, energy, fx, fy, fz, fre1, fre2);
			}
		}
	}
}

void TotalPauli(int i, int j, double rc, double re1, double re2, double *pauli_energy, double *pauli_frc, double *pauli_fre1, double *pauli_fre2)
{
  // a = e coalsence prevention
  // b = size stabilizer
  // c = general
  // d = opposite spin interaction

  // angle dependent 
      double s_energy, s_frc, s_fre1, s_fre2;
      s_energy = s_frc = s_fre1 = s_fre2 = 0;
      PauliElecElec(elec[j].spin == elec[i].spin, rc, re1, re2, &s_energy, &s_frc, &s_fre1, &s_fre2, 0.0, 0.0, 1.0, 0, i, j);

  // angle independent
      double s2_energy, s2_frc, s2_fre1, s2_fre2;
      s2_energy = s2_frc = s2_fre1 = s2_fre2 = 0;
      PauliElecElec(elec[j].spin == elec[i].spin, rc, re1, re2, &s2_energy, &s2_frc, &s2_fre1, &s2_fre2, 0.5, 3.0, 0, 0, i, j); 

      if (rc > s_cutoff) 
      {
		  *pauli_energy = 0;
		  *pauli_frc = 0;
		  *pauli_fre1 = 0;
		  *pauli_fre2 = 0;
			return;
		}

      double x = rc / s_cutoff;
		double spline = unit_cubic_spline(x); 
		double dspline = d_unit_cubic_spline(x) / s_cutoff;

		s_frc = s_frc * spline - s_energy * dspline;
		s_fre1 = spline * s_fre1;
		s_fre2 = spline * s_fre2;
		s_energy = spline * s_energy;

		s2_frc = s2_frc * spline - s2_energy * dspline;
		s2_fre1 = spline * s2_fre1;
		s2_fre2 = spline * s2_fre2;
		s2_energy = spline * s2_energy;

// Put in a very simple 3 body correction.

      if (!isCore(i) && !isCore(j))
      {
         // modifications to Pauli
         double k1 = 1.0; 
         double k2 = -1.0;
        // double k2 = -1.2;

         if (elec[i].spin != elec[j].spin)
         {
           k1 = k2 = 0; // no Pauli
         }

         double rho, r_cut;

         double dx1, dy1, dz1, dx2, dy2, dz2;
         double dot, r1_sq, r2_sq, r1, r2;
         double cos_theta, factor;
         double dCos_dX1, dCos_dY1, dCos_dZ1, dCos_dX2, dCos_dY2, dCos_dZ2;
         double dRadial1_dS1, dRadial2_dS2;
         double s1, s2, radial1, radial2, re1, re2;
         double dS1_dX1, dS1_dY1, dS1_dZ1, dS1_dR1, dS2_dX2, dS2_dY2, dS2_dZ2, dS2_dR2;
         double dF_dX1, dF_dY1, dF_dZ1, dF_dR1, dF_dX2, dF_dY2, dF_dZ2, dF_dR2;

         double repel, total_repel = 1.0;
         double ei_fx[2], ei_fy[2], ei_fz[2], ei_fr[2];
         double ej_fx[2], ej_fy[2], ej_fz[2], ej_fr[2];

         int k;
         for (k = 0; k < 2; k++)
         {
           ei_fx[k] = ei_fy[k] = ei_fz[k] = ei_fr[k] = 0;
           ej_fx[k] = ej_fy[k] = ej_fz[k] = ej_fr[k] = 0;
         }

         double override, total_override = 1.0;

         for (k = 0; k < numnuclei; k++)
         {
           mb_nuc_fx[k] = mb_nuc_fy[k] = mb_nuc_fz[k] = 0;

           dx1 = elec[i].x - nuc[k].x;
           dy1 = elec[i].y - nuc[k].y;
           dz1 = elec[i].z - nuc[k].z;
           dx2 = elec[j].x - nuc[k].x;
           dy2 = elec[j].y - nuc[k].y;
           dz2 = elec[j].z - nuc[k].z;

           dot   = dx1 * dx2 + dy1 * dy2 + dz1 * dz2; 
           r1_sq = dx1 * dx1 + dy1 * dy1 + dz1 * dz1; r1 = sqrt(r1_sq);
           r2_sq = dx2 * dx2 + dy2 * dy2 + dz2 * dz2; r2 = sqrt(r2_sq);
           cos_theta = dot / (r1 * r2);

           if (nuc[k].q < 3)
           {

             // Nearby protons make orbitals more s-like.
             rho = 1.0; r_cut = 2.5; 
             
             re1 = elec[i].r; re2 = elec[j].r;
             s1 = r1 / (rho * re1); s2 = r2 / (rho * re2);
             if (s1 > r_cut && s2 > r_cut) continue;

             radial1 = quintic_spline2(s1, &dRadial1_dS1, 0.5, r_cut);
             radial2 = quintic_spline2(s2, &dRadial2_dS2, 0.5, r_cut);

//radial1 = radial2 = 1; dRadial1_dS1 = dRadial2_dS2 = 0;
            
             SmallRForce(dx1, dy1, dz1, r1, 1.0 / (rho * re1), &dS1_dX1, &dS1_dY1, &dS1_dZ1);
             dS1_dR1 = -r1 / (rho * re1 * re1);
             SmallRForce(dx2, dy2, dz2, r2, 1.0 / (rho * re2), &dS2_dX2, &dS2_dY2, &dS2_dZ2);
             dS2_dR2 = -r2 / (rho * re2 * re2);

             double override = radial1 * radial2;

             // The final override is taken as the product of all these overrides
             total_override *= override;

             // Derivatives of the repulsion saved for later

             if (override == 0.0) continue;
             double c2 = -dRadial1_dS1 / radial1;
             double c3 = -dRadial2_dS2 / radial2;

             dF_dX1 = c2 * dS1_dX1;
             dF_dY1 = c2 * dS1_dY1;
             dF_dZ1 = c2 * dS1_dZ1;
             dF_dR1 = c2 * dS1_dR1;

             dF_dX2 = c3 * dS2_dX2;
             dF_dY2 = c3 * dS2_dY2;
             dF_dZ2 = c3 * dS2_dZ2;
             dF_dR2 = c3 * dS2_dR2;

             ei_fx[0] += dF_dX1;
             ei_fy[0] += dF_dY1;
             ei_fz[0] += dF_dZ1;
             ei_fr[0] += dF_dR1;
 
             ej_fx[0] += dF_dX2;
             ej_fy[0] += dF_dY2;
             ej_fz[0] += dF_dZ2;
             ej_fr[0] += dF_dR2;

             mb_nuc_fx[k] = -(dF_dX1 + dF_dX2);
             mb_nuc_fy[k] = -(dF_dY1 + dF_dY2);
             mb_nuc_fz[k] = -(dF_dZ1 + dF_dZ2);
           }
           else
           {
             // Nearby cores make orbitals more p-like.
             rho = 0.6; r_cut = 2.5; 
             double r_slope = 0;
             
             re1 = elec[i].r; re2 = elec[j].r;
             s1 = r1 / (rho * re1); s2 = r2 / (rho * re2);
             if (s1 > r_cut || s2 > r_cut) continue;

             radial1 = quintic_spline(s1, &dRadial1_dS1, r_slope, r_cut, rd1, rd2, rd3);
             radial2 = quintic_spline(s2, &dRadial2_dS2, r_slope, r_cut, rd1, rd2, rd3);
             //radial1 = cubic_spline(s1, &dRadial1_dS1, r_slope, r_cut);
             //radial2 = cubic_spline(s2, &dRadial2_dS2, r_slope, r_cut);

//radial1 = radial2 = 1; dRadial1_dS1 = dRadial2_dS2 = 0;

             // Derivatives pertaining to cosine angular term
             dCos_dX1 = (dx2 - dx1 * dot / r1_sq) / (r1 * r2);
             dCos_dY1 = (dy2 - dy1 * dot / r1_sq) / (r1 * r2);
             dCos_dZ1 = (dz2 - dz1 * dot / r1_sq) / (r1 * r2);
             dCos_dX2 = (dx1 - dx2 * dot / r2_sq) / (r1 * r2);
             dCos_dY2 = (dy1 - dy2 * dot / r2_sq) / (r1 * r2);
             dCos_dZ2 = (dz1 - dz2 * dot / r2_sq) / (r1 * r2);

             // Angular terms and derivatives 
             double angular = 1 - cos_theta * cos_theta;
             double dAngular = -2 * cos_theta;

             // Derivatives pertaining to radial terms
             SmallRForce(dx1, dy1, dz1, r1, 1.0 / (rho * re1), &dS1_dX1, &dS1_dY1, &dS1_dZ1);
             dS1_dR1 = -r1 / (rho * re1 * re1);
             SmallRForce(dx2, dy2, dz2, r2, 1.0 / (rho * re2), &dS2_dX2, &dS2_dY2, &dS2_dZ2);
             dS2_dR2 = -r2 / (rho * re2 * re2);

             double radial = radial1 * radial2;

             // The final repulsion is taken as the product of all these repulsions
             repel = 1.0 - radial * angular;
             total_repel *= repel;

             // Derivatives of the repulsion saved for later
             double c1 = dAngular * radial;
             double c2 = angular * radial2 * dRadial1_dS1;
             double c3 = angular * radial1 * dRadial2_dS2;

             if (repel == 0.0) continue;
             double invrepel = 1.0 / repel;
             
             dF_dX1 = (c1 * dCos_dX1 + c2 * dS1_dX1) * invrepel;
             dF_dY1 = (c1 * dCos_dY1 + c2 * dS1_dY1) * invrepel;
             dF_dZ1 = (c1 * dCos_dZ1 + c2 * dS1_dZ1) * invrepel;
             dF_dR1 = c2 * dS1_dR1 * invrepel;

             dF_dX2 = (c1 * dCos_dX2 + c3 * dS2_dX2) * invrepel;
             dF_dY2 = (c1 * dCos_dY2 + c3 * dS2_dY2) * invrepel;
             dF_dZ2 = (c1 * dCos_dZ2 + c3 * dS2_dZ2) * invrepel;
             dF_dR2 = c3 * dS2_dR2 * invrepel;

             ei_fx[1] += dF_dX1;
             ei_fy[1] += dF_dY1;
             ei_fz[1] += dF_dZ1;
             ei_fr[1] += dF_dR1;
 
             ej_fx[1] += dF_dX2;
             ej_fy[1] += dF_dY2;
             ej_fz[1] += dF_dZ2;
             ej_fr[1] += dF_dR2;

             mb_nuc_fx[k] = -(dF_dX1 + dF_dX2);
             mb_nuc_fy[k] = -(dF_dY1 + dF_dY2);
             mb_nuc_fz[k] = -(dF_dZ1 + dF_dZ2);
          }
        }

        // chi = total_override * (1 - total_repel)
        // E = ((k2 - k1) chi + k1) pauli 

        // We calculate derivatives using x dy + y dx = x y (dx/x + dy/y)
        double c5 = (s_energy * (k2 - k1)) * total_override; 
        double c6 = c5 * (1 - total_repel);
        double c7 = c5 * -total_repel;

        elec[i].fx += c6 * ei_fx[0] + c7 * ei_fx[1];
        elec[i].fy += c6 * ei_fy[0] + c7 * ei_fy[1];
        elec[i].fz += c6 * ei_fz[0] + c7 * ei_fz[1];
        elec[i].fr += c6 * ei_fr[0] + c7 * ei_fr[1];
        
        elec[j].fx += c6 * ej_fx[0] + c7 * ej_fx[1];
        elec[j].fy += c6 * ej_fy[0] + c7 * ej_fy[1];
        elec[j].fz += c6 * ej_fz[0] + c7 * ej_fz[1];
        elec[j].fr += c6 * ej_fr[0] + c7 * ej_fr[1];
        
        for (k = 0; k < numnuclei; k++)
        {
          if (nuc[k].q == 1)
          {
            nuc[k].fx += mb_nuc_fx[k] * c6;
            nuc[k].fy += mb_nuc_fy[k] * c6;
            nuc[k].fz += mb_nuc_fz[k] * c6;
          }
          else
          {
            nuc[k].fx += mb_nuc_fx[k] * c7;
            nuc[k].fy += mb_nuc_fy[k] * c7;
            nuc[k].fz += mb_nuc_fz[k] * c7;
          }
        }

        // Modify original derivatives for repulsion as well
        double chi = total_override * (1 - total_repel);
        factor = (k2 - k1) * chi + k1;

        s_energy *= factor;
        s_frc *= factor; s_fre1 *= factor; s_fre2 *= factor;
     }
  *pauli_energy = (s_energy + s2_energy);
  *pauli_frc = (s_frc + s2_frc);
  *pauli_fre1 = (s_fre1 + s2_fre1);
  *pauli_fre2 = (s_fre2 + s2_fre2);
}

void UpdateNucNuc(int i, int j, double energy, double fx, double fy, double fz)
{
	nuc[i].energy += 0.5 * energy;
	nuc[j].energy += 0.5 * energy;
	nuc[j].fx += fx;
	nuc[j].fy += fy;
	nuc[j].fz += fz;
	nuc[i].fx += -fx;
	nuc[i].fy += -fy;
	nuc[i].fz += -fz;
}

void UpdateNucElec(int i, int j, double energy, double fx, double fy, double fz, double fre1)
{
	nuc[i].energy += 0.5 * energy;
	elec[j].energy += 0.5 * energy;
	elec[j].fx += fx;
	elec[j].fy += fy;
	elec[j].fz += fz;
	elec[j].fr += fre1;
	nuc[i].fx += -fx;
	nuc[i].fy += -fy;
	nuc[i].fz += -fz;
}

void UpdateElecElec(int i, int j, double energy, double fx, double fy, double fz, double fre1, double fre2)
{
	elec[i].energy += 0.5 * energy;
	elec[j].energy += 0.5 * energy;
	elec[j].fx += fx;
	elec[j].fy += fy;
	elec[j].fz += fz;
	elec[j].fr += fre2;
	elec[i].fx += -fx;
	elec[i].fy += -fy;
	elec[i].fz += -fz;
	elec[i].fr += fre1;
}

void RForce(double dx, double dy, double dz, double rc, double force, double *fx, double *fy, double *fz)
{
  force /= rc;
  *fx = force * dx;
  *fy = force * dy;
  *fz = force * dz;
}

void SmallRForce(double dx, double dy, double dz, double rc, double force, double *fx, double *fy, double *fz)
{
	/* Handles case where rc is small to avoid division by zero */
	if (rc > 0.000001)
	{
		force /= rc;
		*fx = force * dx; *fy = force * dy; *fz = force * dz;
	}
	else
	{
		if (dx != 0) *fx = force / sqrt(1 + (dy * dy + dz * dz) / (dx * dx)); else *fx = 0.0;
		if (dy != 0) *fy = force / sqrt(1 + (dx * dx + dz * dz) / (dy * dy)); else *fy = 0.0;
		if (dz != 0) *fz = force / sqrt(1 + (dx * dx + dy * dy) / (dz * dz)); else *fz = 0.0;
		if (dx < 0) *fx = -*fx;
		if (dy < 0) *fy = -*fy;
		if (dz < 0) *fz = -*fz;
	}
}

void ElecNucNuc(double q, double rc, double *energy, double *frc)
{
	*energy += q / rc;
	*frc += q / (rc * rc);
}

void ElecNucNucDipole(double q1, double px1, double py1, double pz1, double q2, double px2, double py2, double pz2, double dx, double dy, double dz, double *energy, double *fx, double *fy, double *fz)
{
}

void SlaterPotential(double Z, double q, double rc, double *energy, double *frc)
{
  // E = (1 - (1 + Z r) exp(-2 Z r)) / r
  // Interaction of a Slater function with a point charge

  double f, expval, recip_rc;
  recip_rc = 1.0 / rc;
  expval = exp(-2 * Z * rc);
  f = 1 - (1 + Z * rc) * expval;

  *energy += q * f * recip_rc;
  *frc += q * recip_rc * (recip_rc * f - Z * (1 + 2 * Z * rc) * expval);
}

void ElecNucElec(double q, double rc, double re1, double *energy, double *frc, double *fre1, int i, int j)
{
	double a, b, arc;
        double coeff_a;

coeff_a = sqrt(2);
        
        // E = -Z/r Erf(a r / re) 
	// constants: sqrt(2), 2 / sqrt(pi)
	a = coeff_a / re1;
	arc = a * rc;

        // Interaction between nuclear point charge and Gaussian electron 
        double E, dEdr, dEdr1, f, df;

        f = erfoverx1(arc, &df);
        dEdr = -q * a * a * df;
        dEdr1 = q * (a / re1) * (f + arc * df);
	E = q * a * f;

        *energy += E;
        *frc += dEdr; 
        *fre1 += dEdr1;
}

void ElecElecElec(double rc, double re1, double re2, double *energy, double *frc, double *fre1, double *fre2, int i, int j)
{
	double a, b, arc, re, fre;
        double coeff_a;

coeff_a = sqrt(2);

	re = sqrt(re1 * re1 + re2 * re2);

	// constants: sqrt(2), 2 / sqrt(pi)
	a = coeff_a / re;
	arc = a * rc;
        
        // V_elecelec = E * F
        // where E = -Z/r Erf(a r / re) 
        //       F = (1 - (b s + c s^2) exp(-d s^2))
        // and s = r / re

        double E, dEdr, dEdr1, dEdr2, f, df;

	f = erfoverx1(arc, &df);
	dEdr = -a * a * df;
	fre = a * (f + arc * df) / (re * re);
	dEdr1 = fre * re1;
	dEdr2 = fre * re2;

	E = a * f;

        *energy += E;
        *frc += dEdr;
        *fre1 += dEdr1;
        *fre2 += dEdr2;
}

void PauliElecElec(int samespin, double rc, double re1, double re2, double *energy, double *frc, double *fre1, double *fre2, double a, double b, double c, double d, int i, int j)
{
	double ree, rem;
	double S, t1, t2, tt;
	double dSdr1, dSdr2, dSdr;
	double dTdr1, dTdr2, dTdr;
	double O, dOdS, dOdr1, dOdr2;

	/* Scaling factors and exchange parameter */
	ree = re1 * re1 + re2 * re2;
	rem = re1 * re1 - re2 * re2;

/*
        S = (2.82842712474619 / pow((re2 / re1 + re1 / re2), 1.5)) * exp(-rc * rc / ree);
	dSdr1 = S * ((-1.5 / re1) * (rem / ree) + 2 * re1 * rc * rc / (ree * ree));
	dSdr2 = S * ((1.5 / re2) * (rem / ree) + 2 * re2 * rc * rc / (ree * ree));
	dSdr  = S * (-2 * rc / ree);
*/
        double f, dfdr1, dfdr2;
        f = (2.82842712474619 / pow((re2 / re1 + re1 / re2), 1.5));
        dfdr1 = -(1.5 / re1) * (rem / ree) * f;
        dfdr2 = (1.5 / re2) * (rem / ree) * f; 

        double p, dpdr, dpdr1, dpdr2;
        p = sqrt(2 / ree) * rc;
        dpdr = sqrt(2 / ree); 
        dpdr1 = -(re1 / ree) * p; 
        dpdr2 = -(re2 / ree) * p; 
       
        double g, dgdp; 
        double alpha, beta;

        alpha = 0.5; beta = 0.0;
        //alpha = 0.5; beta = 0.0;

        g = exp(-alpha * p * p / (1 + beta * p));
        dgdp = -g * alpha * p * (2 + beta * p) / ((1 + beta * p) * (1 + beta * p));

        S = f * g;
        dSdr1 = dfdr1 * g + f * dgdp * dpdr1;
        dSdr2 = dfdr2 * g + f * dgdp * dpdr2;
        dSdr  = f * dgdp * dpdr;

	t1 =  1.5 * (1 / (re1 * re1) + 1 / (re2 * re2));
	t2 =  2.0 * (3 * ree - 2 * rc * rc) / (ree * ree);
	tt = t1 - t2;

	dTdr1 = -3 / (re1 * re1 * re1) - 12 * re1 / (ree * ree) + 8 * re1 * (-2 * rc * rc + 3 * ree) / (ree * ree * ree);
	dTdr2 = -3 / (re2 * re2 * re2) - 12 * re2 / (ree * ree) + 8 * re2 * (-2 * rc * rc + 3 * ree) / (ree * ree * ree);
	dTdr  = 8 * rc / (ree * ree);

  if (samespin == 1)
  {
    double n_exp = 2; 

    double S2 = S * S;
    double Sn = pow(S, n_exp);
    double one_minus_S2 = 1 - S2;
    double one_minus_Sn = 1 - Sn;

    double repel = a + b * (re1/re2 + re2/re1 - 2);
    double dRepel_dR1 = b * (1.0/re2 - re2 / (re1 * re1));
    double dRepel_dR2 = b * (1.0/re1 - re1 / (re2 * re2)); 

    O = 0.5 * (S2 / one_minus_S2) * (repel * S / (1 - Sn) + c * (1 - S));
    dOdS = -0.5 * repel * S2 * (-3 + S2 - (n_exp - 3) * Sn + (n_exp - 1) * Sn * S2) / (one_minus_S2 * one_minus_S2 * one_minus_Sn * one_minus_Sn) + 0.5 * c * S * (2 + S) / ((1 + S) * (1 + S));
    dOdr1 = 0.5 * dRepel_dR1 * S * S * S / (one_minus_S2 * one_minus_Sn);
    dOdr2 = 0.5 * dRepel_dR2 * S * S * S / (one_minus_S2 * one_minus_Sn);

    double P = d * S2;
    double dPdS = 2 * d * S;

    double dEdr = (dTdr * O  + tt * dOdS * dSdr);
    *fre1 -= dTdr1 * O + tt * (dOdS * dSdr1 + dOdr1);
    *fre2 -= dTdr2 * O + tt * (dOdS * dSdr2 + dOdr2);
    *frc  -= dEdr;
    *energy += tt * O;
  }
  else
  {
  }
}

double cubic_spline(double x, double *fx, double r_slope, double r_cut)
{
/* A cubic spline that starts at (0,0) with slope r_slope,
   continues to (1,1) with zero slope,
   and ends at (r_cut, 0) with zero slope.
*/

  if (x < 1.0)
  {
    *fx = (x - 1) * (-6 * x + r_slope * (3 * x - 1));
    return x * (r_slope * (x - 1) * (x - 1) + (3 - 2 * x) * x);
  }
  else if (x < r_cut)
  {
    *fx = 6 * (x - r_cut) * (x - 1) / ((r_cut - 1) * (r_cut - 1) * (r_cut - 1));
    return (x - r_cut) * (x - r_cut) * (-3 + r_cut + 2 * x) / ((r_cut - 1) * (r_cut - 1) * (r_cut - 1));
  }
  else
  {
    *fx = 0;
    return 0;
  }
}

double cubic_spline2(double x, double *fx, double r_cut)
{
  /* A cubic spline that starts at (0, 0) with zero slope and ends at (r_cut, 1) with zero slope. */
  
  double t = x / r_cut;
  if (t < 1)
  {
    *fx = -6 * t * (t - 1) / r_cut;
    return t * t * (3 - 2 * t);
  }
  else
  {
    *fx = 0;
    return 1;
  }
}

double quintic_spline(double x, double *fx, double r_slope, double r_cut, double rd_1, double rd_2, double rd_3)
{
/* A quintic spline that starts at (0,0) with slope r_slope,
   continues to (1,1) with zero slope,
   and ends at (r_cut, 0) with zero slope,
   and has second derivatives rd_1, rd_2, and rd_3 at those three points.
*/

  if (x < 1.0)
    return unit_quintic_spline(x, fx, 0, 0, 0, rd_1, 1, 1, 0, rd_2);
  else if (x < r_cut)
    return unit_quintic_spline(x, fx, 1, 1, 0, rd_2, r_cut, 0, 0, rd_3);
  else
    {*fx = 0; return 0;}
}
 
double quintic_spline2(double x, double *fx, double y0, double r_cut)
{
/* A quintic spline that starts at (0, 0) and ends at (rcut, 1) with
   zero derivative and zero second derivative, slope y0 at start.
*/

  if (x < r_cut)
     {return unit_quintic_spline(x, fx, 0, 0.5, 0.55, 0, r_cut, 1, 0, 0);}
  else
    {*fx = 0; return 1;}
}

double quintic_spline3(double x, double *fx, double scale_slope, double deriv2, double r_cut, double deriv22)
{
/* Quintic spline to modify potential functions */

//  *fx = 1; return x;

  if (x < r_cut)
    return unit_quintic_spline(x, fx, 0, 0, scale_slope, deriv2, r_cut, r_cut, 1, deriv22);
  else
    {*fx = 1; return x;}

/*
if (scale_slope == 1.0) {*fx = 1; return x;}

double frac = 0.4;

  if (x < frac)
    return unit_quintic_spline(x, fx, 0, 0, 1.385, 0, frac, frac*1.385, 1.1925, 0);
  else if (x < 1.0)
    return unit_quintic_spline(x, fx, frac, frac*1.385, 1.1925, 0, 1, 1, 1, 0);
  else
    {*fx = 1; return x;}
*/
}
    
double unit_quintic_spline(double x, double *fx, double x1, double y1, double dy1, double d2y1, double x2, double y2, double dy2, double d2y2)
{
  double t = (x - x1) / (x2 - x1);
  double t2, t3, t4, t5;
  t2 = t * t; t3 = t * t2; t4 = t2 * t2; t5 = t4 * t;
  dy1 *= (x2 - x1); dy2 *= (x2 - x1);
  d2y1 *= ((x2 - x1) * (x2 - x1)); d2y2 *= ((x2 - x1) * (x2 - x1));

  double h00, h10, h20, h01, h11, h21;
  h00 = -6 * t5 + 15 * t4 - 10 * t3 + 1;
  h10 = -3 * t5 + 8 * t4 - 6 * t3 + t;
  h20 = -0.5 * t5 + 1.5 * t4 - 1.5 * t3 + 0.5 * t2;
  h01 = 1 - h00;
  h11 = -3 * t5 + 7 * t4 - 4 * t3;
  h21 = 0.5 * t5 - t4 + 0.5 * t3;

  double dh00, dh10, dh20, dh01, dh11, dh21;
  dh00 = -30 * t2 + 60 * t3 - 30 * t4;
  dh10 = 1 - 18 * t2 + 32 * t3 - 15 * t4;
  dh20 = t - 4.5 * t2 + 6 * t3 - 2.5 * t4;
  dh01 = 30 * t2 - 60 * t3 + 30 * t4;
  dh11 = -12 * t2 + 28 * t3 - 15 * t4;
  dh21 = 1.5 * t2 - 4 * t3 + 2.5 * t4;  

  *fx = (dh00 * y1 + dh10 * dy1 + dh20 * d2y1 + dh01 * y2 + dh11 * dy2 + dh21 * d2y2) / (x2 - x1);
  return h00 * y1 + h10 * dy1 + h20 * d2y1 + h01 * y2 + h11 * dy2 + h21 * d2y2;
}


double unit_cubic_spline(double x)
{
  return x * x * (2 * x - 3) + 1;
}

double d_unit_cubic_spline(double x)
{
  return 6 * x * x - 6 * x;
}

