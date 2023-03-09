#ifndef GLOBAL_H
#define GLOBAL_H

typedef struct
{
	double x, y, z, q;
	double fx, fy, fz, energy;
} NUCLEUS;

typedef struct 
{
	double x, y, z, r;
	int spin;
	double fx, fy, fz, fr, energy;
} ELECTRON;

extern int numnuclei, numelectrons;
extern NUCLEUS *nuc;
extern ELECTRON *elec;

extern const int INFINITY;

/* global parameters */
enum BoundaryType {PERIODIC, REFLECT, ABSORB};
enum MinimizeType {CONJUGATE_GRADIENT, NEWTON};
enum ThermostatType {NONE, NOSE_HOOVER};
enum ActionType {SINGLE_PT, MINIMIZE, DYNAMICS};
enum OutputType {NONE2, ALL, END};
enum MinimizeResult {NORMAL, FAIL_EXCEED, FAIL_LINSEARCH, FAIL_NODESCENT};

typedef struct
{
  /* Cutoffs */
  double s_cutoff;   /* cutoff for Pauli overlap */

  /* Boundaries */
  double x_bound[2], y_bound[2], z_bound[2];  /* xyz boundaries (left and right) */
  int periodic;

  /* Ewald */
  double ewald_log_precision;
  double ewald_re_cutoff;
  int ewald_autoset;
  double ewald_r_cutoff, ewald_k_cutoff;
  double ewald_nuc_r, ewald_max_re;

  /* Minimization */
  enum MinimizeType min; 

  /* Dynamics */
  enum ThermostatType thermostat; 
  double start_temperature;
  double dt;

  /* Updating */
  int print_every;
  int num_steps;

  /* Action */
  enum ActionType calc;

  /* Output */
  enum OutputType output_position;
  enum OutputType output_velocity;
  enum OutputType output_energy_forces;
  enum OutputType output_restart;

  /* Masses */
  double electron_mass;

  /* Random numbers */
  unsigned int rand_seed;

  /* External field */
  double Ex, Ey, Ez;
  double Efreq;
  double Epacket_duration;

  /* Electron collapse */
  double collapse_move;

} PARAMS;
  
extern PARAMS params;

#endif
