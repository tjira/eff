#include "eff_global.h"
#include "eff_params.h"
#include "eff_access.h"
#include "eff_util.h"
#include "eff_input.h"

typedef struct
{
  char name[255];
  void (*fp)(int, char **);
} PARAM_FUNCTIONS;

const int max_param_functions = 100;
int num_param_functions;
PARAM_FUNCTIONS *param_functions;

void InitializeParameters()
{
  param_functions = (PARAM_FUNCTIONS *) malloc(sizeof(PARAM_FUNCTIONS) * max_param_functions);

  num_param_functions = 0;
  AddParameter("s_cutoff", ParamSCutoff);
  AddParameter("x_bound", ParamXYZBound);
  AddParameter("y_bound", ParamXYZBound);
  AddParameter("z_bound", ParamXYZBound);
  AddParameter("periodic", ParamPeriodic);

  AddParameter("ewald_log_precision", ParamEwald);
  AddParameter("ewald_re_cutoff", ParamEwald);
  AddParameter("ewald_autoset", ParamEwald);
  AddParameter("ewald_r_cutoff", ParamEwald);
  AddParameter("ewald_k_cutoff", ParamEwald);
  AddParameter("ewald_nuc_r", ParamEwald);
  AddParameter("ewald_max_re", ParamEwald);

  AddParameter("min", ParamMin);
  AddParameter("thermostat", ParamThermostat);
  AddParameter("start_temperature", ParamStartTemperature);
  AddParameter("dt", ParamDt);
  AddParameter("print_every", ParamPrintEvery);
  AddParameter("num_steps", ParamNumSteps);
  AddParameter("calc", ParamCalc);
  AddParameter("output_position", ParamOutput);
  AddParameter("output_velocity", ParamOutput);
  AddParameter("output_energy_forces", ParamOutput);
  AddParameter("output_restart", ParamOutput);
  AddParameter("electron_mass", ParamElectronMass);
  AddParameter("rand_seed", ParamRandSeed);

  AddParameter("e_field", ParamEField);
  AddParameter("e_field_freq", ParamEFieldFreq);
  AddParameter("e_field_duration", ParamEFieldPacketDuration);
  
  AddParameter("collapse_move", ParamCollapseMove);
}

void SetDefaultParameters()
{
  /* cutoff for overlap */
  params.s_cutoff = 1000.0;    

  /* boundaries */
  params.x_bound[0] = -INFINITY;
  params.x_bound[1] = INFINITY;

  params.y_bound[0] = -INFINITY;
  params.y_bound[1] = INFINITY;

  params.z_bound[0] = -INFINITY;
  params.z_bound[1] = INFINITY;
  
  params.periodic   = 0;

  /* ewald */
  params.ewald_log_precision = -6.0;
  params.ewald_re_cutoff     = 3.54;
  params.ewald_autoset       = 1;
  params.ewald_r_cutoff      = 7.0;
  params.ewald_k_cutoff      = 8.0;
  params.ewald_nuc_r         = 1e-10;
  params.ewald_max_re        = 4.5;

  /* minimization */
  params.min = CONJUGATE_GRADIENT;

  /* dynamics */
  params.thermostat = NONE;
  params.start_temperature = 0;
  params.dt = 0.005;

  /* updating */
  params.print_every = 50;
  params.num_steps = 10000;

  /* action */
  params.calc = SINGLE_PT;

  /* output */
  params.output_position = ALL;
  params.output_velocity = ALL;
  params.output_energy_forces = NONE2;
  params.output_restart = ALL;

  /* electron mass */
  params.electron_mass = 1.0;

  /* random numbers */
  params.rand_seed = 10000;

  /* external field */
  params.Ex = 0;
  params.Ey = 0;
  params.Ez = 0;
  params.Efreq = 0;
  params.Epacket_duration = 0;

  /* collpase */
  params.collapse_move = 0;
}

void AddParameter(char *name, void (*fp)(int, char **))
{
  strcpy(param_functions[num_param_functions].name, name);
  param_functions[num_param_functions].fp = fp;

  num_param_functions++;
  if (num_param_functions > max_param_functions) error("Added too many parameters.\n");
}

void ExecuteParameter(char *name, int argc, char **argv)
{
  /* Find function with name */
  int i;
  for (i = 0; i < num_param_functions; i++)
  {
    if (string_match(name, param_functions[i].name))
    {
      param_functions[i].fp(argc, argv);
      return;
    }
  }
  error("Parameter function %s not found.\n", name);
}

void ParamSCutoff(int argc, char **argv)
{
  if (argc != 3) error("Usage: s_cutoff = [overlap cutoff].\n");
  sscanf(argv[2], "%lf", &params.s_cutoff);
}

void ParamXYZBound(int argc, char **argv)
{
  if (argc != 4) error("Usage: x/y/z_bound = min max\n");

  double *left_val, *right_val;
  if (string_match(argv[0], "x_bound"))
  {
    left_val  = &params.x_bound[0];      right_val  = &params.x_bound[1];
  }
  else if (string_match(argv[0], "y_bound"))
  {
    left_val  = &params.y_bound[0];      right_val  = &params.y_bound[1];
  }
  else if (string_match(argv[0], "z_bound"))
  {
    left_val  = &params.z_bound[0];      right_val  = &params.z_bound[1];
  }

  ReadValue(argv[2], left_val);
  ReadValue(argv[3], right_val);
}

void ReadValue(char *str, double *val)
{
  if (string_match(str, "-infinity"))
    *val = -INFINITY;
  else if (string_match(str, "infinity"))
    *val = INFINITY;
  else
    sscanf(str, "%lf", val);
}

void ParamPeriodic(int argc, char **argv)
{
  if (argc != 3) error("Usage: periodic = [true/false]\n");
  if (string_match(argv[2], "true") || string_match(argv[2], "TRUE"))
    params.periodic = 1;
  else if (string_match(argv[2], "false") || string_match(argv[2], "FALSE"))
    params.periodic = 0;
}

void ParamEwald(int argc, char **argv)
{
  if (argc != 3) error("Usage: ewald_* = [value]\n");

  if (string_match(argv[0], "ewald_log_precision"))
    sscanf(argv[2], "%lf", &params.ewald_log_precision);
  else if (string_match(argv[0], "ewald_re_cutoff"))
    sscanf(argv[2], "%lf", &params.ewald_re_cutoff);
  else if (string_match(argv[0], "ewald_autoset"))
  {
    if (string_match(argv[2], "true") || string_match(argv[2], "TRUE"))
      params.ewald_autoset = 1;
    else if (string_match(argv[2], "false") || string_match(argv[2], "FALSE"))
      params.ewald_autoset = 0;
    else
      error("Parameter %s not recognized as an autoset parameter.\n", argv[2]);
  }
  else if (string_match(argv[0], "ewald_r_cutoff"))
  {
    sscanf(argv[2], "%lf", &params.ewald_r_cutoff);
    params.ewald_autoset = 0;
  }
  else if (string_match(argv[0], "ewald_k_cutoff"))
  {
    sscanf(argv[2], "%lf", &params.ewald_k_cutoff);
    params.ewald_autoset = 0;
  }
  else if (string_match(argv[0], "ewald_nuc_r"))
    sscanf(argv[2], "%lf", &params.ewald_nuc_r);
  else if (string_match(argv[0], "ewald_max_re"))
    sscanf(argv[2], "%lf", &params.ewald_max_re);
  else
    error("Paramter %s not recognized as an ewald parameter.\n", argv[0]);
}

void ParamMin(int argc, char **argv)
{
  if (argc != 3) error("Usage: min = [conjugate_gradient/newton]\n");
  if (string_match(argv[2], "conjugate_gradient"))
    params.min = CONJUGATE_GRADIENT;
  else if (string_match(argv[2], "newton"))
    params.min = NEWTON;
  else
    error("Unknown minimization option.\n");
}

void ParamThermostat(int argc, char **argv)
{
  if (argc != 3) error("Usage: thermostat = [none/nose_hoover]\n");
  if (string_match(argv[2], "none"))
    params.thermostat = NONE;
  else if (string_match(argv[2], "nose_hoover"))
    params.thermostat = NOSE_HOOVER;
  else
    error("Unknown thermostat option.\n");
}

void ParamStartTemperature(int argc, char **argv)
{
  if (argc != 3) error("Usage: start_temperature = [start_temperature]\n");
  sscanf(argv[2], "%lf", &params.start_temperature);
}

void ParamDt(int argc, char **argv)
{
  if (argc != 3) error("Usage: dt = [dt]\n");
  sscanf(argv[2], "%lf", &params.dt);
}

void ParamPrintEvery(int argc, char **argv)
{
  if (argc != 3) error("Usage: print_every = [print_every]\n");
  sscanf(argv[2], "%i", &params.print_every);
}

void ParamNumSteps(int argc, char **argv)
{
  if (argc != 3) error("Usage: num_steps = [num_steps]\n");
  sscanf(argv[2], "%i", &params.num_steps);
}

void ParamCalc(int argc, char **argv)
{
  if (argc != 3) error("Usage: calc = [single_pt/minimize/dynamics]\n");
  if (string_match(argv[2], "single_pt"))
    params.calc = SINGLE_PT;
  else if (string_match(argv[2], "minimize"))
    params.calc = MINIMIZE;
  else if (string_match(argv[2], "dynamics"))
    params.calc = DYNAMICS;
  else
    error("Unknown thermostat option.\n");
}

void ParamOutput(int argc, char **argv)
{
  if (argc != 3) error("Usage: output_* = [none/all/end]\n");

  enum OutputType *output_val;
  if (string_match(argv[0], "output_position"))
    output_val = &params.output_position;
  else if (string_match(argv[0], "output_velocity"))
    output_val = &params.output_velocity;
  else if (string_match(argv[0], "output_energy_forces"))
    output_val = &params.output_energy_forces;
  else if (string_match(argv[0], "output_restart"))
    output_val = &params.output_restart;
  else
    error("Output type not recognized.\n");

  if (string_match(argv[2], "none"))
    *output_val = NONE2;
  else if (string_match(argv[2], "all"))
    *output_val = ALL;
  else if (string_match(argv[2], "end"))
    *output_val = END;
  else
    error("Unknown output type option.\n");
}

void ParamElectronMass(int argc, char **argv)
{
  if (argc != 3) error("Usage: electron_mass = [electron_mass in amu]\n");
  sscanf(argv[2], "%lf", &params.electron_mass);
}

void ParamRandSeed(int argc, char **argv)
{
  if (argc != 3) error("Usage: rand_seed = [integer random number seed]\n");
  sscanf(argv[2], "%i", &params.rand_seed);
}

void ParamEField(int argc, char **argv)
{
  if (argc != 5) error("Usage: e_field = [Ex Ey Ez]\n");
  sscanf(argv[2], "%lf", &params.Ex);
  sscanf(argv[3], "%lf", &params.Ey);
  sscanf(argv[4], "%lf", &params.Ez);
}

void ParamEFieldFreq(int argc, char **argv)
{
  if (argc != 3) error("Usage: e_field_freq = [freq]\n");
  sscanf(argv[2], "%lf", &params.Efreq);
}

void ParamEFieldPacketDuration(int argc, char **argv)
{
  if (argc != 3) error("Usage: e_field_packet_duration = [duration]\n");
  sscanf(argv[2], "%lf", &params.Epacket_duration);
}

void ParamCollapseMove(int argc, char **argv)
{
  if (argc != 3) error("Usage: collapse_move = [length scale]\n");
  sscanf(argv[2], "%lf", &params.collapse_move);
}
