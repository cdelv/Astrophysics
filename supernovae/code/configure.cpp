#include "header.h"

void configure(char *argv[], Config &data)
{
  data.t_max=std::atof(argv[1]);
  data.dt=std::atof(argv[2]);
  data.E=std::atof(argv[3]);
  data.N=std::atoi(argv[4]);
  data.R_influ=std::atof(argv[5]);
  data.Frames=std::atoi(argv[6]);
  data.N_frag=std::atoi(argv[7]);
  data.M_loss=std::atof(argv[8]);
  data.Exp_E=std::atof(argv[9]);
  init_files();
}
