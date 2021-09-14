#include "header.h"

void configure(char *argv[], Config &data)
{
  data.t_max=std::atof(argv[1]);
  data.dt=std::atof(argv[2]);
  data.N_frag=std::atoi(argv[3]);
  data.R_influ=std::atof(argv[4]);
  data.Frames=std::atoi(argv[5]);
  data.M_loss=std::atof(argv[6]);
  data.Exp_E=std::atof(argv[7]);
  data.TE=std::atof(argv[8]);
  data.G=std::atof(argv[9]);
}
