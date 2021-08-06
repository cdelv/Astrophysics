#include "header.h"


int main(int argc, char *argv[])
{
  std::vector <Cuerpo> star;
  std::vector <Cuerpo> star_dust;
  Config data;

  Configure(data);
  Initialize_system(data,star,star_dust);

  Propagate(data,star,star_dust);

  Results();
  
  return 0;
}
