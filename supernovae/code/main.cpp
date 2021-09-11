#include "header.h"

int main(int argc, char *argv[])
{
  std::vector <Cuerpo> star;
  Config data;

  configure(argv,data);
  Initialize_system(data,star);
 
  Propagate(data,star);

  Results();
  
  return 0;
}
