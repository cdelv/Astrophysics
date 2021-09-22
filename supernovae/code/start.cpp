#include "header.h"

void Initialize_system(Config &data, std::vector <Cuerpo> &star)
{
  char d=' ';
  std::vector<std::string> a;
  std::ifstream Data("data/initial_conditions.txt");
  std::string str;
  std::getline(Data, str); //ignore 1 line

  while(std::getline(Data, str)) //get the position data
    {
      tokenize(str,d,a);
      Cuerpo point(std::atof((a[0]).c_str()), std::atof((a[1]).c_str()), std::atof((a[2]).c_str()), std::atof((a[3]).c_str()),
		   std::atof((a[4]).c_str()), std::atof((a[5]).c_str()), std::atof((a[6]).c_str()), std::atof((a[7]).c_str()));
      star.push_back(point);
    }
  init_files();
}
void init_files(void)
{
  std::ofstream fout;
  fout.open("data/TotalEnergy.dat");
  fout.close();

  fout.open("data/StarEnergy.dat");
  fout.close();

  fout.open("data/Mass.dat");
  fout.close();

  fout.open("data/Baricenter.dat");
  fout.close();

  fout.open("data/COES.dat");
  fout.close();
}
