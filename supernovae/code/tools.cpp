#include "header.h"

void tokenize(std::string &str, char delim, std::vector<std::string> &out)
{
  out.clear();
  size_t start;
  size_t end = 0;
  
  while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
    {
      end = str.find(delim, start);
      out.push_back(str.substr(start, end - start));
    }
}
void init_files(void)
{
  //std::ofstream fout;
  //fout.open("data/animation/EnergiaTotal.dat");
  
}
void Animation(Config &data, std::vector <Cuerpo> &star, int i)
{
  std::ofstream fout;
  fout.open(("data/animation/Frame_star_"+std::to_string(i)+".csv").c_str());
  fout << "X" <<","<< "Y" <<","<< "Z" <<","<< "R" << std::endl;
  for(int i=0; i<2; i++)
    fout << star[i].Getx() <<","<< star[i].Gety()<<","<<star[i].Getz()<<","<<star[i].GetR() << std::endl;
  fout.close();

  std::ofstream fout1;
  fout1.open(("data/animation/Frame_fragment_"+std::to_string(i)+".csv").c_str());
  fout1 << "X" <<","<< "Y" <<","<< "Z" <<","<< "R" << std::endl;
  for(int i=2; i<star.size(); i++)
      fout1 << star[i].Getx() <<","<< star[i].Gety()<<","<<star[i].Getz()<<","<<star[i].GetR() << std::endl;
  fout.close();
  
  //fout.open("EnergiaTotal.dat", std::fstream::in | std::fstream::out | std::fstream::app);
}

void MultiplyVectorByScalar(std::vector<double> &myarray, double myconstant)
{
std::transform(myarray.begin(), myarray.end(), myarray.begin(), [&myconstant](auto& c){return c*myconstant;});
}
