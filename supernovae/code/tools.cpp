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

void print_body(Cuerpo &star)
{
  std::cout << std::fixed;
  std::cout << std::setprecision(2);
  std::cout << star.Getx() <<"\t"<< star.Gety() <<"\t"<< star.Getz() <<"\t"<< star.GetVx() <<"\t"<<  star.GetVy() <<"\t"<< star.GetVz() <<"\t"<< star.GetFx() <<"\t"<< star.GetFy() <<"\t"<< star.GetFz() <<"\t"<< star.Getm() <<"\t"<< star.GetR() <<"\t"<< star.Get_rho() <<"\t"<< star.GetF() <<std::endl;
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
void Progress(double progress)
{
  int barWidth = 60;

  if(progress<=1)
    {
      std::cout << "[";
      int pos = barWidth * progress;
      for (int i = 0; i < barWidth; ++i) {
	if (i < pos) std::cout << "=";
	else if (i == pos) std::cout << ">";
	else std::cout << " ";}
      
      std::cout << "] " << int(progress * 100.0) << " %\r";
      std::cout.flush();
    }
}

