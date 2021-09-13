#include "header.h"

void Results(Config &data, std::vector <Cuerpo> &star, double t)
{
  std::ofstream fout;
  //parametros orbitales

  //velocidad del centro de masa
  
  //Energia del sistema
  fout.open("data/TotalEnergy.dat", std::fstream::in | std::fstream::out | std::fstream::app);
  double a=Kenergy(star);
  double b=Penergy(star);
  fout << t <<","<< a <<","<< b <<","<< a+b <<std::endl;
  fout.close();
  //Energia de las estrellas
  fout.open("data/StarEnergy.dat", std::fstream::in | std::fstream::out | std::fstream::app);
  a=SKenergy(star);
  b=SPenergy(star);
  fout << t <<","<< a <<","<< b <<","<< a+b <<std::endl;
  fout.close();
  // masa de las estrellas
  fout.open("data/Mass.dat", std::fstream::in | std::fstream::out | std::fstream::app);
  fout << t <<","<< star[0].Getm()<<","<<star[1].Getm() <<std::endl;
  fout.close();
  //posicion del centro de masa
  vector3D c, d;
  c=Mcenter(star);
  d=SMcenter(star);
  fout.open("data/Baricenter.dat", std::fstream::in | std::fstream::out | std::fstream::app);
  fout << t <<","<< c.x() <<","<< c.y() <<","<< c.z() <<","<< d.x() <<","<< d.y() <<","<< d.z() <<std::endl;
  fout.close();  
}
double Kenergy(std::vector <Cuerpo> &star)
{
  double sum=0;
  for(auto i: star)
    sum+=i.Getm()*i.GetV2()*0.5;
  
  return sum;
}
double Penergy(std::vector <Cuerpo> &star)
{
  double sum=0; double G=1;
  for(int i=0; i<star.size(); i++)
    for(int j=i+1; j<star.size(); j++)
      {
	sum-=G*star[i].Getm()*star[j].Getm()/(norma(star[i].Getr()-star[j].Getr()));
      }
  
  return sum;
}
double SKenergy(std::vector <Cuerpo> &star)
{
  return star[0].Getm()*star[0].GetV2()*0.5 + star[1].Getm()*star[1].GetV2()*0.5;
}
double SPenergy(std::vector <Cuerpo> &star)
{
  double sum=0; double G=1, dist=0; 

  dist=norma(star[0].Getr()-star[1].Getr());
  
  return -G*star[0].Getm()*star[1].Getm()/dist;
}
vector3D Mcenter(std::vector <Cuerpo> &star)
{
  vector3D c;
  double mt=0;
  for(auto i: star)
    mt+=i.Getm();
  for(auto i: star)
    c+=i.Getm()*i.Getr();
  return c*1.0/mt;
}

vector3D SMcenter(std::vector <Cuerpo> &star)
{
  return (star[0].Getm()*star[0].Getr()+star[1].Getm()*star[1].Getr())/(star[0].Getm()+star[1].Getm());
}
