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
  //COES estrella 1
  std::vector <double> coes;
  rv2coes(star,coes,0);
  fout.open("data/COES1.dat", std::fstream::in | std::fstream::out | std::fstream::app);
  fout << t <<","<< coes[0] <<","<< coes[1] <<","<< coes[2] <<","<< coes[3] <<","<< coes[4] <<","<< coes[5] <<std::endl;
  fout.close();
  //COES estrella 2
  fout.open("data/COES2.dat", std::fstream::in | std::fstream::out | std::fstream::app);
  rv2coes(star,coes,1);
  fout << t <<","<< coes[0] <<","<< coes[1] <<","<< coes[2] <<","<< coes[3] <<","<< coes[4] <<","<< coes[5] <<std::endl;
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
void rv2coes(std::vector <Cuerpo> &star, std::vector <double> &coes, int j)
{
  vector3D c, h, r, v, e, n;
  double i, raan, aop, ta, a, G=1, mu=G*(star[0].Getm()+star[1].Getm());
  double norm_r, norm2_v, norm_h, norm_e, norm_n;
  
  c=SMcenter(star);
  r.cargue(star[j].Getx(),star[j].Gety(),star[j].Getz()); r=r-c; //vector relativo al centro de masa
  v.cargue(star[j].GetVx(),star[j].GetVy(),star[j].GetVz());
  norm_r=norma(r);
  norm2_v=norma2(v);

  //angular momentum
  h=(r^v);
  norm_h=norma(h);

  //inclination
  i=std::acos(h.z()/norm_h);

  //eccentricity vector
  e=((norm2_v-mu/norm_r)*r-(r*v)*v)/mu;
  norm_e=norma(e);

  //node line
  n.cargue(-h.x(),h.y(),1);
  norm_n=norma(n);

  //RAAN
  raan=std::acos(n.x()/norm_n);
  if(n.y()<0)
    raan=2*M_PI-raan;

  //argument of perigee
  aop=std::acos(n*e/(norm_n*norm_e));
  if(e.z()<0)
    aop=2*M_PI-aop;

  //true anomaly
  ta=std::acos(e*r/(norm_r*norm_e));
  if(r*v<0)
    ta=2*M_PI-ta;

  //semi major axis
  a=norm_r*(1+norm_e*std::cos(ta))/(1-norm_e*norm_e);

  coes={a,norm_e,i,ta,aop,raan};
}
