#include "header.h"

void Results(Config &data, std::vector <Cuerpo> &star, double t)
{
  std::ofstream fout;
  
  //Energia del sistema
  fout.open("data/TotalEnergy.dat", std::fstream::in | std::fstream::out | std::fstream::app);
  double a=Kenergy(star);
  double b=Penergy(star,data);
  fout << t <<","<< a <<","<< b <<","<< a+b <<std::endl;
  fout.close();
  //Energia de las estrellas
  fout.open("data/StarEnergy.dat", std::fstream::in | std::fstream::out | std::fstream::app);
  a=SKenergy(star);
  b=SPenergy(star,data);
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
  //parametros orbitales
  std::vector <double> coes;
  rv2coes(star,coes,data.G,0);
  fout.open("data/COES.dat", std::fstream::in | std::fstream::out | std::fstream::app);
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
double Penergy(std::vector <Cuerpo> &star, Config &data)
{
  double sum=0;
  for(int i=0; i<star.size(); i++)
    for(int j=i+1; j<star.size(); j++)
      {
	sum-=data.G*star[i].Getm()*star[j].Getm()/(norma(star[i].Getr()-star[j].Getr()));
      }
  
  return sum;
}
double SKenergy(std::vector <Cuerpo> &star)
{
  return star[0].Getm()*star[0].GetV2()*0.5 + star[1].Getm()*star[1].GetV2()*0.5;
}
double SPenergy(std::vector <Cuerpo> &star,Config &data)
{
  double dist; 

  dist=norma(star[0].Getr()-star[1].Getr());
  
  return -data.G*star[0].Getm()*star[1].Getm()/dist;
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
void rv2coes(std::vector <Cuerpo> &star, std::vector <double> &coes, double G, int j)
{
  double mt=0;
  for(auto k: star)
    mt+=k.Getm();
  
  vector3D h, r, v, e, n, nn;
  double i, raan, aop, ta, a, mu=G*mt;
  double norm_r, norm2_v, norm_h, norm_e, norm_n;
  
  //c=Mcenter(star);
  r.cargue(star[1].Getx()-star[0].Getx(),star[1].Gety()-star[0].Gety(),star[1].Getz()-star[0].Getz()); //reduccion al problema de un cuerpo
  v.cargue(star[1].GetVx()-star[0].GetVx(),star[1].GetVy()-star[0].GetVy(),star[1].GetVz()-star[0].GetVz());
  norm_r=norma(r);
  norm2_v=norma2(v);

  //angular momentum
  h=(r^v);
  norm_h=norma(h);

  //eccentricity vector
  e=((v^h)*1.0/mu)-r*1.0/norm_r;
  //e=((norm2_v-mu/norm_r)*r-(r*v)*v)/mu;
  norm_e=norma(e);

  //node line vector
  nn.cargue(0,0,1);
  n=(nn^h);
  norm_n=norma(n);

  //true anomaly
  ta=std::acos(e*r/(norm_r*norm_e));
  if(r*v<0)
    ta=2*M_PI-ta;

  //inclination
  i=std::acos(h.z()/norm_h);

  //ascending node
  if(norm_n != 0)
    {
      raan=std::acos(n.x()/norm_n);
      if(n.y()<0)
	raan=2*M_PI-raan;
    }
  else
    raan = 0;
  
  //argument of perigee
  if(norm_n !=0)
    {
      aop=std::acos(n*e/(norm_n*norm_e));
      if(e.z()<0)
	aop=2*M_PI-aop;
    }
  else
    aop = 0;


  //semi major axis
  //a=norm_r*(1+norm_e*std::cos(ta))/(1-norm_e*norm_e);
  a=1/((2/norm_r)-norm2_v/mu);

  coes={a,norm_e,i,ta,aop,raan};
}

void CMVelocity(void)
{
  
  std::string str;
  char d=',';
  std::vector<std::string> a;

  double Vx, Vy, Vz, V;
  
  std::ofstream fout;
  fout.open("data/CMVelocity.dat");
  
  int j=0;
  std::ifstream CM("data/Baricenter.dat");
  while(std::getline(CM, str))
    {
      j++;
    }
  CM.close();
  
  double data[j][4]= {0.0};

  j=0;
  std::ifstream CMM("data/Baricenter.dat");
  while(std::getline(CMM, str))
    {
      tokenize(str,d,a);
      data[j][0]=std::atof((a[0]).c_str());
      data[j][1]=std::atof((a[1]).c_str());
      data[j][2]=std::atof((a[2]).c_str());
      data[j][3]=std::atof((a[3]).c_str());
      j++;
    }

  for(int i=0; i<j-1; i++)
    {
      Vx=(data[i+1][1]-data[i][1])/(std::abs(data[i+1][0]-data[i][0]));
      Vx=(data[i+1][2]-data[i][2])/(std::abs(data[i+1][0]-data[i][0]));
      Vx=(data[i+1][3]-data[i][3])/(std::abs(data[i+1][0]-data[i][0]));
      V=std::hypot(data[i+1][1]-data[i][1],data[i+1][2]-data[i][2],data[i+1][3]-data[i][3])/(std::abs(data[i+1][0]-data[i][0]));
      
      fout << data[i+1][0] << "," << Vx << "," << Vy << "," << Vz << "," << V << std::endl;
    }
}
void Energy_derivative(void)
{
  
  std::string str, stre;
  char d=',';
  std::vector<std::string> a;
  std::vector<std::string> b;

  double Es, E;
  
  std::ofstream fout;
  fout.open("data/Energy_derivative.dat");
  
  int j=0;
  std::ifstream En("data/TotalEnergy.dat");
  std::ifstream Ens("data/StarEnergy.dat");
  
  while(std::getline(En, str) and std::getline(Ens, stre))
    {
      j++;
    }
  En.close();
  Ens.close();
  
  double data[j][3]= {0.0};

  j=0;
  std::ifstream Enn("data/TotalEnergy.dat");
  std::ifstream Ensn("data/StarEnergy.dat");
  
  while(std::getline(Enn, str) and std::getline(Ensn, stre))
    {
      tokenize(str,d,a);
      tokenize(stre,d,b);
      data[j][0]=std::atof((a[0]).c_str());
      data[j][1]=std::atof((a[3]).c_str());
      data[j][2]=std::atof((b[3]).c_str());
      j++;
    }

  for(int i=0; i<j-1; i++)
    {
      E=(data[i+1][1]-data[i][1])/(std::abs(data[i+1][0]-data[i][0]));
      Es=(data[i+1][2]-data[i][2])/(std::abs(data[i+1][0]-data[i][0]));
      
      fout << data[i+1][0] << "," << E << "," << Es << std::endl;
    }
}
