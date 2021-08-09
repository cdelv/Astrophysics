#include "header.h"

const double const1 = 0.1786178958448091;      //zeta
const double const3 = -0.06626458266981849;   //Xi
const double const4 = -0.2123418310626054;   //Lambda
const double const2 = (1-2*const4)/2;       //(1-2*Lambda)/2 
const double const5 = 1-2*(const3+const1); //1-2*(Xi+zeta)

void Propagate(Config &data, std::vector <Cuerpo> &star)
{
  double t;
  std::cout << star[0].Getx() <<"\t"<< star[0].Gety() <<"\t"<< star[1].Getx() <<"\t"<<star[1].Gety()<< std::endl;
  for(t=0; t<data.t_max; t+=data.dt)
    {
      integrate(data,star);
      std::cout << star[0].Getx() <<"\t"<< star[0].Gety() <<"\t"<< star[1].Getx() <<"\t"<<star[1].Gety()<< std::endl;
    }
}
void integrate(Config &data, std::vector <Cuerpo> &star)
{
  Colisionador Newton;
  //Muevase con Omelyan PEFRL
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const1);
  Newton.CalculeFuerzas(star);
  for(int i=0; i<star.size(); i++) star[i].Mueva_V(data.dt,const2);
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const3);
  Newton.CalculeFuerzas(star);
  for(int i=0; i<star.size(); i++) star[i].Mueva_V(data.dt,const4);
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const5);
  Newton.CalculeFuerzas(star);
  for(int i=0; i<star.size(); i++) star[i].Mueva_V(data.dt,const4);
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const3);
  Newton.CalculeFuerzas(star);
  for(int i=0; i<star.size(); i++) star[i].Mueva_V(data.dt,const2);
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const1);
}
