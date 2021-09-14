#include "header.h"

Cuerpo::Cuerpo(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0,double R0)
{
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  F.cargue(0,0,0);
  m=m0;
  R=R0;
  rho=3*m/(4*M_PI*std::pow(R,3));
}
void Cuerpo::Mueva_r(double dt, double coeficiente){
  r+=V*dt*coeficiente;
}
void Cuerpo::Mueva_V(double dt, double coeficiente){
  V+=(F+dmdt)*dt*coeficiente/m;
  dmdt.cargue(0,0,0);
}
void Cuerpo::Add_m(double M, vector3D v){
  m+=M;
  dmdt+=M*v;
  Edit_r(std::cbrt(3*m/(4*M_PI*rho)));
}
void Cuerpo::colide(double M, vector3D v)
{
  //V=(m*V+M*v)/(m+M);
  //v.cargue(0,0,0);
  Add_m(M,v);
}
