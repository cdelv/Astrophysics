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
  V+=F*dt*coeficiente/m;
}
void Cuerpo::Add_m(double dM){
  m+=dM;
}
void Cuerpo::collide(double dM, vector3D v)
{
  Add_m(dM);
  V=(m*V+dM*v)/m;
}
void Cuerpo::expell(double dM, vector3D v)
{
  V=(m*V-dM*v)/(m-dM);
  Add_m(-dM);
}
