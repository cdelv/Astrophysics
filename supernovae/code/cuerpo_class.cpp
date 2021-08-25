#include "header.h"

Cuerpo::Cuerpo(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0,double R0)
{
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  m=m0;
  R=R0;
}
void Cuerpo::Mueva_r(double dt, double coeficiente){
  r+=V*dt*coeficiente;
}
void Cuerpo::Mueva_V(double dt, double coeficiente){
  V+=(F*dt*coeficiente)/m;
}
void Cuerpo::Add_m(double M){
  m+=M;
}
void Cuerpo::colide(double M, vector3D v)
{
  V=(m*V+M*v)/(m+M);
}
