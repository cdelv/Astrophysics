#include "header.h"

const double G=1;

void Colisionador::CalculeFuerzas(std::vector<Cuerpo> &star){
  
  //borrar todas las fuerzas
  for(int i=0; i<star.size();i++)
    star[i].BorreFuerza();

  //Calcular las fuerzas entre las estrellas y los fragmentos
  //Los fragmentos no se hacen fuerza entre ellos (por ahora)

  for(int i=0; i<2; i++)
    for(int j=i+1; j<star.size(); j++)
      {
	CalculeFuerzaEntre(star[i], star[j]);
      }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo &Molecula1, Cuerpo &Molecula2){
  vector3D r21=Molecula2.r-Molecula1.r;
  double d=norma(r21);
  //double s=(Molecula1.R+Molecula2.R)-d; //distancia de inter penetracion   if (s>0){ colisionan  F2= K*std::pow(s,1.5)*n;
  vector3D n= r21*(1.0/d);
  vector3D F= -G*Molecula1.m*Molecula2.m*std::pow(d,-2)*n;

  Molecula2.AdicioneFuerza(F); Molecula1.AdicioneFuerza(-1*F); //3 ley de Newton
}
