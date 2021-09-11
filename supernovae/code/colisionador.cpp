#include "header.h"

const double G=1;

void Colisionador::CalculeFuerzas(std::vector<Cuerpo> &star){
  vector3D r21;
  double d, s;
  //borrar todas las fuerzas
  for(int i=0; i<star.size();i++)
    star[i].BorreFuerza();

  //Calcular las fuerzas entre las estrellas y los fragmentos
  //Los fragmentos no se hacen fuerza entre ellos (por ahora)
  
  for(int i=0; i<2; i++)
    for(int j=i+1; j<star.size(); j++)
      {
	int yes=0;
	r21=star[j].r-star[i].r;
	d=norma(r21);
	s=(star[i].R+star[j].R)-d; //distancia de interpenetracion
	if (s>0)
	  {
	    star[i].colide(star[j].m,star[j].V);
	    star.erase(star.begin() + j);
	    j--;
	    yes=1;
	  }
	if(yes==0)
	  CalculeFuerzaEntre(star[i], star[j]);
      }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo &Molecula1, Cuerpo &Molecula2){
  vector3D r21=Molecula2.r-Molecula1.r;
  double d=norma(r21);  
  vector3D n= r21*(1.0/d);
  vector3D F= -G*Molecula1.m*Molecula2.m*std::pow(d,-2)*n;

  Molecula2.AdicioneFuerza(F); Molecula1.AdicioneFuerza(-1*F); //3 ley de Newton
}
