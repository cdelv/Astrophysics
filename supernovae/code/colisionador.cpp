#include "header.h"


void Colisionador::CalculeFuerzas(std::vector<Cuerpo> &star, Config &data){
  vector3D r21, cm;
  double d, s, influ;
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
	cm=(star[0].r*star[0].m+star[1].r*star[1].m)/(star[0].m+star[1].m); //centro de masa de las estrellas
	influ=norma(cm-star[j].r); //distancia al centro de masa
	
	if (s>0) //colision entre fragmentos y estrellas
	  {
	    star[i].colide(star[j].m,star[j].V);
	    star.erase(star.begin() + j);
	    j--;
	    yes=1;
	  }
	
	if(influ>data.R_influ && norma(star[j].V)>std::sqrt(2*data.G*(star[0].m+star[1].m)/influ) && yes==0) //si el fragmento se sale de la esfera de influencia borrarlo
	  {
	    star.erase(star.begin() + j);
	    j--;
	    yes=2;
	  }
	
	if(yes==0) // si el fragmento colisiona o se escapa , no hay que calcular fuerzas
	  CalculeFuerzaEntre(star[i], star[j], data.G);
      }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo &Molecula1, Cuerpo &Molecula2, double G){
  vector3D r21=Molecula2.r-Molecula1.r;
  double d=norma(r21);  
  vector3D n= r21*(1.0/d);
  vector3D F= -G*Molecula1.m*Molecula2.m*std::pow(d,-2)*n;

  Molecula2.AdicioneFuerza(F); Molecula1.AdicioneFuerza(-1.0*F); //3 ley de Newton
}
