#include "header.h"

const double const1 = 0.1786178958448091;      //zeta
const double const3 = -0.06626458266981849;   //Xi
const double const4 = -0.2123418310626054;   //Lambda
const double const2 = (1-2*const4)/2;       //(1-2*Lambda)/2 
const double const5 = 1-2*(const3+const1); //1-2*(Xi+zeta)

void Propagate(Config &data, std::vector <Cuerpo> &star)
{
  double t, tdibujo=0; int i=0;
  Animation(data, star, i); i++; //imprimir posicion inicial para la animacion
  
  for(t=0; t<=data.t_max; t+=data.dt, tdibujo+=data.dt)
    {
      if(tdibujo>data.t_max/data.Frames)
	{
	  Animation(data, star, i); i++;
	  Results(data,star,t);
	  Progress((1.0*i)/data.Frames);
	  tdibujo=0;
	}
      
      if( data.TE < t && t < data.TE+data.dt)
	{
	  explode(data,star);
	}
      
      integrate(data,star);
    }
}
void integrate(Config &data, std::vector <Cuerpo> &star)
{
  Colisionador Newton;
  //Muevase con Omelyan PEFRL
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const1);
  Newton.CalculeFuerzas(star,data);
  for(int i=0; i<star.size(); i++) star[i].Mueva_V(data.dt,const2);
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const3);
  Newton.CalculeFuerzas(star,data);
  for(int i=0; i<star.size(); i++) star[i].Mueva_V(data.dt,const4);
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const5);
  Newton.CalculeFuerzas(star,data);
  for(int i=0; i<star.size(); i++) star[i].Mueva_V(data.dt,const4);
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const3);
  Newton.CalculeFuerzas(star,data);
  for(int i=0; i<star.size(); i++) star[i].Mueva_V(data.dt,const2);
  for(int i=0; i<star.size(); i++) star[i].Mueva_r(data.dt,const1);
}
void explode(Config &data, std::vector <Cuerpo> &star)
{
  std::vector <double> masses;
  std::vector <double> energies;
  std::vector <double> x0;
  std::vector <double> y0;
  std::vector <double> z0;
  vector3D v0;
  vector3D r0;
  double normr0, normE;
  double rr,th,ph;
  
  std::mt19937 gen(4);//semilla 4

  //random distribution for mas and energy of fragments
  std::normal_distribution<> m{5,2};//mu,sigma  //d(gen)
  std::normal_distribution<> e{5,2};//mu,sigma  //d(gen)

  //random distribution for the initialization inside a sphere
  std::uniform_real_distribution<> ur(star[1].GetR()*1.5, star[1].GetR()*2.5); // uniform, unbiased //-std::cbrt(3*data.M_loss/(100*4*M_PI*rho))
  std::uniform_real_distribution<> ut(0, M_PI); // uniform, unbiased
  std::uniform_real_distribution<> up(0, 2*M_PI); // uniform, unbiased

  for(int i=0; i<data.N_frag; i++)
    {
      rr=ur(gen); th=ut(gen); ph=up(gen); //random point in a sphere
      masses.push_back(std::abs(m(gen)));
      energies.push_back(std::abs(e(gen)));
      x0.push_back(rr*std::sin(th)*std::cos(ph)+star[1].Getx());
      y0.push_back(rr*std::sin(th)*std::sin(ph)+star[1].Gety());
      z0.push_back(rr*std::cos(th)+star[1].Getz());
    }
  
  double sum1=0, sum2=0;
  for(auto i: masses)
    sum1+=i;
  for(auto i: energies)
    sum2+=i;

  MultiplyVectorByScalar(masses,1.0/sum1);
  MultiplyVectorByScalar(masses,star[1].Getm()*data.M_loss*1.0/100.0);

  MultiplyVectorByScalar(energies,1.0/sum2);
  MultiplyVectorByScalar(energies,data.Exp_E*1.0);
  
  for(int i=0; i<data.N_frag; i++)
    {
      r0.cargue(x0[i]-star[1].Getx(), y0[i]-star[1].Gety(), z0[i]-star[1].Getz());
      normr0=norma(r0);
      normE=std::sqrt(2*energies[i]/masses[i]);
      v0=r0*(normE*1.0/normr0);
      Cuerpo point(x0[i], y0[i], z0[i], v0.x(), v0.y(), v0.z(), masses[i], std::cbrt(3*masses[i]/(4*M_PI*star[1].Get_rho())));
      star.push_back(point);
      star[1].Add_m(-masses[i],v0);
    }
}
