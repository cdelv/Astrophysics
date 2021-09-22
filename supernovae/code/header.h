#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include<algorithm>
#include <iomanip>


struct Config
{
  double t_max;
  double dt;
  double R_influ;
  int Frames;
  int N_frag;
  double M_loss;
  double Exp_E;
  double TE;
  double G;
  bool sphere;
};

//Vector.cpp
class vector3D{
  double v[3];
 public:
  void   cargue(double x0, double y0, double z0);
  void   show(void);
  // Funciones de salida de componentes
  double x(void){return v[0];};
  double y(void){return v[1];};
  double z(void){return v[2];};
  //Lectura de Elementos
  double & operator[](int i){return v[i];};

  // Operaciones vectoriales
  vector3D    operator= (vector3D v2);
  vector3D    operator+ (vector3D v2);
  vector3D    operator+=(vector3D v2);
  vector3D    operator- (vector3D v2);
  vector3D    operator-=(vector3D v2);
  // Producto por escalar
  vector3D    operator* (double a);
  vector3D    operator*=(double a);
  friend  vector3D    operator* (double a,vector3D v1);	
  // Division por escalar
  vector3D    operator/ (double a);
  // Producto cruz
  vector3D    operator^ (vector3D v2);
  // Producto punto
  double operator* (vector3D v2);
  // Norma 
  friend  double norma2(vector3D v1);    
  friend  double norma(vector3D v1);    
};

//cuerpo_class.cpp
class Cuerpo{
private:
  vector3D  r, V, F;   double m, R, rho=1;
public:
  Cuerpo(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0,double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double coeficiente);
  void Mueva_V(double dt, double coeficiente);
  void collide(double dM, vector3D v);
  void expell(double dM, vector3D v);
  void Add_m(double M);
  void Edit_r(double r){R=r;};
  double Getx(void){return r.x();};
  double Gety(void){return r.y();}; 
  double Getz(void){return r.z();};
  vector3D Getr(void){return r;}; 
  double GetVx(void){return V.x();};
  double GetVy(void){return V.y();};
  double GetVz(void){return V.z();};
  double GetV2(void){return norma2(V);};
  vector3D GetV(void){return V;};
  double GetFx(void){return F.x();};
  double GetFy(void){return F.y();};
  double GetFz(void){return F.z();};
  double Getm(void){return m;};
  double GetR(void){return R;};
  double Get_rho(void){return rho;};
  double GetF(void){return norma(F);};

  friend class Colisionador;
};

class Colisionador{

private:

public:
  void CalculeFuerzas(std::vector<Cuerpo> &star, Config &data);
  void CalculeFuerzaEntre(Cuerpo &Molecula1, Cuerpo &Molecula2, double G);
};

//configure.cpp
void configure(char *argv[], Config &data);

//start.cpp
void Initialize_system(Config &data, std::vector <Cuerpo> &star);
void init_files(void);

//propagate.cpp
void Propagate(Config &data, std::vector <Cuerpo> &star);
void integrate(Config &data, std::vector <Cuerpo> &star);
void explode(Config &data, std::vector <Cuerpo> &star);

//results.cpp
void Results(Config &data, std::vector <Cuerpo> &star, double t);
double Kenergy(std::vector <Cuerpo> &star);
double Penergy(std::vector <Cuerpo> &star,Config &data);
double SKenergy(std::vector <Cuerpo> &star);
double SPenergy(std::vector <Cuerpo> &star,Config &data);
vector3D Mcenter(std::vector <Cuerpo> &star);
vector3D SMcenter(std::vector <Cuerpo> &star);
void rv2coes(std::vector <Cuerpo> &star, std::vector <double> &coes, double G, int j);
void CMVelocity(void);
void Energy_derivative(void);


//tools.cpp
void tokenize(std::string &str, char delim, std::vector<std::string> &out);
void Animation(Config &data, std::vector <Cuerpo> &star, int i);
void MultiplyVectorByScalar(std::vector<double> &v, double k);
void print_body(Cuerpo &star);
void Progress(double progress);


