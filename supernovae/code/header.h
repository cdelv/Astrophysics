#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

struct Config
{
  int x=1;
};

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

class Cuerpo{
private:
  vector3D  r, V, F;   double m, R;
public:
  Cuerpo(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0,double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double coeficiente);
  void Mueva_V(double dt, double coeficiente);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();}; 
  double Getz(void){return r.z();}; 
  double GetVx(void){return V.x();};
  double GetVy(void){return V.y();};
  double GetVz(void){return V.z();};

  friend class Colisionador;
};

  void Configure(Config &data);
  void Initialize_system(Config &data, std::vector <Cuerpo> &star, std::vector <Cuerpo> &star_dust);
  void Propagate(Config &data, std::vector <Cuerpo> &star, std::vector <Cuerpo> &star_dust);
void Results(void);


