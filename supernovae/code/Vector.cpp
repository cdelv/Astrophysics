// Vectores
#include "header.h"
using namespace std;

// Metodos de la clase vector3D
void vector3D::cargue(double x0, double y0, double z0){
  v[0]=x0; v[1]=y0; v[2]=z0;
}
void vector3D::show(void){
  cout << "(" <<v[0]<< "," <<v[1]<< "," <<v[2]<< ")" << endl;
}
vector3D vector3D::operator=(vector3D v2){
  for(int i=0;i<3;i++)
    v[i] = v2.v[i];
  return *this;
}
vector3D vector3D::operator+(vector3D v2){
  vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = v[i] + v2.v[i];
  return total;
}
vector3D vector3D::operator+=(vector3D v2){
  *this = *this + v2;
  return *this;
}
vector3D vector3D::operator*(double a){
  vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = a*v[i];
  return total;
}
vector3D vector3D::operator*=(double a){
  *this = (*this)*a;
  return *this;
}
vector3D vector3D::operator/(double a){
  double inver = 1.0/a;
  vector3D total;
  for(int i=0;i<3;i++)
    total.v[i] = inver*v[i];
  return total;
}
vector3D vector3D::operator-(vector3D v2){
  return *this + v2*(-1); 
}
vector3D vector3D::operator-=(vector3D v2){
  *this = *this - v2;
  return *this;
}
double vector3D::operator*(vector3D v2){
  double p=0;
  for(int i=0;i<3;i++)
    p += v[i]*v2.v[i];
  return p;
}
vector3D vector3D::operator^(vector3D v2){
  vector3D c;
  c.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
  c.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
  c.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];
  return c;
}
vector3D operator*(double a,vector3D v1){
  vector3D total;
  total = v1*a;	
  return total;
}
double norma2(vector3D v1){
  double n=0;
  for(int i=0;i<3;i++)
    n += v1.v[i]*v1.v[i];
  return n;
}
double norma(vector3D v1){
  return sqrt(norma2(v1));
}
