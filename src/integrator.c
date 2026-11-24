#include <math.h>
#include "integrator.h"

void accel(int n,const double *x,const double *m,double *a){
 int d=nd*n;
 for(int i=0;i<d;i++) a[i]=0.0;
 for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){   //n(n-1)/2 pairs, newton 3
  double dr[nd],r2=0.0;
  for(int k=0;k<nd;k++){dr[k]=x[nd*j+k]-x[nd*i+k]; r2+=dr[k]*dr[k];}
  double f=gg/(r2*sqrt(r2));
  for(int k=0;k<nd;k++){
   a[nd*i+k]+=f*m[j]*dr[k];
   a[nd*j+k]-=f*m[i]*dr[k];
  }
 }
}

double energy(int n,const double *x,const double *v,const double *m){
 double ek=0.0,ep=0.0;
 for(int i=0;i<n;i++){
  double v2=0.0;
  for(int k=0;k<nd;k++) v2+=v[nd*i+k]*v[nd*i+k];
  ek+=0.5*m[i]*v2;
  for(int j=i+1;j<n;j++){
   double r2=0.0;
   for(int k=0;k<nd;k++){double q=x[nd*j+k]-x[nd*i+k]; r2+=q*q;}
   ep-=gg*m[i]*m[j]/sqrt(r2);
  }
 }
 return ek+ep;
}

void rk4_step(int n,double *x,double *v,const double *m,double dt,double *w){
 int d=nd*n;
 double *kv1=w,     *kv2=w+d,   *kv3=w+2*d, *kv4=w+3*d;
 double *kx2=w+4*d, *kx3=w+5*d, *kx4=w+6*d, *xt=w+7*d;   //kx1 is just v
 double h=dt, hh=0.5*dt, h6=dt/6.0;

 accel(n,x,m,kv1);
 for(int i=0;i<d;i++){kx2[i]=v[i]+hh*kv1[i]; xt[i]=x[i]+hh*v[i];}
 accel(n,xt,m,kv2);
 for(int i=0;i<d;i++){kx3[i]=v[i]+hh*kv2[i]; xt[i]=x[i]+hh*kx2[i];}
 accel(n,xt,m,kv3);
 for(int i=0;i<d;i++){kx4[i]=v[i]+h*kv3[i];  xt[i]=x[i]+h*kx3[i];}
 accel(n,xt,m,kv4);
 for(int i=0;i<d;i++){
  x[i]+=h6*(v[i]+2.0*kx2[i]+2.0*kx3[i]+kx4[i]);   //uses old v, so x first
  v[i]+=h6*(kv1[i]+2.0*kv2[i]+2.0*kv3[i]+kv4[i]);
 }
}