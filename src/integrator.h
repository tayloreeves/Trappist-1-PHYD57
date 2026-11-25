#ifndef integrator_h
#define integrator_h

//units: au, yr, msun  ->  g = 4 pi^2
//state layout: x[nd*i+k] = coordinate k of body i, same for v and a
//nd = 2 for planar, compile with -Dnd=3 for 3d
//nd must be the same for integrator.c and every file that includes this header

#ifndef nd
#define nd 2
#endif

//g in au^3 / (msun yr^2) = 4 pi^2, defined here so integrator.c and the tests share one value
#define gg      39.47841760435743
#define wsize(n) (8*nd*(n))

//acceleration on every body, newtonian gravity, no softening
void accel(int n,const double *x,const double *m,double *a);

//total energy (kinetic + potential), for accuracy checks
double energy(int n,const double *x,const double *v,const double *m);

//one rk4 step of size dt, advances x and v in place
//w is caller-owned scratch of length wsize(n): one per thread, no allocation inside
void rk4_step(int n,double *x,double *v,const double *m,double dt,double *w);

#endif