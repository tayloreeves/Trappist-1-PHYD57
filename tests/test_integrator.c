#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../src/integrator.h"

//number of bodies in the test: sun + earth
#define nb 2

int main(void)
{
    const double pi=3.14159265358979323846;

    //masses in msun: sun, earth
    double m[nb]={1.0,3.003e-6};
    double mt=m[0]+m[1];

    //circular two-body orbit of separation r
    //relative speed vr = sqrt(g*mt/r), period p = 2*pi*sqrt(r^3/(g*mt))
    double r=1.0;
    double vr=sqrt(gg*mt/r);
    double p=2.0*pi*sqrt(r*r*r/(gg*mt));

    //start with the largest timestep, halve it each pass
    double dt0=p/100.0, err[3];

    printf("%-10s %14s %14s %8s\n","dt","pos error","dE/E","order");
    for(int c=0;c<3;c++)
    {
        //dt halves each pass: dt0, dt0/2, dt0/4
        double dt=dt0/(1<<c);
        long   ns=(long)round(p/dt);

        //state arrays and rk4 workspace (wsize(nb) doubles)
        double x[nd*nb]={0}, v[nd*nb]={0};
        double *w=malloc(wsize(nb)*sizeof(double));

        //barycentric start: centre of mass at rest at the origin
        //sun at -(m_earth/mt)*r, earth at +(m_sun/mt)*r, velocities split the same way
        x[0]=-(m[1]/mt)*r;   v[1]=-(m[1]/mt)*vr;
        x[nd]=(m[0]/mt)*r;   v[nd+1]=(m[0]/mt)*vr;
        double xe0=x[nd], ye0=x[nd+1];

        //reference energy before integrating
        double e0=energy(nb,x,v,m);

        //integrate exactly one period
        for(long s=0;s<ns;s++)
        {
            rk4_step(nb,x,v,m,dt,w);
        }

        //after one period the earth should be back where it started
        //distance from start = numerical error (the true orbit is closed)
        double dx=x[nd]-xe0, dy=x[nd+1]-ye0;
        err[c]=sqrt(dx*dx+dy*dy);

        //relative energy change over the run
        double de=fabs((energy(nb,x,v,m)-e0)/e0);

        //order = log2(error(dt)/error(dt/2)); rk4 should give about 4
        if(c==0)
        {
            printf("%-10.5f %14.3e %14.3e %8s\n",dt,err[c],de,"-");
        }
        else
        {
            printf("%-10.5f %14.3e %14.3e %8.2f\n",dt,err[c],de,log2(err[c-1]/err[c]));
        }
        free(w);
    }
    return 0;
}