#include <math.h>
#include "integrator.h"

//Newtonian gravity, RK4 integrator, acceleration of n bodies, energy of n bodies, RK4 step of n bodies
//n bodies, nd dimensions, x position, v velocity, m mass, a acceleration, w workspace, dt timestep
//a_i = sum_j G m_j (x_j-x_i)/|x_j-x_i|^3, E = sum_i 1/2 m_i v_i^2 - sum_{i<j} G m_i m_j / |x_j-x_i|
void accel(int n,const double *x,const double *m,double *a)
{
    //total number of coordinates and dimensions
    int d=nd*n;

    //initialize acceleration to zero
    for(int i=0;i<d;i++) 
    {
        a[i]=0.0;
    }

    for(int i=0;i<n;i++)
    {    
        //n(n-1)/2 pairs, newton 3rd law, no self-interaction
        for(int j=i+1;j<n;j++)
        {  
            double dr[nd],r2=0.0;

            //separation vector and distance squared
            for(int k=0;k<nd;k++)
            {
                dr[k]=x[nd*j+k]-x[nd*i+k]; 
                r2+=dr[k]*dr[k];
            }

            // force magnitude and acceleration update
            double f=gg/(r2*sqrt(r2));

            //newton's 3rd law: a_i = sum_j G m_j (x_j-x_i)/|x_j-x_i|^3, a_j = -G m_i (x_j-x_i)/|x_j-x_i|^3
            for(int k=0;k<nd;k++)
            {
                a[nd*i+k]+=f*m[j]*dr[k];
                a[nd*j+k]-=f*m[i]*dr[k];
            }
        }   
    }
}


//total energy of n bodies, this should be staying relatively constant, but will drift due to numerical error
//total energy = kinetic + potential
double energy(int n,const double *x,const double *v,const double *m)
{
    //initialize kintetic and potential to zero
    double ek=0.0,ep=0.0;
    for(int i=0;i<n;i++)
    {
        //kinteic energy of i is 1/2*m*v^2
        double v2=0.0;
        for(int k=0;k<nd;k++)
        {
            v2+=v[nd*i+k]*v[nd*i+k];
        }
        ek+=0.5*m[i]*v2;

        //potential energy of each pair of bodies, and sum over them
        for(int j=i+1;j<n;j++)
        {
            double r2=0.0;
            for(int k=0;k<nd;k++)
            {
                double q=x[nd*j+k]-x[nd*i+k];
                r2+=q*q;
            }
            ep-=gg*m[i]*m[j]/sqrt(r2);
        }
    }
    return ek+ep;
}


//RK4 step of n bodies, uses workspace w to store intermediate values
//w is of size 8*nd*n, and is used to store k1,k2,k3,k4 for both position and velocity, as well as intermediate positions
//k1 = v, k2 = v + dt/2 * a(x+dt/2*v), k3 = v + dt/2 * a(x+dt/2*k2), k4 = v + dt * a(x+dt*k3)
//x_{n+1} = x_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4), v_{n+1} = v_n + dt/6 * (a(x_n) + 2*a(x_n+dt/2*v) + 2*a(x_n+dt/2*k2) + a(x_n+dt*k3))
//uses old v for k1, so x is updated first
//one forest-ruth step of size dt, updates x and v in place
//4th order symplectic: alternates drifts (x += c*v*dt) and kicks (v += d*a*dt)
//3 force evaluations per step; w is scratch holding the acceleration (nd*n doubles)
//coefficients use th = 1/(2 - 2^(1/3)):
//  drifts: th/2, (1-th)/2, (1-th)/2, th/2     kicks: th, 1-2*th, th
//the middle kick coefficient is negative, which is normal for this scheme
void fr4_step(int n,double *x,double *v,const double *m,double dt,double *w)
{
    int d=nd*n;
    const double th=1.3512071919596578;

    //drift and kick step sizes
    double c1=0.5*th*dt, c2=0.5*(1.0-th)*dt;
    double d1=th*dt,     d2=(1.0-2.0*th)*dt;

    //drift 1, then kick 1 using the acceleration at the new positions
    for(int i=0;i<d;i++) x[i]+=c1*v[i];
    accel(n,x,m,w);
    for(int i=0;i<d;i++) v[i]+=d1*w[i];

    //drift 2, then kick 2 (negative step)
    for(int i=0;i<d;i++) x[i]+=c2*v[i];
    accel(n,x,m,w);
    for(int i=0;i<d;i++) v[i]+=d2*w[i];

    //drift 3, then kick 3
    for(int i=0;i<d;i++) x[i]+=c2*v[i];
    accel(n,x,m,w);
    for(int i=0;i<d;i++) v[i]+=d1*w[i];

    //final drift
    for(int i=0;i<d;i++) x[i]+=c1*v[i];
}