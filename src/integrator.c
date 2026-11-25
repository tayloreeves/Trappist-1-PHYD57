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
void rk4_step(int n,double *x,double *v,const double *m,double dt,double *w)
{
    int d=nd*n;
    double *kv1=w,     *kv2=w+d,   *kv3=w+2*d, *kv4=w+3*d;
    //kx1 is just v
    double *kx2=w+4*d, *kx3=w+5*d, *kx4=w+6*d, *xt=w+7*d;  
    double h=dt, hh=0.5*dt, h6=dt/6.0;

    //set 1, slope at first point in step
    accel(n,x,m,kv1);

    //set 2, slope at midpoint in step
    for(int i=0;i<d;i++)
    {
        //velocity at midpoint, position at midpoint
        kx2[i]=v[i]+hh*kv1[i]; 
        xt[i]=x[i]+hh*v[i];
    }
    accel(n,xt,m,kv2);

    //set 3, slope at midpoint in step
    for(int i=0;i<d;i++)
    {
        //velocity at midpoint, position at midpoint
        kx3[i]=v[i]+hh*kv2[i];
        xt[i]=x[i]+hh*kx2[i];
    }
    accel(n,xt,m,kv3);

    //set 4, slope at endpoint in step
    for(int i=0;i<d;i++)
    {
        //velocity at endpoint, position at endpoint
        kx4[i]=v[i]+h*kv3[i];
        xt[i]=x[i]+h*kx3[i];
    }
    accel(n,xt,m,kv4);

    //weighted averages of the slopes to update position and velocity
    //x is updated first, then v is updated using the new x
    for(int i=0;i<d;i++)
    {
        x[i]+=h6*(v[i]+2.0*kx2[i]+2.0*kx3[i]+kx4[i]);
        v[i]+=h6*(kv1[i]+2.0*kv2[i]+2.0*kv3[i]+kv4[i]);
    }
}