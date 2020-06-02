#include "spinny.hpp"
#include "conics.hpp"

///////// Orbit Utility functions ///////////

// Forward Kepler Equation
inline double f2M(const double &f,const double &e) {
    const double tanE2 = sqrt((1.-e)/(1.+e))*tan(f/2.);
    const double E = atan(tanE2)*2;
    return E-e*sin(E);
}

// Solve Kepler's Equation using an accelerated Newton's method
// From "Fundamentals of Celestial Mechanics", Danby, 2nd ed., section 6.6
inline double kepler(const double &M,const double &e) {
    const double Ms = fmod(M,2.*M_PI);
    double x=Ms, es,ec, f,fp,dx;
    if (sin(Ms)>0) x += 0.85*e;
    else           x -= 0.85*e;
    for (unsigned nc(0);nc<10;nc++) {
        es=e*sin(x); ec=e*cos(x);
        f = x-es-Ms;
        if (fabs(f)<1e-15) break;
        fp = 1.0-ec; dx = -f/fp;
        dx = -f/(fp + 0.5*dx*es);
        dx = -f/(fp + 0.5*dx*es + dx*dx*ec/6.0);
        x += dx;
    }
    return 2.*atan(sqrt((1.+e)/(1.-e))*tan(x/2.));
}

///////// Simple Spin-Orbit Class /////////

class SimpleSpinOrbit { public:

    double t,a,e,h0, tol;
    double mu,n,P,M0,t0,coef;
    std::vector<double> arr0;
    Physical_Properties phys1,phys2;
    cashkarp_class<SimpleSpinOrbit> ck;

    SimpleSpinOrbit(double a0,double e0,double theta0,double dtheta0,double f0,
            Physical_Properties phys10,Physical_Properties phys20) {
        t = 0;
        a = a0;
        e = e0;
        arr0 = std::vector<double>(2);
        arr0[0] = dtheta0;
        arr0[1] = theta0;
        phys1 = phys10;
        phys2 = phys20;
        mu = phys1.mass+phys2.mass;
        n = sqrt(mu/(a*a*a));
        P = 2.*M_PI/n;
        M0 = f2M(f0,e);
        t0 = -M0/n;
        coef = 1.5*phys2.Ic2*phys1.mass;
        h0 = 1e-5*P;
        tol = 1e-15;
        ck = cashkarp_class<SimpleSpinOrbit>(this,&SimpleSpinOrbit::derivative,tol);
        ck.h0 = h0;
        ck.norm_arr = {1./n,1./(2.*M_PI)};
        //ck.verbose = true;
    }

    // Time derivative function
    void derivative(
            const std::vector<double> &arr,
            const double &time,
            std::vector<double> &deriv) const {
        const double M = n*(time-t0);
        const double f = kepler(M,e);
        const double r = a*(1.-(e*e))/(1.+e*cos(f));
        const double phi = f-arr[1];
        deriv[0] = (coef*sin(2.*phi))/(r*r*r);
        deriv[1] = arr[0];
    }
        
    // Integrate the system in time
    void evolve(double t1) {
        if(t1==t) return;
        ck.evolve(arr0,t,t1);
        arr0[1] = fmod(arr0[1],2.*M_PI);
        if (arr0[1]<0) arr0[1] += 2.*M_PI;
        t = t1;
    }
};
