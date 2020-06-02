#ifndef __CONICS_HPP__
#define __CONICS_HPP__

#include <cmath>
#include <vector>
#include <iostream>

// Namespace for neatness
namespace conics_util {
    inline double norm(const std::vector<double> &v);
    inline double dot(const std::vector<double> &v1,const std::vector<double> &v2);
    inline std::vector<double> cross(const std::vector<double> &u,const std::vector<double> &v);
    inline double sgn(const double &a);
    inline double init_elip(const double &dt,const double &mu,const double &alpha,const double &r0,const double &u);
    inline double init_hyper(const double &dt,const double &mu,const double &alpha,const double &r0,const double &u);
    void stumpff(double x,double &c0,double &c1,double &c2,double &c3);
    void kepler(const double &r0,const double &dt,const double &mu,const double &u,const double &alpha,
                    double &s,double &c0,double &c1,double &c2,double &c3,double &fp);
    double vsep(const std::vector<double> &v1,const std::vector<double> &v2);
}

// Utility functions

inline double conics_util::norm(const std::vector<double> &v) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline double conics_util::dot(const std::vector<double> &v1,const std::vector<double> &v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

inline std::vector<double> conics_util::cross(const std::vector<double> &u,const std::vector<double> &v) {
    std::vector<double> s(3);
    s[0] = u[1]*v[2] - u[2]*v[1];
    s[1] = u[2]*v[0] - u[0]*v[2];
    s[2] = u[0]*v[1] - u[1]*v[0];
    return s;
}

inline double conics_util::sgn(const double &a) {
    if(a<0) return -1.;
    return 1.;
}

// Danby pages 178-180
// Translated from BASIC, because it's not 1988 any more

// Smart initial Kepler value for elliptical motion
inline double conics_util::init_elip(const double &dt,const double &mu,const double &alpha,
        const double &r0,const double &u) {
    const double a = mu/alpha;
    const double en = sqrt(mu/(a*a*a));
    const double ec = 1. - (r0/a);
    const double es = u/(en*a*a);
    const double e = sqrt(ec*ec + es*es);
    const double dt2 = dt - (floor(en*dt/(2.*M_PI))*(2.*M_PI/en));
    const double y = en*dt2 - es;
    const double sigma = conics_util::sgn(es*cos(y) + ec*sin(y));
    const double x = y + sigma*.85*e;
    return x/sqrt(alpha);
}

// Smart initial Kepler value for hyperbolic motion
inline double conics_util::init_hyper(const double &dt,const double &mu,const double &alpha,
        const double &r0,const double &u) {
    const double a = mu/alpha;
    const double en = sqrt(-mu/(a*a*a));
    const double ch = 1. - (r0/a);
    const double sh = u/sqrt(-a*mu);
    const double e = sqrt(ch*ch - sh*sh);
    const double dm = en*dt;
    if(dm<0)
        return -log((-2.*dm + 1.8*e)/(ch-sh))/sqrt(-alpha);
    else
        return  log(( 2.*dm + 1.8*e)/(ch+sh))/sqrt(-alpha);
}

// Stumpff functions for Kepler estimation
void conics_util::stumpff(double x,double &c0,double &c1,double &c2,double &c3) {
    unsigned n = 0;
    const static double xm = .1;
    while(fabs(x)>=xm) {
        n += 1;
        x *= .25;
    }
    c2 = (1.-x*(1.-x*(1.-x*(1.-x*(1.-x*(1.-x/182.)/132.)/ 90.)/56.)/30.)/12.)/2.;
    c3 = (1.-x*(1.-x*(1.-x*(1.-x*(1.-x*(1.-x/210.)/156.)/110.)/72.)/42.)/20.)/6.;
    c1 = 1. - (x*c3);
    c0 = 1. - (x*c2);
    while(n>0) {
        n -= 1;
        c3 = (c2 + c0*c3)*.25;
        c2 = c1*c1*.5;
        c1 = c0*c1;
        c0 = 2.*(c0*c0) - 1.;
    }
    return;
}

// Universal Kepler equation solver
void conics_util::kepler(const double &r0,const double &dt,const double &mu,const double &u,const double &alpha,
        double &s,double &c0,double &c1,double &c2,double &c3,double &fp) {

    const static double tol = 1e-13;

    // Simple initial value for small delta
    if(fabs(dt/r0)<0.2)
        s = (dt/r0) - (dt*dt*u)/(2.*r0*r0*r0);

    // Smart initial value for elliptical motion
    else if(alpha>0)
        s = conics_util::init_elip(dt,mu,alpha,r0,u);

    // Smart initial value for hyperbolic motion
    else
        s = conics_util::init_hyper(dt,mu,alpha,r0,u);

    // Save initial value
    double st = s, x,f,fpp,fppp,ds;

    // Catch parabolas being dumb
    if(s!=s) s = 1;

    // Start with accelerated Newton-Rapson
    for(unsigned nc=0;nc<7;nc++) {
        x = s*s*alpha;
        conics_util::stumpff(x,c0,c1,c2,c3);
        c1 *= s;
        c2 *= s*s;
        c3 *= s*s*s;
        f = r0*c1 + u*c2 + mu*c3 - dt;
        fp = r0*c0 + u*c1 + mu*c2;
        fpp = (-r0*alpha + mu)*c1 + u*c0;
        fppp = (-r0*alpha + mu)*c0 - u*alpha*c1;
        ds = -f/fp;
        ds = -f/(fp + ds*fpp*.5);
        ds = -f/(fp + ds*fpp*.5 + (ds*ds*fppp/6.));
        s += ds;
        if(fabs(ds)<tol) return;
    }

    // If Newton failed, try Laguarre-Conway
    s = st;
    const static double ln = 5.;
    for(unsigned nc=0;nc<100;nc++) {
        x = s*s*alpha;
        conics_util::stumpff(x,c0,c1,c2,c3);
        c1 *= s;
        c2 *= s*s;
        c3 *= s*s*s;
        f = r0*c1 + u*c2 + mu*c3 - dt;
        fp = r0*c0 + u*c1 + mu*c2;
        fpp = (-r0*alpha + mu)*c1 + u*c0;
        ds = -ln*f/(fp + conics_util::sgn(fp)*sqrt(fabs((ln-1.)*(ln-1.)*fp*fp - (ln-1.)*ln*f*fpp)));
        s += ds;
        if(fabs(ds)<tol) return;
    }

    // If that failed, warn, but use what we have
    std::cout<<"# Kepler: WARNING! No convergence! ds = "<<ds<<'\n';
}

// Universal Kepler motion solver
std::vector<double> danbyuni(const std::vector<double> &state0,const double &mu,const double &dt) {

    // Create initial parameters
    const double r0 = sqrt(state0[0]*state0[0] + state0[1]*state0[1] + state0[2]*state0[2]);
    const double v0s = state0[3]*state0[3] + state0[4]*state0[4] + state0[5]*state0[5];
    const double u = state0[0]*state0[3] + state0[1]*state0[4] + state0[2]*state0[5];
    const double alpha = (2.*mu/r0) - v0s;

    // Solve conic motion
    double s,c0,c1,c2,c3,fp;
    conics_util::kepler(r0,dt,mu,u,alpha,s,c0,c1,c2,c3,fp);

    // Use the f and g linear combination to convert to state
    const double f = 1. - (mu/r0)*c2;
    const double g = dt - (mu*c3);
    const double fdot = -(mu/(fp*r0))*c1;
    const double gdot = 1. - (mu/fp)*c2;
    std::vector<double> state(6);
    state[0] = state0[0]*f    + state0[3]*g;
    state[1] = state0[1]*f    + state0[4]*g;
    state[2] = state0[2]*f    + state0[5]*g;
    state[3] = state0[0]*fdot + state0[3]*gdot;
    state[4] = state0[1]*fdot + state0[4]*gdot;
    state[5] = state0[2]*fdot + state0[5]*gdot;
    return state;
}

// Find the angle between vectors

double conics_util::vsep(const std::vector<double> &v1,const std::vector<double> &v2) {

    std::vector<double> u1(3),u2(3),vtemp(3);

    const double dmag1 = conics_util::norm(v1);
    u1[0] = v1[0]/dmag1;
    u1[1] = v1[1]/dmag1;
    u1[2] = v1[2]/dmag1;

    if(dmag1==0.0) return 0.0;

    const double dmag2 = conics_util::norm(v2);
    u2[0] = v2[0]/dmag2;
    u2[1] = v2[1]/dmag2;
    u2[2] = v2[2]/dmag2;

    if(dmag2==0.0) return 0.0;

    if(conics_util::dot(u1,u2)>0.) {
        vtemp[0] = u1[0] - u2[0];
        vtemp[1] = u1[1] - u2[1];
        vtemp[2] = u1[2] - u2[2];
        return 2.00 * asin(0.50 * conics_util::norm(vtemp));
    }

    else if(conics_util::dot(u1,u2)<0.) {
        vtemp[0] = u1[0] + u2[0];
        vtemp[1] = u1[1] + u2[1];
        vtemp[2] = u1[2] + u2[2];
        return M_PI - 2.00 * asin(0.50 * conics_util::norm(vtemp));
    }

    else {
        return 0.5*M_PI;
    }

}

// Convert state to elements
// Replicates CSPICE oscelt_c

std::vector<double> oscelt(const std::vector<double> &state,const double &et,const double &mu) {

    std::vector<double> r(3),v(3),h(3),n(3,0.),e(3),zvec(3,0.),perix(3);
    //r = list(state[:3])
    r[0] = state[0];
    r[1] = state[1];
    r[2] = state[2];
    //v = list(state[3:])
    v[0] = state[3];
    v[1] = state[4];
    v[2] = state[5];

    // Check for non-physical cases. Probably due to user input error

    const double rmag = conics_util::norm(r);
    const double vmag = conics_util::norm(v);
    //assert rmag!=0., "Position vector is zero!"
    //assert vmag!=0., "Velocity vector is zero!"

    h = conics_util::cross(r,v);

    // If the specific angular momentum vector is the zero vector,
    // we have a degenerate orbit and cannot proceed.

    //assert _norm(h)!=0., "Angular momentum is zero!"
    n[0] = -h[1];
    n[1] =  h[0];
    //n[2] = 0.;

    // Computing 2nd power

    const double d1 = (vmag*vmag) - (mu/rmag);
    const double rdv = conics_util::dot(r, v);
    //e = (d1*r - rdv*v)/mu
    e[0] = (d1*r[0] - rdv*v[0])/mu;
    e[1] = (d1*r[1] - rdv*v[1])/mu;
    e[2] = (d1*r[2] - rdv*v[2])/mu;

    // We begin by determining the size and shape of the orbit.

    double ecc = conics_util::norm(e);
    if(fabs(ecc-1.)<1e-10)
        ecc = 1.;
    const double p = conics_util::dot(h, h) / mu;
    const double rp = p / (ecc + 1.);

    // Next, the orientation of the orbit.

    //zvec = [0,0,1]
    zvec[2] = 1.;
    double inc = conics_util::vsep(h, zvec);
    if(fabs(inc)<1e-10) {
        inc = 0.;
        n[0] = 1.;
        n[1] = 0.;
        n[2] = 0.;
    }
    else if(fabs(inc-M_PI)<1e-10) {
        inc = M_PI;
        n[0] = 1.;
        n[1] = 0.;
        n[2] = 0.;
    }

    double lnode = atan2(n[1], n[0]);
    if(lnode<0.)
        lnode += 2.*M_PI;

    double argp;

    if(ecc==0.)
        argp = 0.;

    else {

        // Set the magnitude of ARGP; we'll determine the sign next.

        argp = conics_util::vsep(n, e);

        if(argp!=0.) {

            if((inc==0.) or (inc==M_PI)) {

                // The quadrant of ARGP is determined by the component of E
                // in the direction H x N.

                std::vector<double> xprod = conics_util::cross(h,n);
                const double nxprod = conics_util::norm(xprod);
                xprod[0] /= nxprod;
                xprod[1] /= nxprod;
                xprod[2] /= nxprod;
                if(conics_util::dot(e,xprod)<0.)
                    argp = 2.*M_PI - argp;
            }

            else if(e[2]<0.) {

                // The periapsis is below the reference plane;  the argument
                // of periapsis must be greater than 180 degrees.

                argp = 2.*M_PI - argp;
            }
        }
    }

    if(ecc==0.) {
        // In this case, the argument of periapse is set to zero,
        // so the nu is measured from N.
        const double nn = conics_util::norm(n);
        perix[0] = n[0]/nn;
        perix[1] = n[1]/nn;
        perix[2] = n[2]/nn;
    }

    else {
        perix[0] = e[0]/ecc;
        perix[1] = e[1]/ecc;
        perix[2] = e[2]/ecc;
    }
    
    std::vector<double> periy = conics_util::cross(h,perix);
    const double nperiy = conics_util::norm(periy);
    periy[0] /= nperiy;
    periy[1] /= nperiy;
    periy[2] /= nperiy;
    const double nu = atan2(
            conics_util::dot(r, periy),
            conics_util::dot(r, perix)
            );

    double m0;
    if(ecc<1.) {

        // For improved numerical performance, we compute both the
        // sine and cosine of the eccentric anomaly, then let ATAN2
        // find the eccentric anomaly.

        const double cosea = (ecc + cos(nu)) / (ecc * cos(nu) + 1.);
        const double sinea = (rmag / rp) * sqrt((1. - ecc) / (ecc + 1.)) * sin(nu);
        const double ea = atan2(sinea, cosea);
        m0 = fabs(ea - ecc * sin(ea));
        if(nu<0)
            m0 *= -1.;
        if(m0<0.)
            m0 += 2.*M_PI;
    }
    else if (ecc>1.) {
        const double coshf = (ecc + cos(nu)) / (ecc * cos(nu) + 1.);
        const double ea = acosh(fmax(1.,coshf));
        m0 = fabs(ecc * sinh(ea) - ea);
        if(nu<0)
            m0 *= -1.;
    }
    else {
        const double ea = tan(nu / 2.);
        // Computing 3rd power
        m0 = fabs(ea + (ea*ea*ea)/3.);
        if(mu<0)
            m0 *= -1.;
    }

    std::vector<double> elts(8);
    elts[0] = rp;
    elts[1] = ecc;
    elts[2] = inc;
    elts[3] = lnode;
    elts[4] = argp;
    elts[5] = m0;
    elts[6] = et;
    elts[7] = mu;
    return elts;
}

// Convert elements to state
// Replicates CSPICE conics_c

std::vector<double> conics(const std::vector<double> &elts,const double &et) {

    const double rp    = elts[0];
    const double ecc   = elts[1];
    const double inc   = elts[2];
    const double lnode = elts[3];
    const double argp  = elts[4];
    const double m0    = elts[5];
    const double t0    = elts[6];
    const double mu    = elts[7];

    // First construct the orthonormal basis vectors that span the orbit plane.

    const double cosi = cos(inc);
    const double sini = sin(inc);
    const double cosn = cos(lnode);
    const double sinn = sin(lnode);
    const double cosw = cos(argp);
    const double sinw = sin(argp);
    const double snci = sinn * cosi;
    const double cnci = cosn * cosi;

    const double basisp0 = cosn * cosw - snci * sinw;
    const double basisp1 = sinn * cosw + cnci * sinw;
    const double basisp2 = sini * sinw;

    const double basisq0 = -cosn * sinw - snci * cosw;
    const double basisq1 = -sinn * sinw + cnci * cosw;
    const double basisq2 = sini * cosw;

    // Next construct the state at periapse.

    const double v = sqrt(mu * (ecc + 1.) / rp);
    std::vector<double> pstate(6);
    pstate[0] = rp*basisp0;
    pstate[1] = rp*basisp1;
    pstate[2] = rp*basisp2;
    pstate[3] =  v*basisq0;
    pstate[4] =  v*basisq1;
    pstate[5] =  v*basisq2;

    // Finally compute DT the elapsed time since the epoch of periapse.
    double dt;

    // Ellipses first, since they are the most common.
    if(ecc<1.) {
        const double ainvrs = (1. - ecc) / rp;
        const double n = sqrt(mu * ainvrs) * ainvrs;
        const double period = 2*M_PI / n;
        dt = fmod(et - t0 + (m0 / n), period);
    }

    // Hyperbolas next.
    else if(ecc>1.) {

        // Again, recall that:
        //
        // N ( mean motion ) is given by DSQRT( MU / |A**3| ).
        // But since, |A| = RP / ( ECC - 1 ) ...

        const double ainvrs = (ecc - 1.) / rp;
        const double n = sqrt(mu * ainvrs) * ainvrs;
        dt = et - t0 + m0 / n;
    }

    // Finally, parabolas.
    else {
        const double n = sqrt(mu / (rp * 2.)) / rp;
        dt = et - t0 + (m0 / n);
    }

    // Now propagate the state.
    return danbyuni(pstate, mu, dt);
}

#endif
