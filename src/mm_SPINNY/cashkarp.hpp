#ifndef __CASHKARP_HPP__
#define __CASHKARP_HPP__

#include <iostream>
#include <vector>
#include <cmath>

// Template-based Cash-Karp class
//
// Typical call from within a method of My_Class
//
// rk4_class<My_Class> rk(this,&My_Class::derivative,time_step);
// rk.evolve(array_to_integrate,start_time,end_time);

const static double C[] = {0.,1./5.,3./10.,3./5.,1.,7./8.};
const static double A1[] = {1./5.};
const static double A2[] = {3./40.,9./40.};
const static double A3[] = {3./10.,-9./10.,6./5.};
const static double A4[] = {-11./54.,5./2.,-70./27.,35./27.};
const static double A5[] = {
    1631./55296.,
    175./512.,
    575./13824.,
    44275./110592.,
    253./4096.
};
const static double B2[] = {-3./2.,5./2.};
const static double B3[] = {19./54.,0.,-10./27.,55./54.};
const static double B4[] = {
    2825./27648.,
    0.,
    18575./48384.,
    13525./55296.,
    277./14336.,
    1./4.
};
const static double B5[] = {
    37./378.,
    0.,
    250./621.,
    125./594.,
    0.,
    512./1771.
};

/*template <int n> inline double nth_rootd(double x) {
    const int ebits = 11, fbits = 52;
    typedef long long int64;
    int64& i = (int64&) x;
    const int64 bias = (1 << (ebits-1))-1;
    i = (i - (bias << fbits)) / n + (bias << fbits);
    return x;
}*/

template <class T>
class cashkarp_class { public:

    // Step worked
    bool success;

    // Twiddle factors
    // Tee Hee
    double twiddle1, twiddle2;

    // Quit factors
    double quit1, quit2;

    // Step size
    double h0, h00;
    bool recover;

    // Tolerance
    double tol, tol2,tol3,tol5;

    // Safety factor
    double SF;

    // Number of equations
    unsigned neq;

    // Be verbose
    bool verbose;

    // Integration arrays
    std::vector<double> temp, k1,k2,k3,k4,k5,k6, y1,y2,y3,y4,y5, norm_arr;

    // The object to integrate
    T* obj;

    // The method of the object to integrate
    void (T::*deriv)(
                const std::vector<double> &,
                const double &, 
                std::vector<double> &) 
        const;

    // Statistics
    unsigned
        nfun,
        nstep,
        ne1,
        nu21,
        ne151,
        ne2,
        nu3,
        nu22,
        ne152,
        ne4,
        nsuc;

    // Initializers
    cashkarp_class() { init(); }

    cashkarp_class(
            T* obj0,
            void (T::*deriv0)(
                const std::vector<double> &,
                const double &,     
                std::vector<double> &) 
            const,
            const double &tol0
            ) {
        init();
        obj = obj0;
        deriv = deriv0;
        tol = tol0;
        tol2 = sqrt(tol);
        tol3 = pow(tol,1./3.);
        tol5 = pow(tol,.2);
    }

    void init() {
        twiddle1 = 1.5;
        twiddle2 = 1.1;
        quit1 = 100.0;
        quit2 = 100.0;
        h0 = 1.0;
        tol = 0.0;
        SF = 0.9;
        neq = 0;
        verbose = false;

        nfun = 0;
        nstep = 0;
        ne1 = 0;
        nu21 = 0;
        ne151 = 0;
        ne2 = 0;
        nu3 = 0;
        nu22 = 0;
        ne152 = 0;
        ne4 = 0;
        nsuc = 0;
    }

    // Allocate vectors
    void allocate(const unsigned &n);

    // A single step
    void step(std::vector<double> &y, double &t);

    // Integrate from one time to another
    void evolve(std::vector<double> &y, const double &t_start, const double &t_end);

    // Spit out stats
    void stats();
};

template <class T> 
void cashkarp_class<T>::step(
        std::vector<double> &y,
        double &t
        ) {

    success = false;
    nstep++;

    double h = h0;

    // First step
    (obj->*deriv)(y,t,k1);

    // Second step
    for(unsigned i=0;i<neq;i++) 
        temp[i] = y[i] + A1[0]*h*k1[i];
    (obj->*deriv)(temp,t+C[1]*h,k2);

    nfun += 2;

    // Find first error
    double diff1 = 0, d, ng=0;
    for(unsigned i=0;i<neq;i++) {
        if(norm_arr[i]==0.) continue;
        y1[i] = y[i] + h*k1[i];
        y2[i] = y[i] + B2[0]*h*k1[i] + B2[1]*h*k2[i];
        d = (y2[i]-y1[i])*norm_arr[i];
        diff1 += d*d;
        ng += 1.;
    }
    double err1 = sqrt(sqrt(diff1/ng));
    double E1 = err1/tol2;
    
    // Bail if terrible
    if(E1>twiddle1*quit1) {
        ne1++;
        double esttol = E1/quit1;
        h0 = fmax(.2,SF/esttol)*h0;
        if(verbose) {
            std::cout<<"# Cash Karp: Bailing on E1_1 = "<<E1<<"\n";
            std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";
        }
        return;
    }
        
    // Third step
    for(unsigned i=0;i<neq;i++)
        temp[i] = y[i] + A2[0]*h*k1[i] + A2[1]*h*k2[i];
    (obj->*deriv)(temp,t+C[2]*h,k3);
    
    // Fourth step
    for(unsigned i=0;i<neq;i++)
        temp[i] = y[i] + A3[0]*h*k1[i] + A3[1]*h*k2[i] + A3[2]*h*k3[i];
    (obj->*deriv)(temp,t+C[3]*h,k4);

    nfun += 2;

    // Compute third order solution
    double diff2 = 0;
    for(unsigned i=0;i<neq;i++) {
        if(norm_arr[i]==0.) continue;
        y3[i] = y[i];
        y3[i] += B3[0]*h*k1[i];
        //y3[i] += B3[1]*h*k2[i];
        y3[i] += B3[2]*h*k3[i];
        y3[i] += B3[3]*h*k4[i];
        d = (y3[i]-y2[i])*norm_arr[i];
        diff2 += d*d;
    }
    double err2 = sqrt(cbrt(diff2/ng));//nth_rootd<6>(diff2);
    double E2 = err2/tol3;

    if(E2>twiddle2*quit2) {
        // Try a lower order solution
        if(E1<1.) {
            // Check the error of the second order solution
            double E15 = 0.;
            for(unsigned i=0;i<neq;i++) {
                d = h*.1*(k2[i]-k1[i]);
                E15 += d*d;
            }
            E15 = sqrt(E15);
            if(E15<tol) {
                // Accept the second order solution
                nu21++;
                for(unsigned i=0;i<neq;i++)
                    y[i] += .2*h*k1[i];
                t += .2*h;
                //h0 = .2*h;
                h0 = fmax(.2,SF/E15)*h0;
                if(verbose) {
                    std::cout<<"# Cash Karp: Using second-order solution 1\n";
                    std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";
                }
                return;
            } else {
                // Bail
                ne151++;
                //h0 = h/5.;
                h0 = fmax(.2,SF/E15)*h0;
                if(verbose) {                                
                    std::cout<<"# Cash Karp: Bailing on E15_1 = "<<E15<<"\n";
                    std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";
                }
                return;
            }
        } else {
            // Bail
            ne2++;
            double esttol = E2/quit2;
            h0 = fmax(.2,SF/esttol)*h0;
            if(verbose) {
                std::cout<<"# Cash Karp: Bailing on E2_1 = "<<E2<<"\n";
                std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";
            }
            return;
        }
    }

    // Fifth step
    for(unsigned i=0;i<neq;i++) {
        temp[i] = y[i];
        temp[i] += A4[0]*h*k1[i];
        temp[i] += A4[1]*h*k2[i];
        temp[i] += A4[2]*h*k3[i];
        temp[i] += A4[3]*h*k4[i];
    }
    (obj->*deriv)(temp,t+C[4]*h,k5);

    // Sixth step
    for(unsigned i=0;i<neq;i++) {
        temp[i] = y[i];
        temp[i] += A5[0]*h*k1[i];
        temp[i] += A5[1]*h*k2[i];
        temp[i] += A5[2]*h*k3[i];
        temp[i] += A5[3]*h*k4[i];
        temp[i] += A5[4]*h*k5[i];
    }
    (obj->*deriv)(temp,t+C[5]*h,k6);

    nfun += 2;

    // Compute order 4 and 5 solutions
    double diff4 = 0;
    for(unsigned i=0;i<neq;i++) {

        if(norm_arr[i]==0.) continue;
        
        // Order 4
        y4[i] = y[i];
        y4[i] += B4[0]*h*k1[i];
        //y4[i] += B4[1]*h*k2[i];
        y4[i] += B4[2]*h*k3[i];
        y4[i] += B4[3]*h*k4[i];
        y4[i] += B4[4]*h*k5[i];
        y4[i] += B4[5]*h*k6[i];

        // Order 5
        y5[i] = y[i];
        y5[i] += B5[0]*h*k1[i];
        //y5[i] += B5[1]*h*k2[i];
        y5[i] += B5[2]*h*k3[i];
        y5[i] += B5[3]*h*k4[i];
        y5[i] += B5[4]*h*k5[i];
        y5[i] += B5[5]*h*k6[i];

        d = (y5[i]-y4[i])*norm_arr[i];
        diff4 += d*d;
    }
    double ERR4 = pow(diff4/ng,.1);//nth_rootd<10>(diff4);
    double E4 = ERR4/tol5;

    if(E4>1.) {
        // Readjust the twiddle factors
        if(E1/quit1<twiddle1)
            twiddle1 = fmax(1.1,E1/quit1);
        if(E2/quit2<twiddle2)
            twiddle2 = fmax(1.1,E2/quit2);
        if(E2<1.) {
            // Check the accuracy of the third order solution
            double E35 = 0.;
            for(unsigned i=0;i<neq;i++) {
                d = .1*h*(k1[i]-2.*k3[i]+k4[i]);
                E35 += d*d;
            }
            if(E35<tol) {
                // Accept the third order solution
                nu3++;
                for(unsigned i=0;i<neq;i++)
                    y[i] += .6*h*k3[i];
                t += .6*h;
                //h0 = .6*h;
                h0 = fmax(.6,SF/E4)*h0;
                if(verbose) {
                    std::cout<<"# Cash Karp: Using third-order solution\n";
                    std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";
                }
                return;
            } else {
                // Try a lower order solution
                if(E1<1) {
                    // Check the accuracy of the second order solution
                    double E15 = 0.;
                    for(unsigned i=0;i<neq;i++) {
                        d = h*.1*(k2[i]-k1[i]);
                        E15 += d*d;
                    }
                    E15 = sqrt(E15);
                    if(E15<tol) {
                        nu22++;
                        // Accept the second order solution
                        for(unsigned i=0;i<neq;i++)
                            y[i] += .2*h*k1[i];
                        t += .2*h;
                        //h0 = .2*h;
                        h0 = fmax(.2,SF/E15)*h0;
                        if(verbose) {
                            std::cout<<"# Cash Karp: Using second-order solution 2\n";
                            std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";
                        }
                        return;
                    } else {
                        // Bail
                        ne152++;
                        //h0 = h/5.;
                        h0 = fmax(.2,SF/E15)*h0;
                        if(verbose) {
                            std::cout<<"# Cash Karp: Bailing on E15_2 = "<<E15<<"\n";
                            std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";
                        }
                        return;
                    }
                } else {
                    // Bail
                    ne4++;
                    h0 = fmax(.2,SF/E4)*h0;
                    if(verbose) {
                        std::cout<<"# Cash Karp: Bailing on E4 = "<<E4<<"\n";
                        std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";
                    }
                    return;
                }
            }
        }
    }

    // Accept order 5 solution
    success = true;
    nsuc++;
    for(unsigned i=0;i<neq;i++) {
        y[i] = y5[i];
    }
    t += h;
    //std::cout<<"# Cash Karp: Using h0 = "<<h0<<"\n";

    h0 = fmin(5.0,SF/E4)*h;
    //std::cout<<"# Cash Karp: Setting h0 = "<<h0<<"\n";

    double Q1 = E1/E4;
    double Q2 = E2/E4;
    if(Q1>quit1)
        Q1 = fmin(Q1,10.*quit1);
    else
        Q1 = fmax(Q1,2.*quit1/3.);
    quit1 = fmax(1,fmin(10000,Q1));
    if(Q2>quit2)
        Q2 = fmin(Q2,10.*quit2);
    else
        Q2 = fmax(Q2,2.*quit2/3.);
    quit2 = fmax(1,fmin(10000,Q2));

}

//  Allocate arrays; do as infreqently as possible
template <class T>
void cashkarp_class<T>::allocate(const unsigned &n) {
    neq = n;
    temp = std::vector<double>(neq);
    k1 = std::vector<double>(neq);
    k2 = std::vector<double>(neq);
    k3 = std::vector<double>(neq);
    k4 = std::vector<double>(neq);
    k5 = std::vector<double>(neq);
    k6 = std::vector<double>(neq);
    y1 = std::vector<double>(neq);
    y2 = std::vector<double>(neq);
    y3 = std::vector<double>(neq);
    y4 = std::vector<double>(neq);
    y5 = std::vector<double>(neq);
    if(norm_arr.size()!=neq)
        norm_arr = std::vector<double>(neq,1.);
}

template <class T> 
void cashkarp_class<T>::evolve(
        std::vector<double> &y,
        const double &t_start,
        const double &t_end
        ) {

    //  Allocate arrays if necessary
    if(neq!=y.size())
        allocate(y.size());

    // Start with the default step size
    h0 = fabs(h0);
    if (t_start>t_end) 
        h0 *= -1;
    double t = t_start;

    // Step through until done
    //double nreset = 0;
    double minh = 1.0e-5;
    while(t!=t_end) {

        // Reduce step size if necessary
        if ((t_end-t)/h0<1.) {
            h00 = h0;
            h0 = t_end-t;
            recover = true;
        }
        else {
            recover = false;
        }

        // The magic
        step(y,t);

        // If we had reduced the step size, and it worked, recover
        if(recover and success) {
            h0 = h00;
        }
        //if(success) {
        //    std::cout<<"Step\n";
        //    std::cout<<"Timestep "<<h0<<"\n";
        //}
        //if(!success) {
        //    nreset++;
        //    //std::cout<<"# Fail\n";
        //    //std::cout<<"Timestep "<<h0<<"\n";
        //}
        //if(nreset>10000.) {
        //    y[0] = NAN;
        //    t = t_end;
        //    std::cout<<"Quitting integration";
        //}
        if(fabs(h0) < minh){
            y[0] = NAN;
            t = t_end;
            std::cout<<"Equations of motion are stiff. Quitting integration.\n";
        }

    }

}

template <class T>
void cashkarp_class<T>::stats() {
    //nstep = 0;
    //ne1 = 0;
    //nu21 = 0;
    //ne151 = 0;
    //ne2 = 0;
    //nu3 = 0;
    //nu22 = 0;
    //ne152 = 0;
    //ne4 = 0;
    //nsuc = 0;
    double psuc = 100.*(float)(nsuc)/(float)(nstep);
    unsigned nfail = ne1+ne151+ne2+ne152+ne4;
    double pfail = 100.*(float)(nfail)/(float)(nstep);
    unsigned ninter = nu21+nu3+nu22;
    double pinter = 100.*(float)(ninter)/(float)(nstep);
    //std::cout<<std::fixed;
    std::cout<<"# Cash-Karp Statistics:\n#\n";
    std::cout<<"# Number of function calls:   "<<nfun<<"\n";
    std::cout<<"# Number of steps:            "<<nstep<<"\n";
    std::cout<<"# Number of successful steps: "<<nsuc<<", "<<psuc<<"%\n";
    std::cout<<"# Number of failed steps:     "<<nfail<<", "<<pfail<<"%\n";
    std::cout<<"# Number of reduced steps:    "<<ninter<<", "<<pinter<<"%\n";
}

#endif
