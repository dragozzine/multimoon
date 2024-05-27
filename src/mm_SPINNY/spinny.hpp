#ifndef __SPINNY_HPP__
#define __SPINNY_HPP__
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <map>
#include "cashkarp.hpp"

//##### Object Physical Properties #####

class Physical_Properties { public:

    double mass;
    double J2R2,C22R2;
    double ellip_a,ellip_b,ellip_c;
    std::vector<double> I,iI;
    double Ic0,Ic1,Ic2;

    Physical_Properties() {
        I = std::vector<double>(3,0);
        iI = std::vector<double>(3,0);
    }

    Physical_Properties(double mass0,std::string mode,
            std::vector<double> I=std::vector<double>(3,0),
            double a=0,double b=0,double c=0,double J2R2=0,double C22R2=0) {
        if(mode=="I") {
            set_I(mass,I);
        }
        else if(mode=="ellip") {
            set_ellip(mass,a,b,c);
        }
        else if(mode=="grav") {
            set_grav(mass,J2R2,C22R2,c);
        }
    }

    void set_I(double mass0,const std::vector<double> &I0) {
        mass = mass0;
        for(unsigned i=0;i<3;i++)
            I[i] = I0[i];
        J2R2 = (I[2] - I[0]/2.0 - I[1]/2.0)/mass;
        C22R2 = (I[1]/4.0 - I[0]/4.0)/mass;
        ellip_a = sqrt(5.0/(2.0*mass)*(I[1] + I[2] - I[0]));
        ellip_b = sqrt(5.0/(2.0*mass)*(I[0] + I[2] - I[1]));
        ellip_c = sqrt(5.0/(2.0*mass)*(I[0] + I[1] - I[2]));
        set_Ic();
    }

    void set_ellip(double mass0,double a,double b,double c) {
        mass = mass0;
        ellip_a = a;
        ellip_b = b;
        ellip_c = c;
        J2R2 = 0.2*(0.5*a*a + 0.5*b*b - c*c);
        C22R2 = 0.05*(a*a - b*b);
        I[0] = mass/5.0*(b*b + c*c);
        I[1] = mass/5.0*(a*a + c*c);
        I[2] = mass/5.0*(a*a + b*b);
        set_Ic();
    }

    void set_grav(double mass0,double J2R20,double C22R20,double c) {
        mass = mass0;
        J2R2 = J2R20;
        C22R2 = C22R20;
        ellip_a = sqrt(5.0*J2R2 + c*c + 10.0*C22R2);
        ellip_b = sqrt(5.0*J2R2 + c*c - 10.0*C22R2);
        ellip_c = c;
        double a = ellip_a;
        double b = ellip_b;
        I[0] = mass/5.0*(b*b + c*c);
        I[1] = mass/5.0*(a*a + c*c);
        I[2] = mass/5.0*(a*a + b*b);
        set_Ic();
    }

    void set_Ic() {
        if(fabs(I[0])+fabs(I[1])+fabs(I[2])==0.) {
            Ic0 = 0.;
            Ic1 = 0.;
            Ic2 = 0.;
            iI = std::vector<double>(3,0.);
        } else {
            Ic0 = (I[1]-I[2])/I[0];
            Ic1 = (I[2]-I[0])/I[1];
            Ic2 = (I[0]-I[1])/I[2];
            iI[0] = 1./I[0];
            iI[1] = 1./I[1];
            iI[2] = 1./I[2];
        }
    }

    void print_line() const {
        //print('[{:8.3f} {:8.3f} {:8.3f}] {:8.3f} {:8.3f} [{:8.3f} {:8.3f} {:8.3f}]'.format(
        std::cout<<"["<<ellip_a<<" "<<ellip_b<<" "<<ellip_c<<"]";
        std::cout<<" "<<J2R2<<" "<<C22R2<<" ";
        std::cout<<"["<<I[0]<<" "<<I[1]<<" "<<I[2]<<"]\n";
    }

};

//##### Spinny Itself #####

class Spinny { public:

    // Number of objects
    unsigned nobj;

    // Names
    std::vector<std::string> names;
    std::map<std::string,unsigned> idx;

    // Physical Properties
    std::vector<Physical_Properties> phys;

    // The integration array
    std::vector<double> arr0;

    // Do the objects have spin
    std::vector<bool> hasspin;

    // The time
    double t;

    // The integrator
    cashkarp_class<Spinny> ck;

    // Generic init function
    void init() {
        nobj = 0;
    }

    // Base contructor
    Spinny() { init(); }

    // Fancy constructor
    Spinny(double t0,double tol,double h00=1) {
        init();
        t = t0;
        ck = cashkarp_class<Spinny>(this,&Spinny::derivative,tol);
        ck.h0 = h00;
    }

    // Add an object to the system, with spin
    void add_object(const std::string &name,const Physical_Properties &phys,
            const std::vector<double> &state,
            const std::vector<double> &spin,
            const std::vector<double> &quaternion);

    // Add an object to the system, without spin
    void add_object(const std::string &name,const Physical_Properties &phys0,
            const std::vector<double> &state);

    // Extract information
    std::vector<double> get_position(unsigned i) const;
    std::vector<double> get_velocity(unsigned i) const;
    std::vector<double> get_state(unsigned i) const;
    std::vector<double> get_state(unsigned i,unsigned k) const;
    std::vector<double> get_spin(unsigned i) const;
    std::vector<double> get_quaternion(unsigned i) const;

    // Calculate the barycentre
    void calc_bary(double &ub,std::vector<double> &rb,std::vector<double> &vb) const;

    // Move the barycentre to the origin
    void move2bary();

    // The time derivative
    void derivative(const std::vector<double> &arr,const double &time,std::vector<double> &deriv) const;

    // Point-particle gravitation
    void particle_grav(const std::vector<double> &arr,const unsigned &i,const unsigned &j,std::vector<double> &a) const;

    void q2mat(const std::vector<double> &arr,const unsigned &i,double rot[3][3]) const;
    void transpose(double mat1[3][3],double mat2[3][3]) const;

    // Gravity from oblate objects
    void oblate_grav(const std::vector<double> &arr,const unsigned &i,const unsigned &j,double rot[3][3],
            std::vector<double> &darr) const;

    // Evovle the system in time
    void evolve(const double &t1);

    // Normalize things
    void norm();
};

// Add an object to the system, with spin
void Spinny::add_object(const std::string &name,const Physical_Properties &phys0,const std::vector<double> &state,
        const std::vector<double> &spin,const std::vector<double> &quaternion) {

    // Save the name and physical properties
    names.push_back(name);
    idx[name] = names.size()-1;
    phys.push_back(phys0);

    // Add the state vector
    arr0.insert(arr0.end(), state.begin(), state.end());

    // Add the spin vector
    arr0.insert(arr0.end(), spin.begin(), spin.end());

    // Add the orientation quaternion
    arr0.insert(arr0.end(), quaternion.begin(), quaternion.end());

    // Do we need to spin math on the object
    if(phys0.J2R2==0 and phys0.C22R2==0)
        hasspin.push_back(false);
    else
        hasspin.push_back(true);

    nobj++;
}

// Add an object to the system, without spin
void Spinny::add_object(const std::string &name,const Physical_Properties &phys0,const std::vector<double> &state) {
    std::vector<double> spin(3,0.);
    std::vector<double> quaternion(4,0.);
    quaternion[0] = 1.;
    add_object(name,phys0,state,spin,quaternion);
}

// Get position vector
std::vector<double> Spinny::get_position(unsigned i) const {
    std::vector<double> vec(3);
    for(unsigned j=0;j<3;j++)
        vec[j] = arr0[i*13 + j];
    return vec;
}

// Get absolute velocity vector
std::vector<double> Spinny::get_velocity(unsigned i) const {
    std::vector<double> vec(3);
    for(unsigned j=0;j<3;j++)
        vec[j] = arr0[i*13 + j + 3];
    return vec;
}

// Get relative state vector
std::vector<double> Spinny::get_state(unsigned i,unsigned k) const {
    std::vector<double> vec(6);
    for(unsigned j=0;j<6;j++)
        vec[j] = arr0[i*13 + j] - arr0[k*13 + j];
    return vec;
}

// Get state vector
std::vector<double> Spinny::get_state(unsigned i) const {
    std::vector<double> vec(6);
    for(unsigned j=0;j<6;j++)
        vec[j] = arr0[i*13 + j];
    return vec;
}

// Get spin vector
std::vector<double> Spinny::get_spin(unsigned i) const {
    std::vector<double> vec(3);
    for(unsigned j=0;j<3;j++)
        vec[j] = arr0[i*13 + j + 6];
    return vec;
}

// Get orientation quaternion
std::vector<double> Spinny::get_quaternion(unsigned i) const {
    std::vector<double> vec(4);
    for(unsigned j=0;j<4;j++)
        vec[j] = arr0[i*13 + j + 9];
    return vec;
}

// Calculate the barycentre
void Spinny::calc_bary(double &ub,std::vector<double> &rb,std::vector<double> &vb) const {
    ub = 0.;
    rb = std::vector<double>(3,0.);
    vb = std::vector<double>(3,0.);
    unsigned i,i2,j;
    for(i=0;i<nobj;i++) {
        i2 = i*13;
        ub += phys[i].mass;
        for(j=0;j<3;j++) {
            rb[j] += arr0[i2+j  ]*phys[i].mass;
            vb[j] += arr0[i2+j+3]*phys[i].mass;
        }
    }
    for(j=0;j<3;j++) {
        rb[j] /= ub;
        vb[j] /= ub;
    }
}

// Move the barycentre
void Spinny::move2bary() {
    double ub;
    std::vector<double> rb,vb;
    calc_bary(ub,rb,vb);
    unsigned i,i2,j;
    for(i=0;i<nobj;i++) {
        i2 = i*13;
        for(j=0;j<3;j++) {
            arr0[i2+j  ] -= rb[j];
            arr0[i2+j+3] -= vb[j];
        }
    }
}

// The time derivative function
void Spinny::derivative(
        const std::vector<double> &arr,
        const double &t,
        std::vector<double> &darr) const {

    // Decompress arrays
    //r,v,s,q = [],[],[],[]
    /*for i in range(self.nobj):
            i2 = 13*i
            //r.append(arr[i2  :i2+ 3])

            v.append(arr[i2+3:i2+ 6])
            s.append(arr[i2+6:i2+ 9])
            q.append(arr[i2+9:i2+13])
        r = np.asarray(r)
        v = np.asarray(v)
        s = np.asarray(s)
        q = np.asarray(q)*/

    //std::vector<double> a(nobj*3);
    //std::vector<double> T(nobj*3);
    //std::vector<double> O,dO;
    //O.reserve(3); dO.reserve(3);
    double rot[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double O[3];
    unsigned i,i2,j;

    for(i=0;i<nobj*13;i++) darr[i]=0;

    // Evaulate forces and torques
    for(i=0;i<nobj;i++) {

        // Check if a non-perturbing particle
        if(phys[i].mass==0) continue;

        // If both of those terms are zero, treat object i as a particle
        //has_oblate = (self.phys[i].J2R2+self.phys[i].C22R2 > 0)

        if(hasspin[i]) {
            q2mat(arr,i,rot);
        }

        for(j=0;j<nobj;j++) {

            if(i==j) continue;

            // Full oblate gravity
            else if(hasspin[i]) {
                oblate_grav(arr,i,j,rot,darr);
            }

            // Point-particle gravity
            else {
                particle_grav(arr,i,j,darr);
            }

            //### Can add other physics (i.e. tides, rings) here ###
        }
    }

    double qr,qi,qj,qk;

    // Convert to appropriate derivatives
    for(i=0;i<nobj;i++) {
        i2 = 13*i;
        //i3 = 3*i;

        for(j=0;j<3;j++)
            darr[i2+j  ] = arr[i2+j+3];

        if(hasspin[i]) {

            O[0] = arr[i2+6];
            O[1] = arr[i2+7];
            O[2] = arr[i2+8];

            // Euler's equations of motion
            //
            // Note that we stored the torque in darr[i2+6..8],
            // and are now converting to differential spin
            //
            // BP 11/10/21: Flipped the signs for Euler's equations. These are now consistent with the equations
            //              analytically dervived. This configuration is the only one to conserve total L.

            darr[i2+6] =  phys[i].Ic0*O[1]*O[2] + (darr[i2+6]*phys[i].iI[0]);
            darr[i2+7] =  phys[i].Ic1*O[2]*O[0] + (darr[i2+7]*phys[i].iI[1]);
            darr[i2+8] =  phys[i].Ic2*O[0]*O[1] + (darr[i2+8]*phys[i].iI[2]);

            // Matrix operator for dqdt
            /*Omat = np.array(
                    [[  0.,-O[0],-O[1],-O[2]],
                    [O[0],   0., O[2],-O[1]],
                    [O[1],-O[2],   0., O[0]],
                    [O[2], O[1],-O[0],   0.]]
                    )
            dq = 0.5*np.dot(Omat,np.transpose(q[i]))*/
            qr = arr[i2+ 9];
            qi = arr[i2+10];
            qj = arr[i2+11];
            qk = arr[i2+12];

            // Original code
            darr[i2+ 9] = .5*(        -O[0]*qi -O[1]*qj -O[2]*qk);
            darr[i2+10] = .5*(O[0]*qr          +O[2]*qj -O[1]*qk);
            darr[i2+11] = .5*(O[1]*qr -O[2]*qi          +O[0]*qk);
            darr[i2+12] = .5*(O[2]*qr +O[1]*qi -O[0]*qj         );


            /*darr[i2  :i2+ 3] = v[i]
            darr[i2+3:i2+ 6] = a[i]
            darr[i2+6:i2+ 9] = dO
            darr[i2+9:i2+13] = dq*/

        }
    }
}

// Point-particle gravitation of object i on object j
void Spinny::particle_grav(const std::vector<double> &arr,const unsigned &i,const unsigned &j,
        std::vector<double> &darr) const {
    const unsigned i2 = 13*i, j2 = 13*j;
    const double dr0 = arr[j2  ]-arr[i2  ];
    const double dr1 = arr[j2+1]-arr[i2+1];
    const double dr2 = arr[j2+2]-arr[i2+2];
    const double drs = (dr0*dr0)+(dr1*dr1)+(dr2*dr2);
    const double coef = -phys[i].mass/(drs*sqrt(drs));
    darr[j2+3] += coef*dr0;
    darr[j2+4] += coef*dr1;
    darr[j2+5] += coef*dr2;
}

// Convert a quaternion to a rotation matrix
void Spinny::q2mat(const std::vector<double> &arr,const unsigned &i,double rot[3][3]) const {
    // This should convert the body frame to the world frame (pretty sure not vice versa?)
    //qr,qi,qj,qk = q
    const double qr = arr[13*i +  9];
    const double qi = arr[13*i + 10];
    const double qj = arr[13*i + 11]; 
    const double qk = arr[13*i + 12];
    rot[0][0]=1.-2.*(qj*qj+qk*qk); rot[0][1]=   2.*(qi*qj-qk*qr); rot[0][2]=   2.*(qi*qk+qj*qr);
    rot[1][0]=   2.*(qi*qj+qk*qr); rot[1][1]=1.-2.*(qi*qi+qk*qk); rot[1][2]=   2.*(qj*qk-qi*qr);
    rot[2][0]=   2.*(qi*qk-qj*qr); rot[2][1]=   2.*(qj*qk+qi*qr); rot[2][2]=1.-2.*(qi*qi+qj*qj);
}

// Transpose a matrix
void Spinny::transpose(double mat1[3][3],double mat2[3][3]) const {
    mat2[0][0] = mat1[0][0]; mat2[0][1] = mat1[1][0]; mat2[0][2] = mat1[2][0];
    mat2[1][0] = mat1[0][1]; mat2[1][1] = mat1[1][1]; mat2[1][2] = mat1[2][1];
    mat2[2][0] = mat1[0][2]; mat2[2][1] = mat1[1][2]; mat2[2][2] = mat1[2][2];
}

// Gravity from oblate objects
void Spinny::oblate_grav(const std::vector<double> &arr,const unsigned &i,const unsigned &j,
        double rot[3][3],std::vector<double> &darr) const {

    const double i2=13*i, j2=13*j;
    //double goodrot[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

    // Rotate to body frame
    //dr0 = r[j]-r[i] in the world frame
    //rot = q2mat(q[i]).T
    // THIS CODE DOESN'T DO THE TRANSPOSE!
    //dr[= np.dot(rot,dr0)
    const double dr00 = arr[j2  ]-arr[i2  ];
    const double dr01 = arr[j2+1]-arr[i2+1];
    const double dr02 = arr[j2+2]-arr[i2+2];

    //transpose(rot,goodrot);
    const double dr10 = rot[0][0]*dr00 + rot[1][0]*dr01 + rot[2][0]*dr02;
    const double dr11 = rot[0][1]*dr00 + rot[1][1]*dr01 + rot[2][1]*dr02;
    const double dr12 = rot[0][2]*dr00 + rot[1][2]*dr01 + rot[2][2]*dr02;

    //const double dr10 = goodrot[0][0]*dr00 + goodrot[1][0]*dr01 + goodrot[2][0]*dr02;
    //const double dr11 = goodrot[0][1]*dr00 + goodrot[1][1]*dr01 + goodrot[2][1]*dr02;
    //const double dr12 = goodrot[0][2]*dr00 + goodrot[1][2]*dr01 + goodrot[2][2]*dr02;

    //COMPARE TO BELOW:
    // Rotate force back to world frame
    //aw = np.dot(rot.T,ab)
    // ab is also in the world frame
    //    darr[j2+3] += rot[0][0]*ab0 + rot[0][1]*ab1 + rot[0][2]*ab2;
    //    darr[j2+4] += rot[1][0]*ab0 + rot[1][1]*ab1 + rot[1][2]*ab2;
    //    darr[j2+5] += rot[2][0]*ab0 + rot[2][1]*ab1 + rot[2][2]*ab2;



    // Start with point-particle force
    //dr2 = np.sum(np.square(dr))
    //dr3 = dr2*np.sqrt(dr2)
    //ab = (-self.phys[i].mass*dr)/dr3
    const double dr2 = (dr10*dr10)+(dr11*dr11)+(dr12*dr12);
    const double dr3 = dr2*sqrt(dr2);
    const double coef = -phys[i].mass/dr3;
    //double ab0 = coef*dr10;
    //double ab1 = coef*dr11;
    //double ab2 = coef*dr12;

    // Add oblate force
    const double dr5 = dr3*dr2;
    const double dr7 = dr3*dr2*dr2;

    //A,B,C = self.phys[i].I
    const double A = phys[i].I[0];
    const double B = phys[i].I[1];
    const double C = phys[i].I[2];
    //ab0 += ((B+C-2.*A)*dr10)/dr5;
    //ab1 += ((C+A-2.*B)*dr11)/dr5;
    //ab2 += ((A+B-2.*C)*dr12)/dr5;
    const double f = (B+C-2.*A)*dr10*dr10 + (C+A-2.*B)*dr11*dr11 + (A+B-2.*C)*dr12*dr12;
    //const double f = 0;
    //const double ab0 = (coef*dr10 + ((B+C-2.*A)*dr10)/dr5);
    //const double ab1 = (coef*dr11 + ((C+A-2.*B)*dr11)/dr5);
    //const double ab2 = (coef*dr12 + ((A+B-2.*C)*dr12)/dr5);
    const double ab0 = (coef*dr10 + ((B+C-2.*A)*dr10)/dr5) - ((5.0/2.0)*f*dr10/dr7);
    const double ab1 = (coef*dr11 + ((C+A-2.*B)*dr11)/dr5) - ((5.0/2.0)*f*dr11/dr7);
    const double ab2 = (coef*dr12 + ((A+B-2.*C)*dr12)/dr5) - ((5.0/2.0)*f*dr12/dr7);

    // Rotate force back to world frame
    //aw = np.dot(rot.T,ab)
    // TRANPOSE OF GOODROT IS JUST ROT
    darr[j2+3] += rot[0][0]*ab0 + rot[0][1]*ab1 + rot[0][2]*ab2;
    darr[j2+4] += rot[1][0]*ab0 + rot[1][1]*ab1 + rot[1][2]*ab2;
    darr[j2+5] += rot[2][0]*ab0 + rot[2][1]*ab1 + rot[2][2]*ab2;

    // Calculate torques in i's body frame and temporarily store them in darr
    // BP 09/05/23: Removed that coefficient and sticking with just the mass. Let's hope this one's right.... sigh
    const double coef2 = phys[j].mass;
    //const double coef2 = phys[j].mass*rat;
    // BP 11/10/21: Switched = for +=. This allows for calculations of torques from an arbitrary number of bodies.
    darr[i2+6] += coef2*(3.*(C-B)*dr11*dr12)/dr5;
    darr[i2+7] += coef2*(3.*(A-C)*dr12*dr10)/dr5;
    darr[i2+8] += coef2*(3.*(B-A)*dr10*dr11)/dr5;

}

// Integrate the system in time
void Spinny::evolve(const double &t1) {
    if(t1==t) return;
    norm();
    ck.evolve(arr0,t,t1);
    t = t1;
}

// Normalize
void Spinny::norm() {

    // Integration normalization
    ck.norm_arr = std::vector<double>(arr0.size(),0);
    for(unsigned i=0;i<nobj;i++) {
        unsigned i2 = 13*i;
        double ir = 1./sqrt(arr0[i2  ]*arr0[i2  ] + arr0[i2+1]*arr0[i2+1] + arr0[i2+2]*arr0[i2+2]);
        double iv = 1./sqrt(arr0[i2+3]*arr0[i2+3] + arr0[i2+4]*arr0[i2+4] + arr0[i2+5]*arr0[i2+5]);
        ck.norm_arr[i2  ] = ir;
        ck.norm_arr[i2+1] = ir;
        ck.norm_arr[i2+2] = ir;
        ck.norm_arr[i2+3] = iv;
        ck.norm_arr[i2+4] = iv;
        ck.norm_arr[i2+5] = iv;
        if(hasspin[i]) {
            double is = 1./sqrt(arr0[i2+6]*arr0[i2+6] + arr0[i2+7]*arr0[i2+7] + arr0[i2+8]*arr0[i2+8]);
            ck.norm_arr[i2+ 6] = is;
            ck.norm_arr[i2+ 7] = is;
            ck.norm_arr[i2+ 8] = is;
            ck.norm_arr[i2+ 9] = 1.;
            ck.norm_arr[i2+10] = 1.;
            ck.norm_arr[i2+11] = 1.;
            ck.norm_arr[i2+12] = 1.;
        }
    }

    // Normalize quaternions
    for(unsigned i=0;i<nobj;i++) {
        unsigned i2 = 13*i;
        if(not hasspin[i]) continue;
        double iq = 1./sqrt(arr0[i2+9]*arr0[i2+9] + arr0[i2+10]*arr0[i2+10]
                + arr0[i2+11]*arr0[i2+11] + arr0[i2+12]*arr0[i2+12] );
        arr0[i2+ 9] *= iq;
        arr0[i2+10] *= iq;
        arr0[i2+11] *= iq;
        arr0[i2+12] *= iq;
    }
}

#endif
