# distutils: language = c++
# distutils: extra_compile_args = -O3

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map

import numpy as np
#cimport numpy as np

cdef extern from "cashkarp.hpp":

    cdef cppclass cashkarp_class[T]:

        cashkarp_class() except +

        void stats()

cdef extern from "spinny.hpp":

    cdef cppclass cpp_Physical_Properties "Physical_Properties":
        double mass, J2R2,C22R2, ellip_a,ellip_b,ellip_c #, Ic0,Ic1,Ic2
        vector[double] I,iI
        cpp_Physical_Properties() except +
        void set_I(double,vector[double])
        void set_ellip(double,double,double,double)
        void set_grav(double,double,double,double)
        void print_line()

    cdef cppclass Spinny:
    
        # Number of objects
        #unsigned nobj;

        #// Names
        #std::vector<std::string> names;
        map[string,unsigned] idx

        #// Physical Properties
        #std::vector<cpp_Physical_Properties> phys;

        # The integration array
        vector[double] arr0

        # Do the objects have spin
        #vector[bool] hasspin

        # The time
        double t

        # The integrator
        cashkarp_class[Spinny] ck

        #// Generic init function
        #void init() {
        #    nobj = 0;
        #}

        # Base contructor
        Spinny() except +

        # Fancy constructor
        Spinny(double t0,double tol,double h00) except +

        # Add an object to the system, with spin
        void add_object(string name,cpp_Physical_Properties phys, vector[double] state,
                vector[double] spin, vector[double] quaternion)

        # Add an object to the system, without spin
        void add_object(string name,cpp_Physical_Properties phys0,vector[double] state)

        # Extract information
        vector[double] get_position(unsigned i)
        vector[double] get_velocity(unsigned i)
        vector[double] get_state(unsigned i)
        vector[double] get_state(unsigned i,unsigned k)
        vector[double] get_spin(unsigned i)
        vector[double] get_quaternion(unsigned i)

        # Calculate the barycentre
        #void calc_bary(double &ub,std::vector<double> &rb,std::vector<double> &vb) const;

        # Move the barycentre to the origin
        void move2bary()

        #// The time derivative
        #void derivative(const std::vector<double> &arr,const double &time,std::vector<double> &deriv) const;

        #// Point-particle gravitation
        #void particle_grav(const std::vector<double> &arr,const unsigned &i,const unsigned &j,std::vector<double> &a) const;

        void q2mat(vector[double] arr, unsigned i,double rot[3][3])
        #void transpose(double mat1[3][3],double mat2[3][3]) const;

        #// Gravity from oblate objects
        #void oblate_grav(const std::vector<double> &arr,const unsigned &i,const unsigned &j,double rot[3][3],
        #        std::vector<double> &darr) const;

        # Evovle the system in time
        void evolve(double t1)

        #// Normalize things
        #void norm();


cdef class Physical_Properties:
    
    cdef cpp_Physical_Properties *thisptr
    
    def __cinit__(self, mass=0, mode=None, I=0, a=0, b=0, c=0, J2R2=0, C22R2=0):
        self.thisptr = new cpp_Physical_Properties()
        if self.thisptr is NULL:
            raise MemoryError()
        if mode is None:
            return
        mode = mode.lower()
        if mode=="i": 
            self.set_I(mass,I)
        elif mode=="ellip":
            self.set_ellip(mass,a,b,c)
        elif mode=="grav":
            self.set_grav(mass,J2R2,C22R2,c)
        
    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    def set_I(self, double mass0,vector[double] I0):
        self.thisptr.set_I(mass0, I0)

    def set_ellip(self, double mass0,double a,double b,double c):
        self.thisptr.set_ellip(mass0, a, b, c)

    def set_grav(self, double mass0,double J2R20,double C22R20,double c):
        self.thisptr.set_grav(mass0, J2R20, C22R20, c)

    def print_line(self):
        self.thisptr.print_line()

    @property
    def mass(self): return self.thisptr.mass

    @property
    def J2R2(self): return self.thisptr.J2R2

    @property
    def C22R2(self): return self.thisptr.C22R2

    @property
    def ellip_a(self): return self.thisptr.ellip_a

    @property        
    def ellip_b(self): return self.thisptr.ellip_b

    @property        
    def ellip_c(self): return self.thisptr.ellip_c

    @property
    def I(self): return np.array(self.thisptr.I)

    @property
    def iI(self): return np.array(self.thisptr.iI)


cdef class Spinny_System:

    cdef Spinny *thisptr
    
    def __cinit__(self, double t0=0, double tol=1e-9, double h0=1):
        self.thisptr = new Spinny(t0,tol,h0)
        if self.thisptr is NULL:
            raise MemoryError()
    
    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    def add_object(self, name,Physical_Properties phys0,vector[double] state, spin = None, quaternion = None):
        try: name = name.encode('UTF-8')
        except: pass
        if spin is None:
            self.thisptr.add_object(name,phys0.thisptr[0],state)
        else:
            self.thisptr.add_object(name,phys0.thisptr[0],state,spin,quaternion)

    def move2bary(self):
        self.thisptr.move2bary()

    def evolve(self, double t):
        self.thisptr.evolve(t)

    def get_state(self, unsigned i, j=None):
        if j is None:
            return np.array(self.thisptr.get_state(i))
        else:
            return np.array(self.thisptr.get_state(i,j))

    def get_spin(self, i):
        return np.array(self.thisptr.get_spin(i))

    def get_quaternion(self,i):
        return np.array(self.thisptr.get_quaternion(i))

    def idx(self, name):
        try: name = name.encode('UTF-8')
        except: pass
        return int(self.thisptr.idx[name])

    def q2mat(self, arr, i):
        cdef double rot[3][3];
        self.thisptr.q2mat(arr,i,rot)
        return np.array(rot)
    
    def ck_stats(self):
        self.thisptr.ck.stats()

    @property
    def t(self):
        return self.thisptr.t

    @property
    def arr0(self):
        return np.array(self.thisptr.arr0)


cdef extern from "simplespin.hpp":

    cdef cppclass cpp_SimpleSpinOrbit "SimpleSpinOrbit":

        double n,M0,t,mu,P

        cashkarp_class[cpp_SimpleSpinOrbit] ck

        vector[double] arr0

        cpp_SimpleSpinOrbit(double a0,double e0,double theta0,double dtheta0,double f0,
                cpp_Physical_Properties phys10,cpp_Physical_Properties phys20) except +

        void evolve(double t1)

cdef class SimpleSpinOrbit:

    cdef cpp_SimpleSpinOrbit *thisptr

    def __cinit__(self, double a0,double e0,double theta0,double dtheta0,double f0,
            Physical_Properties phys10,Physical_Properties phys20):
        self.thisptr = new cpp_SimpleSpinOrbit(a0,e0,theta0,dtheta0,f0,
                            phys10.thisptr[0],phys20.thisptr[0])
        if self.thisptr is NULL:
            raise MemoryError()

    def evolve(self, double t1):
        self.thisptr.evolve(t1)

    def set_w(self, double w):
        self.thisptr.arr0[0] = w

    def set_t(self, double t):
        self.thisptr.t = t

    def ck_stats(self):
        self.thisptr.ck.stats()

    @property
    def n(self):
        return self.thisptr.n

    @property
    def P(self):
        return self.thisptr.P

    @property
    def M0(self):
        return self.thisptr.M0
    
    @property
    def t(self):
        return self.thisptr.t
    
    @property
    def mu(self):
        return self.thisptr.mu

    @property
    def arr0(self):
        return np.array(self.thisptr.arr0)
