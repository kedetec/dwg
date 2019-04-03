#
#
import numpy as np
import cmath

class Dwg:

    def __init__(self):
        """
        input parameters from file
        """
        self.nlay = 3
        self.thickness = [2.5, 0.3 ,4.0]
        self.nindx = [3.3, 3.4, 3.3]

        self.wavelength = 1.5
        self.nig = 10.0

        #init guess
        self.nreff = 0.0
        self.nieff = 0.0
     
    def init(self):
        """
        initlize, set some parameters
        """
        self.k0 = 2.0e6*np.pi/self.wavelength
        
        nmax = self.nindx[0]
        nmin = self.nindx[0]
        for i in range(self.nlay):
            if self.nindx[i] > nmax:
                nmax = self.nindx[i]
            if self.nindx[i] < nmin:
                nmin = self.nindx[i]
                
        
        self.dkz = self.k0 * (nmax - nmin)/self.nig
        self.kz = self.k0 * nmax

        self.kzi = [i for i in range(4)]
        #self.kzi = [i for i in range(self.nlay)]
        self.c0 = complex(0.0,0.0)
        self.c1 = complex(1.0, 0.0)
        self.tmatrx = np.array([[self.c0,self.c0],[self.c0,self.c0]])
        self.mmatrx = np.array([[self.c0,self.c0],[self.c0,self.c0]])
        self.smatrx = np.array([[self.c0,self.c0],[self.c0,self.c0]])
        self.cox = []
        t0 = 0.0
        self.cox.append(t0)
        for i in range(1,self.nlay):
            t0 = t0 + self.thickness[i]*1.0e-6
            self.cox.append(t0)
            
    def print(self):
        """
        print("layers  ","Thickness ","RefractiveIndx")
        for i in range(self.nlay):
            print(i+1,self.thickness[i],self.nindx[i])      
            print(self.cox[i])
        """
        print('k0:',self.k0)
        print('dkz:',self.dkz)

    def initgus(self):
        self.kz = self.kz - self.dkz
        c0 = self.funct(self.kz)
        self.kz = self.kz - self.dkz
        c1 = self.funct(self.kz)

        """
        iorr = 1: pure imaginary function
               2: pure real function
               3: complex            function    
        """
        if c0.real * c0.imag == 0:
            if c0.real == 0:
                iorr = 1
            else:
                iorr = 2
        else:
            if abs(c0.imag/c0.real) > 1.0e8:
                iorr = 1
            elif abs(c0.real/c0.imag) > 1.0e8:
                iorr = 2
            else:
                iorr = 3
                
        if iorr == 1:
            x0 = c0.imag
            x1 = c1.imag
        elif iorr == 2:
            x0 = c0.real
            x1 = c1.real
        else:
            x0 = abs(c0)
            x1 = abs(c1)
    
        for i in range(2,int(self.nig)-1):
            self.kz = self.kz - self.dkz
            c2 = self.funct(self.kz)
            
            if c2.real*c2.imag == 0:
                if c2.real == 0:
                    x2 = c2.imag
                else:
                    x2 = c2.real
            else:
                if abs(c2.imag/c2.real) > 1.0e8:
                    x2 = c2.imag
                elif abs(c2.real/c2.imag) > 1.0e8:
                    x2 = c2.real
                else:
                    x2 = abs(c2)
            
            dx0 = x0 - x1
            dx1 = x1 - x2

            ans = self.c0
            if x0*x1 < 0.0:
                ans = (self.kz+self.dkz+self.dkz)/self.k0
                  
            if (dx0*dx1) < 0.0:
                if x1<0.0 and dx0<0.0:
                    ans = (self.kz+self.dkz)/self.k0
                elif x1>0.0 and dx0>0.0:
                    ans = (self.kz+self.dkz)/self.k0
            if abs(ans) != 0.0:
                print("Answer: ",ans)
                print('x0:',x0,'x1:',x1,'x2:',x2)

            x0 = x1
            x1 = x2

    def funct(self,kz):

        c0 = complex(0.0,0.0)
        c1 = complex(1.0,0.0)
        cj = complex(0.0,1.0)
        su = np.array([[c1,c0],[c0,c1]])
        sl = np.array([[c1,c0],[c0,c1]])
        tm0 = np.array([[c0,c0],[c0,c0]])
        tm1 = np.array([[c0,c0],[c0,c0]])
        mm1 = np.array([[c0,c0],[c0,c0]])
        mm2 = np.array([[c0,c0],[c0,c0]])

        self.cox[0] = 2.65e-6
        self.cox[1] = 0.15e-6
        self.cox[2] = -0.15e-6
        #self.cox[3] = -4.15e-6
        
        for i in range(self.nlay):
            self.kzi[i] = cmath.sqrt((self.k0*self.nindx[i])**2 - kz**2)
            if self.kzi[i].imag < 0:
                self.kzi[i] = self.kzi[i]
            else:
                self.kzi[i] = -self.kzi[i]
            
        tm0[0][0] = c0 #cmath.exp(cj*self.kzi[0]*self.cox[0])
        tm0[0][1] = cmath.exp(-cj*self.kzi[0]*self.cox[1])
        tm0[1][0] = c0 #cmath.exp(cj*self.kzi[0]*self.cox[0])* self.kzi[0]
        tm0[1][1] = -cmath.exp(-cj*self.kzi[0]*self.cox[1])* self.kzi[0]

        mm1[0][0] = cmath.exp(cj*self.kzi[1]*self.cox[1])
        mm1[0][1] = cmath.exp(-cj*self.kzi[1]*self.cox[1])
        mm1[1][0] = cmath.exp(cj*self.kzi[1]*self.cox[1])*self.kzi[1]
        mm1[1][1] = -cmath.exp(-cj*self.kzi[1]*self.cox[1])*self.kzi[1]
        
        tm1[0][0] = cmath.exp(cj*self.kzi[1]*self.cox[2])
        tm1[0][1] = cmath.exp(-cj*self.kzi[1]*self.cox[2])
        tm1[1][0] = cmath.exp(cj*self.kzi[1]*self.cox[2])*self.kzi[1]
        tm1[1][1] = -cmath.exp(-cj*self.kzi[1]*self.cox[2])*self.kzi[1]
        
        mm2[0][0] = cmath.exp(cj*self.kzi[2]*self.cox[2])
        mm2[0][1] = c0 #cmath.exp(-cj*self.kzi[2]*self.cox[1])
        mm2[1][0] = cmath.exp(cj*self.kzi[2]*self.cox[2])*self.kzi[2]
        mm2[1][1] = c0 #-cmath.exp(cj*self.kzi[2]*self.cox[1])*self.kzi[2]

        sl = np.dot(np.linalg.inv(tm1), mm2)
        su = np.dot(np.linalg.inv(mm1), tm0)

        
        re = sl[0][0]*su[1][1]-sl[1][0]*su[0][1]

        return re

            

dg = Dwg()
dg.init()
dg.initgus()
        
    
