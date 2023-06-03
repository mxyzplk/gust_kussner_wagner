import numpy as np
import math

class Airspeed:
    def __init__(self):
        self.a0       = 340.294     # sonic speed at sea level ISA 
        self.p        = None        # static pressura (Pa)
        self.p0       = 101325      # static sea level pressure ISA (Pa)
        self.rho      = None        # air density (kg/m3)
        self.rho0     = 1.225       # air density sea level (kg/m3)
        self.rhostd   = None
        self.t0       = 288.15      # standard temperature ISA (K)
        self.t        = None
        self.tstd     = None
        self.a        = None        # sonic speed
        self.isa      = None        # ISA
        self.dynp     = None        # dynamic pressura (Pa)
        self.eas      = None        # equivalent airspeed
        self.alt      = None        # altitude
        self.mach     = None        # mach
        self.cas      = None        # calibrated airspeed
        self.tas      = None        # true airspeed
        self.R        = 287.05287
        self.mu0      = 1.7894*10**(-5)
        self.m0       = 28.96442
        self.R1       = 8314.32
        self.gamma    = 1.4
        self.mu       = None

    def atmos(self, alt, isa):

        self.alt = alt
        self.isa = isa

        if (self.alt <= 11000):
            B = -0.0065
            A = 288.15
            C = 8.9619638
            D = -0.00020216125
            E = 5.2558797
            I = 1.04884
            J = -0.000023659414
            L = 4.2558797
            self.tstd = A + B * self.alt
            self.p = (C + D * self.alt) ** E
            self.rhostd = (I + J * self.alt) ** L
        elif (self.alt > 11000 <= 20000):
            B = 0
            A = 216.65
            F = 128244.5
            G = -0.00015768852
            M = 2.06214
            N = -0.00015768852
            self.tstd = A + B * self.alt
            self.p = F * math.exp(G * self.alt)
            self.rhostd = M * math.exp(N ** self.alt)
        elif (self.alt > 20000 <= 32000):
            B = 0.001
            A = 196.65
            C = 0.70551848
            D = 0.0000035876861
            E = -34.163218
            I = 0.9726309
            J = 0.000004946
            L = -35.163218
            self.tstd = A + B * self.alt
            self.pi = (C + D * self.alt) ** E
            self.rhostd = (I + J * self.alt) ** L
        elif (self.alt > 32000 <= 47000):
            B = 0.0028
            A = 13900000
            C = 0.34926867
            D = 0.000007033098
            E = -12.201149
            I = 0.84392929
            J = 0.000016993902
            L = -13.201149
            self.tstd = A + B * self.alt
            self.p = (C + D * self.alt) ** E
            self.rhostd = (I + J * self.alt) ** L
        elif (self.alt > 47000 <= 50000):
            B = 0
            A = 270.65
            F = 41828.42
            G = -0.00012622656
            M = 0.53839563
            N = -0.00012622656
            self.tstd = A + B * self.alt
            self.p = F * math.exp(G * self.alt)
            self.rhostd = M * math.exp(N *self.alt)

        self.t = self.tstd + self.isa
        self.rho = self.rhostd / (1.0 + self.isa/self.tstd)
        self.a = math.sqrt(self.gamma*self.R*self.t)
        self.mu = self.mu0 * ((self.t/self.t0)**1.5)*(self.t0+120)/(self.t+120)


    def evaluate_airpeed(self, alt, isa, speed, type):
        
        Airspeed.atmos(alt, isa)
        
        if (type == "e"):
            self.eas = speed
            self.mach = self.eas/(self.a0*math.sqrt(self.p/self.p0))
            self.tas = self.mach*self.a
            self.cas = math.sqr(2 * self.a0**2 / (self.gamma - 1) * ((1 + (self.p/self.p0) * ((1 + (self.gamma - 1) / 2 * ((self.eas / self.a0)**2) / (self.p/self.p0))**(self.gamma / (self.gamma - 1)) - 1))**((self.gamma - 1) / self.gamma) - 1))
        elif (type == "c"):
            self.cas = speed
            self.eas = math.sqr(2 * (self.p/self.p0) * self.a0**2 / (self.gamma - 1) * ((1 + 1 / (self.p/self.p0) * ((1 + (self.gamma - 1) / 2 * (self.cas / self.a0)**2)**(self.gamma / (self.gamma - 1)) - 1))**((self.gamma - 1) / self.gamma) - 1))
            self.mach = self.eas/(self.a0*math.sqrt(self.p/self.p0))
            self.tas = self.mach*self.a
        elif (type == "t"):
            self.tas = speed
            self.mach = self.tas/self.a
            self.eas = self.a0 * self.mach * math.sqr(self.p/self.p0)
            self.cas = math.sqr(2 * self.a0**2 / (self.gamma - 1) * ((1 + (self.p/self.p0) * ((1 + (self.gamma - 1) / 2 * ((self.eas / self.a0)**2) / (self.p/self.p0))**(self.gamma / (self.gamma - 1)) - 1))**((self.gamma - 1) / self.gamma) - 1))
        elif (type == "m"):
            self.mach = speed
            self.tas = self.mach*self.a
            self.eas = self.a0 * self.mach * math.sqr(self.p/self.p0)
            self.cas = math.sqr(2 * self.a0**2 / (self.gamma - 1) * ((1 + (self.p/self.p0) * ((1 + (self.gamma - 1) / 2 * ((self.eas / self.a0)**2) / (self.p/self.p0))**(self.gamma / (self.gamma - 1)) - 1))**((self.gamma - 1) / self.gamma) - 1))

        self.dynp = 0.5*self.rho0*self.eas**2
        
