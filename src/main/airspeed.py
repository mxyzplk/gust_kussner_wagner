import math


class Airspeed:
    def __init__(self):
        self.a0 = 340.294     # sonic speed at sea level ISA
        self.p = None        # static pressura (Pa)
        self.p0 = 101325      # static sea level pressure ISA (Pa)
        self.rho = None        # air density (kg/m3)
        self.rho0 = 1.225       # air density sea level (kg/m3)
        self.rhostd = None
        self.t0 = 288.15      # standard temperature ISA (K)
        self.t = None
        self.tstd = None
        self.a = None        # sonic speed
        self.isa = None        # ISA
        self.dynp = None        # dynamic pressura (Pa)
        self.eas = None        # equivalent airspeed
        self.alt = None        # altitude
        self.mach = None        # mach
        self.cas = None        # calibrated airspeed
        self.tas = None        # true airspeed
        self.R = 287.05287
        self.mu0 = 1.7894*10**(-5)
        self.m0 = 28.96442
        self.R1 = 8314.32
        self.gamma = 1.4
        self.mu = None

    def atmos(self, alt, isa):

        self.alt = alt
        self.isa = isa

        if self.alt <= 11000:
            b = -0.0065
            a = 288.15
            c = 8.9619638
            d = -0.00020216125
            e = 5.2558797
            i = 1.04884
            j = -0.000023659414
            l = 4.2558797
            self.tstd = a + b * self.alt
            self.p = (c + d * self.alt) ** e
            self.rhostd = (i + j * self.alt) ** l
        elif self.alt > 11000 <= 20000:
            b = 0
            a = 216.65
            f = 128244.5
            g = -0.00015768852
            m = 2.06214
            n = -0.00015768852
            self.tstd = a + b * self.alt
            self.p = f * math.exp(g * self.alt)
            self.rhostd = m * math.exp(n ** self.alt)
        elif (self.alt > 20000 <= 32000):
            b = 0.001
            a = 196.65
            c = 0.70551848
            d = 0.0000035876861
            e = -34.163218
            i = 0.9726309
            j = 0.000004946
            l = -35.163218
            self.tstd = a + b * self.alt
            self.p = (c + d * self.alt) ** e
            self.rhostd = (i + j * self.alt) ** l
        elif (self.alt > 32000 <= 47000):
            b = 0.0028
            a = 13900000
            c = 0.34926867
            d = 0.000007033098
            e = -12.201149
            i = 0.84392929
            j = 0.000016993902
            l = -13.201149
            self.tstd = a + b * self.alt
            self.p = (c + d * self.alt) ** e
            self.rhostd = (i + j * self.alt) ** l
        elif (self.alt > 47000 <= 50000):
            b = 0
            a = 270.65
            f = 41828.42
            g = -0.00012622656
            m = 0.53839563
            n = -0.00012622656
            self.tstd = a + b * self.alt
            self.p = f * math.exp(g * self.alt)
            self.rhostd = m * math.exp(n *self.alt)


    def evaluate_airpeed(self, alt, isa, speed, type):
        
        Airspeed.atmos(self, alt, isa)
        
        if (type == "e"):
            self.eas = speed
            self.mach = self.eas/(self.a0*math.sqrt(self.p/self.p0))
            self.tas = self.mach*self.a
            self.cas = math.sqrt(2 * self.a0**2 / (self.gamma - 1) * ((1 + (self.p/self.p0) * ((1 + (self.gamma - 1) / 2 * ((self.eas / self.a0)**2) / (self.p/self.p0))**(self.gamma / (self.gamma - 1)) - 1))**((self.gamma - 1) / self.gamma) - 1))
        elif (type == "c"):
            self.cas = speed
            self.eas = math.sqrt(2 * (self.p/self.p0) * self.a0**2 / (self.gamma - 1) * ((1 + 1 / (self.p/self.p0) * ((1 + (self.gamma - 1) / 2 * (self.cas / self.a0)**2)**(self.gamma / (self.gamma - 1)) - 1))**((self.gamma - 1) / self.gamma) - 1))
            self.mach = self.eas/(self.a0*math.sqrt(self.p/self.p0))
            self.tas = self.mach*self.a
        elif (type == "t"):
            self.tas = speed
            self.mach = self.tas/self.a
            self.eas = self.a0 * self.mach * math.sqrt(self.p/self.p0)
            self.cas = math.sqrt(2 * self.a0**2 / (self.gamma - 1) * ((1 + (self.p/self.p0) * ((1 + (self.gamma - 1) / 2 * ((self.eas / self.a0)**2) / (self.p/self.p0))**(self.gamma / (self.gamma - 1)) - 1))**((self.gamma - 1) / self.gamma) - 1))
        elif (type == "m"):
            self.mach = speed
            self.tas = self.mach*self.a
            self.eas = self.a0 * self.mach * math.sqrt(self.p/self.p0)
            self.cas = math.sqrt(2 * self.a0**2 / (self.gamma - 1) * ((1 + (self.p/self.p0) * ((1 + (self.gamma - 1) / 2 * ((self.eas / self.a0)**2) / (self.p/self.p0))**(self.gamma / (self.gamma - 1)) - 1))**((self.gamma - 1) / self.gamma) - 1))

        self.dynp = 0.5*self.rho0*self.eas**2
        
