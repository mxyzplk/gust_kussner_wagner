import math


class Gust:
    def __init__(self, ude, nchords, mac, dt):
        self.ude = ude
        self.time = 0
        self.tas = None
        self.s = None
        self.u = None
        self.np = None
        self.nchords = nchords
        self.mac = mac
        self.dt = dt
        self.G = self.nchords * self.mac
        
    def init_arrays(self):
        self.s = []
        self.s.append(0)
        self.u = []
        self.u.append(0)      
        self.np = 0        
        
    def increment_arrays(self, tas, u):
        self.s.append(self.s[self.np] + (tas * self.dt / self.mac))
        self.time.append(self.time[self.np] + self.dt)
        self.u.append(u)
        self.np = self.np + 1      

    def evaluate_1_cosine_u_gust(self, s):
        return self.ude * 0.5 * (1 - math.cos(2 * math.pi * s / self.G))




    
