import numpy as np
import math


class Gust:
    def __init__(self, ude, nchords, mac, dt, tas):
        self.ude = ude
        self.time = None
        self.tas = None
        self.s = None
        self.u = None
        self.nchords = nchords
        self.mac = mac
        self.np = 0
        self.dt = dt
        self.G = self.nchords * self.mac
        self.dudt = None
        self.duds = None

        self.init_time()
        self.init_tas(tas)
        self.init_gust()

    def init_time(self):
        self.time = []
        self.time.append(0)
        self.s = []
        self.s.append(0)

    def init_tas(self, v0):
        self.tas = []
        self.tas.append(v0)

    def init_gust(self):
        self.u = []
        self.u.append(0)

    def evaluate_1_cosine_u_gust(self, tas):
        self.tas.append(tas)
        self.np = self.np + 1
        self.time.append(self.time[self.np-1]+self.dt)
        mean_tas = 0.5*(float(self.tas[self.np]) + float(self.tas[self.np-1]))
        self.s.append(float(self.s[self.np-1]) + self.dt*mean_tas/self.mac)
        self.u.append(self.ude * 0.5 * (1-math.cos(2*math.pi*self.s[self.np]/self.G)))

    def evaluate_u_gust(self, t):
        x = np.asarray(self.time)
        y = np.asarray(self.u)
        return np.interp(t, x, y)




    
