import numpy as np


class Rk4:
    def __init__(self, nvar, y0):
        self.nvar = int(nvar) 
        self.ya = np. empty(self.nvar)
        self.k = np.empty((4, self.nvar))
        self.y = np.empty((4, self.nvar))               
        self.y0 = np.empty(self.nvar)
        
        for i in range(self.nvar):
            self.yo[i] = float(y0[i])
        

    def rk(self, rkstep, dt, time, dydt):
        #
        # RK 4th Order for an ODE
        # F = integral [ f(t) dt ]
        # F' = f(t)
        #
        for i in range(self.nvar):
            if rkstep == 1:
                if time == 0:
                    self.ya[i] = self.y0[i]
                else:
                    self.ya[i] = self.y[3]
                    
                self.k[rkstep - 1, i] = dydt[i]
                self.y[rkstep - 1, i] = self.ya[i] + 0.5 * dt * self.k[rkstep - 1, i]
                
            elif rkstep == 2:
                self.k[rkstep - 1, i] = dydt[i]
                self.y[rkstep - 1, i] = self.ya[i] + 0.5 * dt * self.k[rkstep - 1, i]
                
            elif rkstep == 3:
                self.k[rkstep - 1, i] = dydt[i]
                self.y[rkstep - 1, i] = self.ya[i] + dt * self.k[rkstep - 1, i]
                
            elif rkstep == 4:
                self.k[rkstep - 1, i] = dydt[i]
                self.y[rkstep - 1, i] = self.ya[i] + (dt / 6) * (self.k[0, i] + 2 * self.k[1, i] + 2 * self.k[2, i] + self.k[3, i])
