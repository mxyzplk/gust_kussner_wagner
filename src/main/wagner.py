# Wagner function: the transient lift response due unit change in angle of attack'
# Source: Hoblit

import rk
import numpy as np
import math

class Wagner:
    
    def __init__(self, tas, mac):
        self.mac = mac
        # Kussner Function Coefficients
        self.b = np.empty(3)
        self.beta = np.empty(2)
        self.define_coefficients()
        # Runge Kutta
        y0 = [0, 0] 
        self.rk = rk.Rk4(2, y0)
        # Kg(0)
        self.wagner0 = self.wagner_step_function(0, tas)


    def define_coefficients(self):
        self.b[0] = 1.0
        self.b[1] = -0.165
        self.b[2] = -0.335
        self.beta[0] = 0.09
        self.beta[1] = 0.6


    def wagner_step_function(self, time, tas):
        return self.b[0] + self.b[1] * math.exp(-1.0 * self.beta[0] * (time * tas) / self.mac ) \
                         + self.b[2] * math.exp(-1.0 * self.beta[1] * (time * tas) / self.mac )
    

    def wagner_rkstep(self, time, dt, rkstep, angle, tas):
        dydt = []
        dydt.append(angle * self.beta[0] * (tas / self.mac) * (math.exp(self.beta[0] * ((time * tas) / self.mac))))
        dydt.append(angle * self.beta[1] * (tas / self.mac) * (math.exp(self.beta[1] * ((time * tas) / self.mac))))
        
        self.rk.rk(rkstep, dt, time, dydt)
    
    
    def wagner_convolution_integral(self, time, tas, angle, rkstep):                    
        integral = angle * self.wagner0 + self.b[1] * (math.exp(-1.0 * self.beta[0] * (tas / self.mac) * time)) * self.rk.y[rkstep - 1, 0] \
                                         + self.b[2] * (math.exp(-1.0 * self.beta[1] * (tas / self.mac) * time)) * self.rk.y[rkstep - 1, 1]
        
        return integral

    def eval_rkstep(self,  time, dt, rkstep, angle, tas):
        self.wagner_rkstep(time, dt, rkstep, angle, tas)
        integral =  self.wagner_convolution_integral(time, tas, angle, rkstep)
        
        return integral