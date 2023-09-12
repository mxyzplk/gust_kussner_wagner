# Kussner function: the transient lift response due to a 1 - cosine gust
# Source: Hoblit

import rk
import numpy as np
import math


class Kussner:

    def __init__(self, tas, mac):
        self.mac = mac
        # Kussner Function Coefficients
        self.b = np.empty(4)
        self.beta = np.empty(3)
        self.define_coefficients()
        # Runge Kutta
        y0 = [0, 0, 0]        
        self.rk = rk.Rk4(3, y0)
        # Kg(0)
        self.kussner0 = self.kussner_step_function(0, tas)

    def define_coefficients(self):
        # Aspect Ratio: inf | Mach 0
        self.b[0] = 1.0
        self.b[1] = -0.236
        self.b[2] = -0.513
        self.b[3] = -0.171
        self.beta[0] = 0.116
        self.beta[1] = 0.728
        self.beta[2] = 4.84

    def kussner_step_function(self, time, tas):
        return self.b[0] + self.b[1] * math.exp(-1.0 * self.beta[0] * (time * tas) / self.mac ) \
                         + self.b[2] * math.exp(-1.0 * self.beta[1] * (time * tas) / self.mac ) \
                         + self.b[3] * math.exp(-1.0 * self.beta[2] * (time * tas) / self.mac )   


    def kussner_rkstep(self, time, dt, rkstep, angle, tas):
        dydt = []
        #print("[kussner] tas: {:10.4f} time: {:10.4f} angle: {:10.4f} rkstep: {:10.4f}".format(tas, time, angle, rkstep))
        dydt.append(angle * self.beta[0] * (tas / self.mac) * (math.exp(self.beta[0] * ((time * tas) / self.mac))))
        dydt.append(angle * self.beta[1] * (tas / self.mac) * (math.exp(self.beta[1] * ((time * tas) / self.mac))))
        dydt.append(angle * self.beta[2] * (tas / self.mac) * (math.exp(self.beta[2] * ((time * tas) / self.mac))))
        
        self.rk.rk(rkstep, dt, time, dydt)


    def kussner_convolution_integral(self, time, tas, angle, rkstep):                    
        integral = angle * self.kussner0 - self.b[1] * (math.exp(-1.0 * self.beta[0] * (tas / self.mac) * time)) * self.rk.y[rkstep - 1, 0] \
                                         - self.b[2] * (math.exp(-1.0 * self.beta[1] * (tas / self.mac) * time)) * self.rk.y[rkstep - 1, 1] \
                                         - self.b[3] * (math.exp(-1.0 * self.beta[2] * (tas / self.mac) * time)) * self.rk.y[rkstep - 1, 2]
        
        return integral
    
    def eval_rkstep(self,  time, dt, rkstep, angle, tas):
        self.kussner_rkstep(time, dt, rkstep, angle, tas)
        integral =  self.kussner_convolution_integral(time, tas, angle, rkstep)
        
        return integral
    
