# Wagner function: the transient lift response due unit change in angle of attack'
# Source: Lomax - Structural Loads Analysis'

import numpy as np
import math

class Wagner:
    
    def __init__(self):
        self.b = np.empty(4)
        self.beta = np.empty(3)


    def define_coefficients(self):
        self.b[0] = 1.0
        self.b[1] = -0.165
        self.b[2] = -0.335
        self.b[3] = 0.0
        self.beta[0] = 0.09
        self.beta[1] = 0.6
        self.beta[2] = 0.0


    def evaluate_wagner_step_function(self, time, tas, mac):
        return self.b[0] + self.b[1] * math.exp(-1.0 * self.beta[0] * (time * tas) / mac ) \
                         + self.b[2] * math.exp(-1.0 * self.beta[1] * (time * tas) / mac ) \
                         + self.b[3] * math.exp(-1.0 * self.beta[2] * (time * tas) / mac )
    

    def evaluate_wagner_impulsive_function(self, s, time, tas, mac):
        v = []
        v.append(-1.0 * self.b[1] * self.beta[0] * (tas / mac) * (math.exp(-1.0 * self.beta[0] * s)) * (math.exp(1.0 * self.beta[0] * ((time * tas) / mac))))
        v.append(-1.0 * self.b[2] * self.beta[1] * (tas / mac) * (math.exp(-1.0 * self.beta[1] * s)) * (math.exp(1.0 * self.beta[1] * ((time * tas) / mac))))
        v.append(-1.0 * self.b[3] * self.beta[2] * (tas / mac) * (math.exp(-1.0 * self.beta[2] * s)) * (math.exp(1.0 * self.beta[2] * ((time * tas) / mac))))
        
        return v