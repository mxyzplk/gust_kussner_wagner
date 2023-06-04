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


    def evaluate_wagner_function(self, s):
        return self.b[0] + self.b[1] * math.exp(-1.0 * self.beta[0] * s) \
                         + self.b[2] * math.exp(-1.0 * self.beta[1] * s) \
                         + self.b[3] * math.exp(-1.0 * self.beta[2] * s)
    
