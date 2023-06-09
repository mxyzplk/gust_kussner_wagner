# Küsser function: the transient lift response due to a sharp-edge gust
# Source: Lomax - Structural Loads Analysis
#         Wright & Cooper - Introduction to Aircraft Aeroelasticity and Loads

import numpy as np
import math


class Kussner:

    def __init__(self):
        self.b = np.empty(4)
        self.beta = np.empty(3)

    def define_coefficients(self):
        self.b[0] = 1.0
        self.b[1] = -0.236
        self.b[2] = -0.513
        self.b[3] = -0.171
        self.beta[0] = 0.116
        self.beta[1] = 0.728
        self.beta[2] = 4.84

    def evaluate_kussner_step_function(self, time, tas, mac):
        return self.b[0] + self.b[1] * math.exp(-1.0 * self.beta[0] * (time * tas) / mac ) \
                         + self.b[2] * math.exp(-1.0 * self.beta[1] * (time * tas) / mac ) \
                         + self.b[3] * math.exp(-1.0 * self.beta[2] * (time * tas) / mac )   


    def evaluate_kussner_impulsive_function(self, s, time, tas, mac):
        v = []
        v.append(-1.0 * self.b[1] * self.beta[0] * (tas / mac) * (math.exp(-1.0 * self.beta[0] * s)) * (math.exp(1.0 * self.beta[0] * ((time * tas) / mac))))
        v.append(-1.0 * self.b[2] * self.beta[1] * (tas / mac) * (math.exp(-1.0 * self.beta[1] * s)) * (math.exp(1.0 * self.beta[1] * ((time * tas) / mac))))
        v.append(-1.0 * self.b[3] * self.beta[2] * (tas / mac) * (math.exp(-1.0 * self.beta[2] * s)) * (math.exp(1.0 * self.beta[2] * ((time * tas) / mac))))
        
        return v
        
    def evaluate_ag_lomax(self, s, g):
        # Lomax - Table 5.6
        c = np.empty(3)
        k = g / (2 * math.pi)

        for i in range(3):
            c[i] = 2 * math.pi ** 2 * self.b[i + 1] * (-1.0 * k * self.beta[i] + \
                                                       math.cos(s / k) - math.exp(-1.0 * self.beta[i] * s)) / \
                                                       (g * g * self.beta[i] * self.beta[i] + 4 * math.pi ** math.pi)

        return -0.5 * self.beta[0] * (1.0 - math.cos(s / k)) - c[0] - c[1] - c[2]



