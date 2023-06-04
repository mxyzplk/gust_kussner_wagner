# Küsser function: the transient lift response due to a sharp-edge gust
# Source: Lomax - Structural Loads Analysis

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

    def evaluate_kussner_function(self, s):
        return self.b[0] + self.b[1] * math.exp(-1.0 * self.beta[0] * s) \
            + self.b[2] * math.exp(-1.0 * self.beta[1] * s) \
            + self.b[3] * math.exp(-1.0 * self.beta[2] * s)

    def evaluate_ag_lomax(self, s, g):
        # Lomax - Table 5.6
        c = np.empty(3)
        k = g / (2 * math.pi)

        for i in range(3):
            c[i] = 2 * math.pi ** 2 * self.b[i + 1] * (-1.0 * k * self.beta[i] + \
                                                       math.cos(s / k) - math.exp(-1.0 * self.beta[i] * s)) / \
                                                       (g * g * self.beta[i] * self.beta[i] + 4 * math.pi ** math.pi)

        return -0.5 * self.beta[0] * (1.0 - math.cos(s / k)) - c[0] - c[1] - c[2]
