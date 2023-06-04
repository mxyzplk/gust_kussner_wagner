import numpy as np


class Rk4:
    def __init__(self, nvars):
        self.yo = np.empty(int(nvars))
        self.y = np.empty(int(nvars), 4)
        self.dydt = np.empty(int(nvars), 4)
        self.nvars = int(nvars)
        self.k = np.empty(int(nvars), 4)

    def define_initial_condition(self, funcs):
        for i in range(len(funcs)):
            self.yo[i] = funcs[i]

    def rk(self, rkstep, dydt, step):
        # F = integral [ f(t) dt ]
        # F' = f(t)
        for i in range(self.nvars):
            if rkstep == 1:
                self.k[i, 1] = step * dydt[i]
            elif rkstep == 2:
                self.k[i, 2] = step * dydt[i]
            elif rkstep == 3:
                self.k[i, 3] = step * dydt[i]
            elif rkstep == 4:
                self.k[i, 4] = step * dydt[i]
                self.y[i] = self.yo[i] + (1 / 6) * (self.k[i, 1] + 2 * self.k[i, 2] + 2 * self.k[i, 3] + self.k[i, 4])
