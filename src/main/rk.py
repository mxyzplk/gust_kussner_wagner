import numpy as np


class Rk4:
    def __init__(self, nvars):
        self.y0 = np.empty(int(nvars))
        self.y = np.empty(int(nvars))
        self.ya = np.empty(int(nvars))
        self.dydt = np.empty(int(nvars), 4)
        self.nvars = int(nvars)
        self.k = np.empty(int(nvars), 4)

    def define_initial_condition(self, funcs):
        for i in range(len(funcs)):
            self.y0[i] = funcs[i]

    def rk(self, rkstep, dydt, tstep, time):
        # RK 4th Order for an ODE
        # F = integral [ f(t) dt ]
        # F' = f(t)
        for i in range(self.nvars):
            if rkstep == 1:
                if time == 0:
                    self.ya[i] = self.y0[i]
                else
                    self.ya[i] = self.y[i]

                self.k[i, 1] = dydt[i]
                self.y[i] = self.ya[i] + 0.5 * tstep * self.k[i, 1]
            elif rkstep == 2:
                self.k[i, 2] = dydt[i]
                self.y[i] = self.ya[i] + 0.5 * tstep * self.k[i, 2]
            elif rkstep == 3:
                self.k[i, 3] = dydt[i]
                self.y[i] = self.ya[i] + 0.5 * tstep * self.k[i, 3]
            elif rkstep == 4:
                self.k[i, 4] = dydt[i]
                self.y[i] = self.ya[i] + (tstep / 6) * (self.k[i, 1] + 2 * self.k[i, 2] + 2 * self.k[i, 3] + self.k[i, 4])
