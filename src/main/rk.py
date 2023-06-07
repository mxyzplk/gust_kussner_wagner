import numpy as np


class Rk4:
    def __init__(self, nvars):
        self.y0 = np.empty(int(nvars))
        self.y = np.empty(int(nvars))
        self.ya = np.empty(int(nvars))
        self.dydt = np.empty(int(nvars), 4)
        self.nvars = int(nvars)
        self.k = np.empty(int(nvars), 4)


    def rk(self, rkstep, tstep, time):
        #
        # RK 4th Order for an ODE
        # F = integral [ f(t) dt ]
        # F' = f(t)
        #
        for rkstep in range(self.nvars):
            if rkstep == 0:
                if time == 0:
                    self.ya[rkstep] = self.y0[rkstep]
                else:
                    self.ya[rkstep] = self.y[rkstep]

                self.k[rkstep, rkstep] = self.dydt[rkstep]
                self.y[rkstep] = self.ya[rkstep] + 0.5 * tstep * self.k[rkstep, 0]
            elif rkstep == 1:
                self.k[rkstep, rkstep] = self.dydt[rkstep]
                self.y[rkstep] = self.ya[rkstep] + 0.5 * tstep * self.k[rkstep, 1]
            elif rkstep == 2:
                self.k[rkstep, rkstep] = self.dydt[rkstep]
                self.y[rkstep] = self.ya[rkstep] + 0.5 * tstep * self.k[rkstep, 2]
            elif rkstep == 3:
                self.k[rkstep, rkstep] =self.dydt[rkstep]
                self.y[rkstep] = self.ya[rkstep] + (tstep / 6) * (self.k[rkstep, 0] + 2 * self.k[rkstep, 1] + 2 * self.k[rkstep, 2] + self.k[rkstep, 3])
