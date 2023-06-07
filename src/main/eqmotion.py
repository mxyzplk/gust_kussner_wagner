import numpy as np
import math
import kussner
import wagner
import gust
import mass
import aircraft
import airspeed
import rk


class Dof1:
    def __init__(self):
        # 1 DOF Equation:
        # wp(s) = (q*Sw*dCLda/weight)ae(s)
        
        self.v0 = None
        self.h0 = None
        self.aoa0 = None
        self.ae = None
        self.ag = None
        self.ac = None
        self.gust_obj = None
        self.kussner_obj = kussner.Kussner()
        self.wagner_obj = wagner.Wagner()
        self.mass_obj = mass.Mass()
        self.ac_obj = aircraft.Aircraft()
        self.airspeed_obj = airspeed.Airspeed()

    def set_gust(self, ude, nchords, mac, dt):
        self.gust_obj = gust.Gust(ude, nchords, mac, dt)

    def initial_conditions_dof1(self, v0, h0, aoa0, weight):
        self.v0 = v0
        self.h0 = h0
        self.aoa0 = aoa0
        self.mass_obj.def_inertia(weight)

    def init_a_vectors(self):
        self.ae = []
        self.ag = []
        self.ac = []
        self.ae.append(0)
        self.ag.append(0)
        self.ac.append(0)

    def eq_motion_dof1(self):
        # Results file
        f = open("dof1.txt", "w")

        rkdof1 = rk.Rk4(6)
        
        kussner_impulsive = np.empty(3)
        #
        # 0: kussner derivative function part 1
        # 1: kussner derivative function part 2
        # 2: kussner derivative function part 3
        # 3: wagner derivative function part 1
        # 4: wagner derivative function part 2
        # 5: vertical acceleratiom
        #
        rkdof1.dydt[:] = 0  # derivatives values in t = 0
        rkdof1.y0[:] = 0 # initial condition
    
        while(True):
                        
            for rkstep in range(4):
                if rkstep == 0: # RK step 1
                    time = self.gust_obj.time
                elif rkstep == 1: # RK step 2
                    time = self.gust_obj.time + 0.5 * self.gust_obj.dt
                elif rkstep == 2: # RK step 3
                    time = self.gust_obj.time + 0.5 * self.gust_obj.dt
                elif rkstep == 3: # RK step 4
                    time = self.gust_obj.time + self.gust_obj.dt
                
                # evaluate a_gust(t)
                s = time * self.v0 / self.gust_obj.mac
                u = self.gust_obj.evaluate_1_cosine_u_gust(s)
                a_gust = math.atan(u / self.v0)
                
                # evaluate kussner impulsive function
                kussner_impulsive = self.kussner_obj.kussner_impulsive_func(s, time, self.v0, self.gust_obj.mac)
                
                # evaluate kussner step function in s = 0:
                kussner_zero = self.kussner_obj.kussner_step_func(0)
                
                rkdof1.rk(rkstep, self.gust_obj.dt, time)
                
                
                
                # Next time step
                if rkstep == 3: 
                    self.gust_obj.time = self.gust_obj.time + self.gust_obj.dt
                    self.gust_obj.increment_arrays(self.v0, u)
                    



        