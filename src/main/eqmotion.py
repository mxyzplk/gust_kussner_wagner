import numpy as np
import math
import kussner
import wagner
import gust
import mass
import aircraft
import airspeed

class Dof1:
    def __init__(self, Ude, nchords, mac, dt):
        # 1 DOF Equation:
        # wp(s) = (q*Sw*dCLda/weight)ae(s)
        
        self.v0            = None
        self.h0            = None
        self.ae            = None
        self.ag            = None
        self.ac            = None
        self.gust_obj      = gust.Gust(Ude, nchords, mac, dt)
        self.kussner_obj   = kussner.Kussner()
        self.wagner_obj    = wagner.Wagner()
        self.mass_obj      = mass.Mass()
        self.ac_obj        = aircraft.Aircraft()
        self.airspeed_obj  = airspeed.Airspeed()



    def initial_conditions(self, v0, h0):
        self.v0 = v0
        self.h0 = h0


    def init_a_vectors(self):
        self.ae = []
        self.ag = []
        self.ac = []
        self.ae.append(0)
        self.ag.append(0)
        self.ac.append(0)


    def eq_motion_dof1(self):
        