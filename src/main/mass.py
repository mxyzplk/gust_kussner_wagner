import numpy as np
import math

class Mass:
    def __init__(self):
        self.weight    = None
        self.Ixx       = None
        self.Iyy       = None
        self.Izz       = None
        self.Ixy       = None
        self.Iyz       = None
        self.Ixz       = None
        self.x         = None
        self.y         = None
        self.z         = None
        self.A         = None
        self.C         = None
        self.D         = None
        self.E         = None
        self.F         = None


    def def_inertia(self, weight, ixx, iyy, izz, ixy, iyz, ixz):
        self.weight    = weight
        self.Ixx       = ixx
        self.Iyy       = iyy
        self.Izz       = izz
        self.Ixy       = ixy
        self.Iyz       = iyz
        self.Ixz       = ixz

        self.A         = (( self.Iyy-self.Izz ) / self.Ixx) - (( self.Ixz * self.Ixz)/ (self.Ixx * self.Izz))
        self.C         = ( (self.Ixz * self.Ixz ) / ( self.Ixx * self.Izz )) + ( self.Ixx - self.Iyy ) / self.Izz
        self.D         = ( self.Ixz * ( self.Iyy - self.Izz )) / ( self.Ixx * self.Izz ) - ( self.Ixz / self.Izz)
        self.F         = ( self.Ixz * (self.Ixx - self.Iyy ))/(self.Ixx * self.Izz) + (self.Ixz / self.Ixx)
        self.E         = 1.0 / (1.0 - ((self.Ixz * self.Ixz) / (self.Ixx * self.Izz ))) 


    def def_cg(self, x, y ,z):
        self.x         = x
        self.y         = y
        self.z         = z

    
