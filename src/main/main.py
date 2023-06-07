from eqmotion import Dof1
import math

# Lomax:
# S.I. Units
sw = 1951 * 0.3048 * 0.3048     # wing area
altitude = 20000 * 0.3048       # altitude
weight = 252000 * 0.45359237    # weight
dclda = 3.2397 * 180 / math.pi  # Cl_a [de]
tas = 782.3 * 0.3048            # true airspeed
dt = 0.001                      # Time step
ude = 15.24                     # Design Gust
mac = 199.7 * 0.0254            # mean aerodynamic chord
nchords = 25                    # number of chords

dof1sim = Dof1()

dof1sim.set_gust(ude, nchords, mac, dt, tas)

