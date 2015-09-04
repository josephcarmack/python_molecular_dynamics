# code for simulating dipole-dipole colloidal particle interactions in the pres-
# ence of a constant external electric field.

from numpy import *
import functions as f

# parameters
dim = 2                             # number of spacial dimensions
N = 2                               # number of colloidal particles
pos = zeros((N,dim))                # array for storing particle positions
vel = zeros((N,dim))                # array for storing particle velocities
force = zeros((N,dim))              # array to store particle forces
efield = zeros(dim)                 # array to store efield vector
eMag = 10.0                         # magnitude of the efield
rad = 5.0                           # particle radius

# integration scheme parameters
dt = 0.01                           # time step
steps = 1000                        # number of steps to integrate
total_time = dt*steps               # total simulation time

# particle interaction parameters
FR = 1.0                            # repulsive force constant

# simulation domain
L = 100.0                           # box length
box = zeros((dim))                  # initialize box
for dimension in box:
    dimension = L                   # set box dimensions

# initialize system assuming electric field is on and particles instantly polar-
# ize and their dipoles rotate to align with the electric field direction

efield[0] = eMag                    # set direction and magnitude of efield
pos = f.pos_init(pos,dim,L)         # initialize positions
vel = f.vel_init(vel,dim,2.0)       # initialize velocities


print "this is the end of the script!"
