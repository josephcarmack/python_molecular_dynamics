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
rad = 1.0                           # particle radius (microns)
density = 5.0                       # particle density (micro-grams/micron^3)
vol = 4.0/3.0*pi*rad**3             # particle volume (microns^3)
mass = density*vol                  # particle mass (micro-grams)
vel_spread = 2.0                    # std 4 velocity generator (micron/microsec)

# integration scheme parameters
dt = 0.01                           # time step (micro-seconds)
steps = 1000                        # number of steps to integrate
total_time = dt*steps               # total simulation time

# particle interaction parameters
sigma = 1.0                         # the leonard-jones parameter particle-part-
                                    # icle interactions.

# simulation domain
L = 100.0                           # box length (microns)
box = zeros(dim)                    # initialize box
for i in range(len(box)):
    box[i] = L                      # set box dimensions

# initialize system assuming electric field is on and particles instantly polar-
# ize and their dipoles rotate to align with the electric field direction

efield[0] = eMag                    # set direction and magnitude of efield
# pos = f.pos_init(pos,dim,L)       # initialize positions
pos[0] = array([50.0,45.0])
pos[1] = array([50.0,60.0])
vel = f.vel_init(vel,dim,vel_spread)# initialize velocities

# EVOLVE SYSTEM

# for step in range(steps):

    # zero out forces

for i in range(N):
    for j in range(i+1,N):

        # get positions of particle i and j
        ri = pos[i]
        rj = pos[j]

        # get separation distance and unit vector that points in the direct-
        # ion of the force on particle i for the i-j particle pair.
        rij_mag,rij_unit = f.calc_rij_pbc(ri,rj,dim,box)

        # get force on particle i from particle j
        fmag = f.calc_repulsion_force(rij_mag,rad,sigma)
        force[i] = force[i] + fmag*rij_unit

        # applying Newton's third law to get force on particle j from part-
        # icle j
        force[j] = force[j] - fmag*rij_unit

# update the velocities and then positions based on the calculated forces
for i in range(N):
    # get the net acceleration on particle i from its net force
    ai = 1/mass*force[i]
    # update the velocity of particle i using its acceleration
    vel[i] = vel[i] + dt*ai
    # update the position of particle i using its updated velocity
    pos[i] = pos[i] + dt*vel[i]

print pos

print "this is the end of the script!"
