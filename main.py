# code for simulating dipole-dipole colloidal particle interactions in the pres-
# ence of a constant external electric field.

from numpy import *
import functions as f
import os

# delete old vtk files from previous runs
if os.path.exists('./vtkoutput'):
    os.system('rm -r vtkoutput')

# parameters
dim = 2                             # number of spacial dimensions
N = 25                               # number of colloidal particles
pos = zeros((N,dim))                # array for storing particle positions
vel = zeros((N,dim))                # array for storing particle velocities
force = zeros((N,dim))              # array to store particle forces
efield = zeros(dim)                 # array to store efield vector
eMag = 10.0                         # magnitude of the efield
rad = 1.0                           # particle radius
density = 3.0/(4.0*pi)*10           # particle density
vol = 4.0/3.0*pi*rad**3             # particle volume
mass = density*vol                  # particle mass
kb = 1.0                            # Boltzmann Constant
temp = 1.0                          # temperature in Kelvin
max_vel = 10.0                      # maximum particle velocity allowed
vel_spread = sqrt(kb*temp/(mass))   # std 4 velocity generator (micron/microsec)

# integration scheme parameters
dt = 0.001                          # time step (mili-seconds)
steps = 100000                      # number of steps to integrate
total_time = dt*steps               # total simulation time
print 'total time = '+str(total_time)+' time units'

# file io parameters
skip = 100                          # how many steps to skip between outputs
description = 'particle'            # a description of the type of vtk data

# particle interaction parameters
eps = 1.0                           # the leonard-jones energy parameter
sigma = 1.0                         # the leonard-jones length parameter

# simulation domain
L = 50.0                             # box length (microns)
box = zeros(dim)                    # initialize box
ao = 4.0*rad                        # lattice constant
for i in range(len(box)):
    box[i] = L                      # set box dimensions

# initialize system assuming electric field is on and particles instantly polar-
# ize and their dipoles rotate to align with the electric field direction

efield[0] = eMag                    # set direction and magnitude of efield

# randomely initialize positions
# pos = f.pos_init_rand(pos,dim,box,rad)

# initialize pos on lattice w/noise
pos,box,eflag = f.pos_init_lat(pos,dim,box,rad,ao,N)
if eflag:
    quit()

# initialize velocities
vel = f.vel_init(vel,dim,vel_spread)

# write initial particle positions to vtk file
f.write_point_data_vtk(description,pos,N,dim,1)

# EVOLVE SYSTEM

for step in range(steps):

    # zero out forces
    force = force*0.0

    for i in range(N):
        for j in range(i+1,N):

            # get positions of particle i and j
            ri = pos[i]
            rj = pos[j]

            # get separation distance and unit vector that points in the direct-
            # ion of the force on particle i for the i-j particle pair.
            rij_mag,rij_unit = f.calc_rij_pbc(ri,rj,dim,box)

            # get force on particle i from particle j
            fmag = f.calc_repulsion_force(rij_mag,rad,sigma,eps)
            force[i] = force[i] + fmag*rij_unit

            # applying Newton's third law to get force on particle j from part-
            # icle j
            force[j] = force[j] - fmag*rij_unit

            if rij_mag < rad:
                print 'particle separation distance less than 2.1*rad',
                print ' for part pair '+str(i)+','+str(j)
                print 'part i pos = '+str(ri)+', part j pos = '+str(rj)
                print 'rij_unit = '+str(rij_unit)
                print 'step = '+str(step)
                quit()

    # update the velocities and then positions based on the calculated forces
    for i in range(N):
        # get the net acceleration on particle i from its net force
        ai = 1/mass*force[i]
        # update the velocity of particle i using its acceleration
        vel[i] = vel[i] + dt*ai
        # enforce maximum particle velocity
        vel_mag = linalg.norm(vel[i])
        if vel_mag > max_vel:
            vel[i] = max_vel
        # update the position of particle i using its updated velocity
        pos[i] = pos[i] + dt*vel[i]
        # implement periodic boundary conditions
        f.enforce_pbc(pos[i],dim,box)

    # write output ever skip steps
    if step > 0:
        if (step % skip) == 0:
            f.write_point_data_vtk(description,pos,N,dim,step)

# write final step output
f.write_point_data_vtk(description,pos,N,dim,steps)
