# code for simulating dipole-dipole colloidal particle interactions in the pres-
# ence of a constant external electric field.

from numpy import *
import functions as f
import os
import time as t

# setup timer for profiling
time_lapse = t.time()

# delete old vtk files from previous runs
if os.path.exists('./vtkoutput'):
    os.system('rm -r vtkoutput')

# parameters
dim = 2                             # number of spacial dimensions
N = 9                              # number of colloidal particles
pos = zeros((N,dim))                # array for storing particle positions
vel = zeros((N,dim))                # array for storing particle velocities
force = zeros((N,dim))              # array to store particle forces
efield = zeros(dim)                 # array to store efield vector
eMag = 1.0                          # magnitude of the efield
rad = 1.0                           # particle radius
density = 3.0/(4.0*pi)*10           # particle density
vol = 4.0/3.0*pi*rad**3             # particle volume
mass = density*vol                  # particle mass
kb = 1.0                            # Boltzmann Constant
temp = 1.0                          # temperature in Kelvin
max_vel = 10.0                      # maximum particle velocity allowed
vel_spread = sqrt(kb*temp/(mass))   # std 4 velocity generator (micron/microsec)
drag_coef = 3.0                     # viscosity coefficient
brownian_str = 15.0                 # strength of brownian motion

# integration scheme parameters
dt = 0.001                          # time step (mili-seconds)
steps = 100000                       # number of steps to integrate
total_time = dt*steps               # total simulation time
print 'total time = '+str(total_time)+' time units'

# file io parameters
skip = 100                          # how many steps to skip between outputs
description = 'particle'            # a description of the type of vtk data

# particle interaction parameters
eps = 1.0                           # the leonard-jones energy parameter
sigma = 1.0                         # the leonard-jones length parameter
rcut = 7*rad                        # cut off radius

# simulation domain
L = 15.0                            # box length (microns)
box = zeros(dim)                    # initialize box
ao = 4.0*rad                        # lattice constant
for i in range(len(box)):
    box[i] = L                      # set box dimensions

# initialize system assuming electric field is on and particles instantly polar-
# ize and their dipoles rotate to align with the electric field direction

efield[dim-1] = eMag                # set direction and magnitude of efield

# dipole interaction testing intitialization
#pos,box = f.pos_init_dipole_test(pos,box)

# randomely initialize positions
pos = f.pos_init_rand(pos,dim,box,rad)

# initialize pos on lattice w/noise
#pos,box,eflag = f.pos_init_lat(pos,dim,box,rad,ao,N)
#if eflag:
#    quit()

# initialize velocities
#vel = f.vel_init(vel,dim,vel_spread)

# write initial particle positions to vtk file
f.write_point_data_vtk(description,pos,N,dim,1)

# print some timing information
print "initialization time: %s seconds" % (str(t.time()-time_lapse))
time_lapse = t.time()
force_time = 0
vel_pos_update_time = 0
write_output_time = 0

# EVOLVE SYSTEM
for step in range(steps):

    # setup timer 
    t1 = t.time()

    # zero out forces
    force = force*0.0
    # calculate pair-wise forces
    for i in range(N):
        for j in range(i+1,N):

            # get positions of particle i and j
            ri = pos[i]
            rj = pos[j]

            # get separation distance and unit vector that points in the direct-
            # ion of the force on particle i for the i-j particle pair.
            rij_mag,rij_unit = f.calc_rij_pbc(ri,rj,dim,box)

            if (rij_mag < rcut):
                # get force on particle i from particle j
                fmag = f.calc_repulsion_force(rij_mag,rad,sigma,eps)
                force[i] = force[i] + fmag*rij_unit

                # applying Newton's third law to get force on particle j from part-
                # icle j
                force[j] = force[j] - fmag*rij_unit

                # add in dipole force
#                fdip = f.calc_dipole_force(rij_unit,rij_mag,rad,efield,eMag)
#                force[i] = force[i] + fdip
#                force[j] = force[j] - fdip

            if rij_mag < rad:
                print 'particle separation distance less than 2.1*rad',
                print 'for part pair '+str(i)+','+str(j)
                print 'part i pos = '+str(ri)+', part j pos = '+str(rj)
                print 'rij_unit = '+str(rij_unit)
                print 'step = '+str(step)
                quit()
    
    # calculate drag forces and brownian forces
    for i in range(N):
        force[i] = force[i] - drag_coef*vel[i]
        force[i] = force[i] + brownian_str*random.normal(0,1,dim)

    # calculate time spent on force calculation
    t2 = t.time()
    force_time = force_time + t2-t1

    # update the velocities and then positions based on the calculated forces
    for i in range(N):
        # get the net acceleration on particle i from its net force
        ai = 1/mass*force[i]
        # update the velocity of particle i using its acceleration
        vel[i] = vel[i] + dt*ai
        # enforce maximum particle velocity
#        vel_mag = linalg.norm(vel[i])
#        if vel_mag > max_vel:
#            vel[i] = max_vel
        # update the position of particle i using its updated velocity
        pos[i] = pos[i] + dt*vel[i]
        # implement periodic boundary conditions
        f.enforce_pbc(pos[i],dim,box)

    # calculate time spent on pos and vel update
    t3 = t.time()
    vel_pos_update_time = vel_pos_update_time + t3-t2

    # write output every 'skip' steps
    if step > 0:
        if (step % skip) == 0:
            f.write_point_data_vtk(description,pos,N,dim,step)

    # calculate time spent on writing output
    t4 = t.time()
    write_output_time = write_output_time + t4-t3

# write final step output
f.write_point_data_vtk(description,pos,N,dim,steps)

# calculate time spent simulating
print "simulation time: %s seconds" % (str(t.time()-time_lapse))
print "time calculating forces: %s seconds" % (str(force_time))
print "time updating vel/pos: %s seconds" % (str(vel_pos_update_time))
print "time writing output: %s seconds" % (str(write_output_time))
