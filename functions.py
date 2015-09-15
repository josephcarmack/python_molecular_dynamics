def pos_init_rand(pos,dim,box,rad):
    """
    This function populates the position array with random positions from
    uniformly distributed throughout the box. It also does not allow part-
    cles to be placed too close to each other (thus avoiding particles ov-
    erlapping).
    """
    from random import seed, uniform
    import numpy as np

    ran = np.zeros(dim)              # array to store random numbers
    seed()

    for i in range(len(pos)):

        tooClose = True

        # check to make sure position is not too close to other particles
        while tooClose:

            tooClose = False
            # get a random position
            for j in range(dim):
                ran[j] = uniform(0.0,box[j])
            # assign position
            pos[i] = ran
            # check to see if random position is too close to other positions
            for k in range(i):
                ri = pos[i]
                rj = pos[k]
                rmag,runit = calc_rij_pbc(ri,rj,dim,box)
                if rmag < 2.0*rad:
                    tooClose = True
                    break

    return pos;

def pos_init_lat(pos,dim,box,rad,ao,N):
    """
    This function initializes the particle positions on a regular lattice with a
    small amount of random noise a fraction of the lattice constant 'ao'. Since
    the number of particles is supplied as an input, this function also calcu-
    lates and returns the box size. This requires N to be a clean root of order
    dim (ie N = 4,9,16,25,... for dim = 2 or N = 8,27,64,125,... for dim = 3)
    """
    import numpy as np
    from random import random

    # check to make sure input arguments are valid
    eflag = False
    if ao < 2.0*rad:
        eflag = True
        print 'ao is too small in relation to the particle radius'

    clean_root_check = N % int(round(N**(1.0/dim)))
    if clean_root_check != 0:
        eflag = True
        print 'N is not a clean root'

    # determine the box size based on ao,N, and dim.
    if dim == 1:
        box[0] = float(N*ao)
        print str(dim)+'-dimensional simulation with box length of '+str(L)+'.'
    elif dim == 2:
        Nd = int(np.sqrt(N))
        L = float(Nd*ao)
        box = np.array([L,L])
        print str(dim)+'-dimensional simulation with box length of '+str(L)+'.'
    else:
        Nd = int(round(N**(1.0/dim)))
        L = float(Nd*ao)
        box = np.array([L,L,L])
        print str(dim)+'-dimensional simulation with box length of '+str(L)+'.'

    noise_strength = 0.15*ao

    # assign positions with noise
    temp = np.zeros(dim)

    if dim == 1:
        for i in range(N):
            # generate random noise
            ran = random
            if ran < 0.5:
                ran = -1.0
            else:
                ran = 1.0
            # assign position
            pos[i] = ao/2.0 + i*ao + ran*noise_strength
    elif dim == 2:
        for j in range(Nd):
            for i in range(Nd):
                noise_dir = random_unit_vec(dim)
                ind = j*Nd + i
                xpos = ao/2.0 + i*ao
                ypos = ao/2.0 + j*ao
                pos[ind] = np.array([xpos,ypos]) + noise_strength*noise_dir
    else:
        for k in range(Nd):
            for j in range(Nd):
                for i in range(Nd):
                    noise_dir = random_unit_vec(dim)
                    ind = j*Nd + i
                    xpos = ao/2.0 + i*ao
                    ypos = ao/2.0 + j*ao
                    zpos = ao/2.0 + k*ao
                    pos[ind] = np.array([xpos,ypos,zpos]) + noise_strength*noise_dir

    return pos,box,eflag;

def calc_rij_pbc(ri,rj,dim,box):
    """
    This function takes the vector positions of two particles, i and j, and cal-
    culates their separation distance assuming periodic boundary conditions
    (pbc). It also calculates a unit vector, 'rij_unit' that points in the dir-
    ection of the force on particle i (with position ri).

    INPUTS:
    ri = position of particle i
    rj = position of particle j
    dim = the dimension of the simulation (i.e. 1-,2-, or 3-dimensions)
    box = array with 'dim' components. each component is the length of the box
            in the respective dimension.

    OUTPUTS:
    rij_mag = the magnitude of the separation vector with pbc.
    rij_unit = separtion vector with pbc made into a unit vector.

    """
    import numpy as n

    # initialize outputs
    rij_unit = n.zeros(dim)
    rij_mag = 0.0

    # loop over each component of the position vectors
    for d in range(dim):
        # get r1d - r2d, where d is the dimension, this points in the direction
        # of the repulsive force
        rij_unit[d] = ri[d] - rj[d]
        # get the magnitude, the sign indicates the direction
        rd_mag = abs(rij_unit[d])
        # if the magnitude is bigger than box[d]/2 then there is a smaller sep-
        # aration distance due to periodic boundary conditions.
        if rd_mag > box[d]/2.0:
            rd_mag = box[d] - rd_mag # use periodic image of rj[d]
            rij_unit[d] = -rij_unit[d] # switch the direction of the force

        rij_mag = rij_mag + rd_mag**2 # add up the square of the d components

    # take the sqrt of the sum of the squares of the components to get the mag
    rij_mag = n.sqrt(rij_mag)

    # calculate unit vector with pbc
    unit_mag = n.linalg.norm(rij_unit)
    rij_unit = 1/unit_mag*rij_unit

    return rij_mag,rij_unit;

def vel_init(vel,dim,sigma):
    """
    This function populates the velocity array with random velocities from
    a Gaussian distribution.
    """
    from random import seed, gauss
    import numpy as np

    ran = np.zeros(dim)              # array to store random numbers
    seed()

    for i in range(len(vel)):

        for j in range(dim):
            ran[j] = gauss(0.0,sigma)

        vel[i] = ran

    return vel;

def calc_repulsion_force(rij,rad,sigma,eps):
    """
    This function takes the separation distance between two particles, rij, the
    particle radius, rad, and the repulsive force constant, sigma, and calcu-
    lates the repulsion force between the two. It returns this force value.

    if we are calculating the force on particle i then rij = |ri - rj| and the
    force will point in the same direction as the vector ri - rj.
    """

    rep_force = 48.0*eps*(sigma/(rij-rad))**(13)

    return rep_force;

def write_point_data_vtk(description,numpy_array,N,dim,step):
    """
    This function takes a numpy array that stores 1-, 2-, or 3-dimensional point
    data as floats, and writes it to a vtkfile with the 'step' as part of the
    filename. The dimension is indicated by the input variable 'dim'.

    INPUTS:

    description = a string that will be the first part of the vtk filename.
    numpy_array = numpy array with shape = (N,dim).
    N = a scalar equal to the number of points in the numpy_array
    dim = scalar that is equal to 1, 2, or 3.
    step = a scalar that represents a step from a loop iteration variable.

    OUTPUTS:

    filename_step.vtk = a vtk file with the corresponding point data

    """
    from numpy import c_,zeros,array
    import os

    # if dim < 3, pad numpy_array with zeros making it 3-dimensional
    if dim == 1:
        numpy_array = c_[numpy_array,zeros((N,2))] # add 2 col of zeros
    if dim == 2:
        numpy_array = c_[numpy_array,zeros(N)] # add 1 col of zeros

    # construct the full filename
    file_name = description + '_' + str(step) + '.vtk'

    # open the file for writing
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dest_dir = os.path.join(script_dir, 'vtkoutput')
    try:
        os.makedirs(dest_dir)
    except OSError:
        pass # already exists
    path = os.path.join(dest_dir, file_name)
    with open(path, 'wb') as f:

        # write vtk header
        f.write('# vtk DataFile Version 3.1\n')
        f.write('VTK file containing ' + description + ' data\n' )
        f.write('ASCII\n')
        f.write('DATASET POLYDATA\n')
        f.write('\n')
        f.write('POINTS\t' + str(N) + '\tfloat\n')

        # write the point data
        for i in range(N):
            f.write('{:>20}'.format(str(numpy_array[i,0]))),
            f.write('{:>20}'.format(str(numpy_array[i,1]))),
            f.write('{:>20}'.format(str(numpy_array[i,2])))
            f.write('\n')

    # close file
    f.close()

    return 0;

def enforce_pbc(pos,dim,box):
    """
    This function takes a particles position and enforces periodic boundary
    conditions according to the 'box' parameter.

    INPUTS:

    pos = numpy array of shape (dim)
    dim = scalar value equal to 1, 2, or 3 representing the dimensionality of
            the simulation (i.e. if dim = 2, pos has 2 components)
    box = numpy array with of shape (dim) that contains the simulation domain
            size in the respective dimensions (i.e. if dim = 2, then
            box = [Lx,Ly])

    OUTPUTS:

    pos = numpy array of shape (dim) altered to comply with periodic boundary
            conditions.
    """
    for d in range(dim):
        if pos[d] > box[d]:
            pos[d] = pos[d] - box[d]
        if pos[d] < 0.0:
            pos[d] = pos[d] + box[d]

    return pos;

def random_unit_vec(dim):
    """
    This function generates a unit vector of dimension = dim that points in a
    random direction.
    """
    import numpy as np
    import random as r

    r.seed()
    vec = np.zeros(dim)
    for d in range(dim):
        vec[d] = r.uniform(-1.0,1.0)

    vec = 1.0/np.linalg.norm(vec)*vec

    return vec;
