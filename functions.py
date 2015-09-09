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

def pos_init(pos,dim,L):
    """
    This function populates the position array with random positions from
    uniformly distributed throughout the box.
    """
    from random import seed, uniform
    import numpy as np

    ran = np.zeros(dim)              # array to store random numbers
    seed()

    for i in range(len(pos)):

        for j in range(dim):
            ran[j] = uniform(0.0,L)

        pos[i] = ran

    return pos;

def calc_repulsion_force(rij,rad,sigma):
    """
    This function takes the separation distance between two particles, rij, the
    particle radius, rad, and the repulsive force constant, sigma, and calcu-
    lates the repulsion force between the two. It returns this force value.

    if we are calculating the force on particle i then rij = |ri - rj| and the
    force will point in the same direction as the vector ri - rj.
    """

    rep_force = sigma/(rij-2*rad)**7

    return rep_force;

def calc_rij_pbc(ri,rj,dim,box):
    """
    This function takes the vector positions of two particles, i and j, and cal-
    culates their separation distance assuming periodic boundary conditions
    (pbc). It also calculates the unit vector of the of the rij_pbc vector,
    'rij_unit'.

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

    # implement periodic boundary conditions
    rij_pbc = n.zeros(dim)
    rij_mag = 0.0
    for d in range(dim):
        rij_pbc[d] = ri[d] - rj[d]
        if abs(rij_pbc[d]) > box[d]/2:
            rij_pbc[d] = rij_pbc[d] + box[d]
        rij_mag = rij_mag + rij_pbc[d]**2

    rij_mag = n.sqrt(rij_mag)

    # calculate unit vector with pbc
    rij_unit = 1/rij_mag*rij_pbc

    return rij_mag,rij_unit;
