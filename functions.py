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

def calc_repulsion_force(rij,FR):
    """
    This function takes the separation distance between two particles, rij and
    the repulsive force constant, FR, and calculates the repulsion force between
    the two. It returns this force value.

    if we are calculating the force on particle i then rij = |ri - rj| and the
    force will point in the same direction as the vector ri - rj.
    """

    rep_force = FR/rij**7

    return rep_force;
