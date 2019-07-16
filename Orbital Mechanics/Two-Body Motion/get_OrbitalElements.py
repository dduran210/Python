import numpy as np


def get_SpecificAngularMomentum_vec(r, v):
    return np.cross(r, v)



def get_2BodyParameters_scalar():
    p =  0
    e = 0
    a = 0
    b = 0

    r_p = p / (1 + e)
    ra = p / (1 - e)




