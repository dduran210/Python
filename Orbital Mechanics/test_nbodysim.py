import numpy as np
import matplotlib.pyplot as plt

class particle:
    def __init__(self, x=np.zeros((6, 1)), m=5.9722e24):
        self.x = x
        self.m = m
        self.A = np.zeros(6)

    def r_1to2(self, p2):
        return p2.x[0:3] - self.x[0:3]

    def set_x(self, x):
        self.x = x

    def generate_A(self, lst_particles):
        G = (6.67408e-11) / (1000 ** 3)

        # self.A = np.hstack((np.zeros(3), np.identity(3)))

        ddx = np.zeros((3, 1))
        for k_, p in enumerate(lst_particles):
            r = self.r_1to2(p)
            r_mag = np.linalg.norm(r)
            ddx += ((G * p.m) / (r_mag) ** 3) * r

        self.A = lambda x, t: np.array([[0, 0, 0, 1, 0, 0],
                                        [0, 0, 0, 0, 1, 0],
                                        [0, 0, 0, 0, 0, 1],
                                        [ddx[0], 0, 0, 0, 0, 0],
                                        [0, ddx[1], 0, 0, 0, 0],
                                        [0, 0, ddx[2], 0, 0, 0]]).dot(np.vstack((np.ones((3, 1)), x[3:].reshape(-1, 1))))




def RV2COE(r_IJK, v_IJK):
    mu_ = 398600.4418 # km^3 / sec^2

    K = np.array([0, 0, 1])

    h = np.cross(r_IJK, v_IJK)
    n = np.cross(K, h)

    r_mag = np.linalg.norm(r_IJK)
    v_mag = np.linalg.norm(v_IJK)
    h_mag = np.linalg.norm(h)
    n_mag = np.linalg.norm(n)

    e = ((v_mag**2 - mu_ / r_mag) * r_IJK - np.dot(r_IJK, v_IJK) * v_IJK) / mu_
    e_mag = np.linalg.norm(e)


    if e_mag != 1.0:
        epsilon = v_mag ** 2 / 2 - mu_ / r_mag
        a = -mu_ / (2 * epsilon)
    else:
        a = float("inf")

    p = h_mag ** 2 / mu_


    i = np.rad2deg(np.arccos(h[2] / h_mag))
    omega = np.rad2deg(np.arccos(n[0] / n_mag))
    if n[1] < 1:
        omega = 360 - omega

    w = np.rad2deg(np.arccos(np.dot(n, e) / (n_mag * e_mag)))
    if e[2] < 0:
        w = 360 - w

    v = np.rad2deg(np.arccos(np.dot(e, r_IJK) / (e_mag * r_mag)))



    return p, a, e, i, omega, w, v

def COE2RV(p, e, i, RAAN, omega, nu):
    mu_ = 398600.4418  # km^3 / sec^2
    # i_rad = np.deg2rad(i)
    # RAAN_rad = np.deg2rad(RAAN)
    # omega_rad = np.deg2rad(omega)
    nu_rad = np.deg2rad(nu)

    r_PQW = np.array([[p * np.cos(nu_rad) / (1 + e * np.cos(nu_rad))],
                      [p * np.sin(nu_rad) / (1 + e * np.cos(nu_rad))],
                      [0]])
    v_PQW = np.array([[-np.sqrt(mu_ / p) * np.sin(nu_rad)],
                      [np.sqrt(mu_ / p) * (e + np.cos(nu_rad))],
                      [0]])


    T = np.matmul(ROT(3, -RAAN), np.matmul(ROT(1, -i), ROT(3, -omega)))

    r_IJK = T.dot(r_PQW)
    v_IJK = T.dot(v_PQW)
    return r_IJK, v_IJK

def ROT(i_axis=1, theta_deg=0.0):
    theta_rad = np.deg2rad(theta_deg)
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)

    if i_axis == 1:
        mat_rot = np.array([[1, 0, 0],
                            [0, c, s],
                            [0, -s, c]])

    elif i_axis == 2:
        mat_rot = np.array([c, 0, s],
                           [0, 1, 0],
                           [-s, 0, c])

    elif i_axis == 3:
        mat_rot = np.array([[c, s, 0],
                            [-s, c, 0],
                            [0, 0, 1]])

    return mat_rot

def RK4_step(G, y, t, dt):
	k1 = G(y, t)
	k2 = G(y + 0.5 * k1 * dt, t + 0.5 * dt)
	k3 = G(y + 0.5 * k2 * dt, t + 0.5 * dt)
	k4 = G(y + k3 * dt, t + dt)

	#return dt * G(y,t)
	return dt * (k1 + 2*k2 + 2*k3 + k4) /6

def RV2COE_test():
    r_IJK = np.array([6524.834, 6862.875, 6448.296])
    v_IJK = np.array([4.901327, 5.533756, -1.976341])
    RV2COE(r_IJK, v_IJK)

def COE2RV_test():
    p = 11067.790
    e = 0.83285
    i = 87.87
    RAAN = 227.89
    omega = 53.38
    nu = 92.335
    r_IJK, v_IJK = COE2RV(p, e, i, RAAN, omega, nu)

def f_xdotdot(lst_particles, i_):
    G = (6.67408e-11) / (1000 ** 3)  # km^3 / (kg * s)
    xdotdot = np.zeros((3, 1))
    p_i = lst_particles[i_]
    for k_, p in enumerate(lst_particles):
        if i_ != k_:
            r = p_i.r_1to2(p)
            r_mag = np.linalg.norm(r)
            xdotdot += ((G * p.m) / (r_mag) ** 3) * r


    # check_x1 = xdotdot[0, 0]
    # check_x2 = xdotdot[1, 0]
    # check_x3 = xdotdot[2, 0]
    # check_x = p_i.x[0:3]
    A = lambda x, t: np.array([[0, 0, 0, 1, 0, 0],
                               [0, 0, 0, 0, 1, 0],
                               [0, 0, 0, 0, 0, 1],
                               [xdotdot[0, 0], 0, 0, 0, 0, 0],
                               [0, xdotdot[1, 0], 0, 0, 0, 0],
                               [0, 0, xdotdot[2, 0], 0, 0, 0]]).dot(np.vstack((np.ones((3, 1)), x[3:].reshape(-1, 1))))

    # check = np.array([[0, 0, 0, 1, 0, 0],
    #                            [0, 0, 0, 0, 1, 0],
    #                            [0, 0, 0, 0, 0, 1],
    #                            [xdotdot[0, 0], 0, 0, 0, 0, 0],
    #                            [0, xdotdot[1, 0], 0, 0, 0, 0],
    #                            [0, 0, xdotdot[2, 0], 0, 0, 0]]).dot(np.vstack((p_i.x[0:3], np.ones((3, 1)))))

    return A

def nbodysim_test():
    G = (6.67408e-11) / (1000 ** 3)# km^3 / (kg * s)

    p = 11067.790
    e = 0.83285
    i = 87.87
    RAAN = 227.89
    omega = 53.38
    nu = 92.335
    r_IJK, v_IJK = COE2RV(p, e, i, RAAN, omega, nu)

    T = np.matmul(ROT(3, -RAAN), np.matmul(ROT(1, -i), ROT(3, -omega)))
    T_inv = np.linalg.inv(T)
    r_PQW = T_inv.dot(r_IJK)
    v_PQW = T_inv.dot(v_IJK)
    x0 = np.vstack((r_PQW, v_PQW))

    p_earth = particle()
    p_i = particle(x=np.vstack((r_PQW, v_PQW)), m=1000)



    # xdotdot = lambda x: ((-G * p_earth.m) / (np.linalg.norm(x[0:3]) ** 3))
    # A = lambda x, t: np.array([[0, 0, 0, 1, 0, 0],
    #                            [0, 0, 0, 0, 1, 0],
    #                            [0, 0, 0, 0, 0, 1],
    #                            [xdotdot(x), 0, 0, 0, 0, 0],
    #                            [0, xdotdot(x), 0, 0, 0, 0],
    #                            [0, 0, xdotdot(x), 0, 0, 0]]).dot(x)






    t_start = 0
    t_stop = 5000
    dt = 0.1
    t_vec = np.arange(t_start, t_stop, dt)
    y = np.zeros((6, t_vec.size))
    y[:, 0] = x0[:, 0]

    lst_particles = [p_i, p_earth]


    for k in np.arange(1, t_vec.size):
        t = t_vec[k - 1]
        #
        # check_x = A(y[:, k - 1].reshape(-1, 1), t)

        # check_x1_step = RK4_step(A, y[:, k - 1], t, dt)
        # check_x2_step = RK4_step(A_fx, y[:, k - 1].reshape(-1, 1), t, dt)

        # y[:, k] = y[:, k - 1] + RK4_step(A, y[:, k - 1], t, dt)

        A = f_xdotdot(lst_particles, 0)
        temp = y[:, k - 1].reshape(-1, 1) + RK4_step(A, y[:, k - 1].reshape(-1, 1), t, dt)
        y[:, k] = temp[:, 0]
        p_i.x = y[:, k].reshape(-1, 1)




    fig = plt.figure()



    x_2body, y_2body = PQW_xy(p, e)

    plt.plot(x_2body, y_2body)
    plt.plot(y[0, :], y[1, :])
    plt.show()


def plot(y, t_vec):
    fig = plt.figure()
    plt.plot(y[0, :], y[1, :])
    # plt.show()

def PQW_xy(p, e):
    N = 1000
    x = np.zeros(N)
    y = np.zeros(N)
    k = 0
    for nu_rad in np.linspace(0, 2 * np.pi, N):
        r = p / (1 + e * np.cos(nu_rad))
        x[k] = r * np.cos(nu_rad)
        y[k] = r * np.sin(nu_rad)
        k += 1

    return x, y
    # fig = plt.figure()
    # plt.plot(x, y)
    # plt.show()



if __name__ == '__main__':

    # RV2COE_test()
    # COE2RV_test()
    nbodysim_test()