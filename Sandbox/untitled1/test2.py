import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

big_number = 10000.0
width = 0.0001

ts = np.linspace(0, 10, 2000)
def f(X, t):
    dx0 = X[1]
    dx1 = -9.83
    dx1 += big_number / (1 + np.exp(X[0]/width))
    return [dx0, dx1]

with np.errstate(over='ignore'):
    # don't print overflow warning messages from exp(),
    # and limit the step size so that the solver doesn't step too deep
    # into the forbidden region
    X = odeint(f, [2, 0], ts, hmax=0.01)
plt.plot(ts, X[:, 0])
plt.show()