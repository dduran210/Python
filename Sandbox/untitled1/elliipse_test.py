import numpy as np
import matplotlib.pyplot as plt

a = 10.0
b = 5.0


pi = np.pi
nVec = 1000
t = np.arange(0, 2 * pi, (2 * pi) / nVec)

x0 = 0
y0 = 0

x = x0 + a * np.cos(t)
y = y0 + b * np.sin(t)

plt.figure()
plt.plot(x, y)

plt.show()