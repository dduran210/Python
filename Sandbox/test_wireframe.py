import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

theta, phi = np.linspace(0, 2 * np.pi, 13), np.linspace(0, np.pi, 7)
THETA, PHI = np.meshgrid(theta, phi)
R = 1.0
X = R * np.sin(PHI) * np.cos(THETA)
Y = R * np.sin(PHI) * np.sin(THETA)
Z = R * np.cos(PHI)
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
plot = ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1, linewidth=0, color='white', antialiased=False, alpha=1, frameon=False)

plt.show()