#from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

file = 'results/tree_test_run.txt'

data = np.loadtxt(file)

f = np.asarray(data[0])
x = np.asarray(data[1])
y = np.asarray(data[2])
z = np.asarray(data[3])

colormap = plt.cm.plasma
normalize = colors.PowerNorm(gamma=0.3, vmin=f.min(), vmax=f.max())

fig1 = plt.figure(1)
ax = fig1.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c=f, cmap=colormap, norm=normalize, s=5)
ax.set_xlim3d(-100, 100)
ax.set_ylim3d(-100, 100)
ax.set_zlim3d(-100, 100)
ax.set_axis_off()

fig2 = plt.figure(2)
bx = fig2.add_subplot(111, projection='3d')
bx.scatter(x, y, z, c=f, cmap=colormap, norm=normalize, s=5)
bx.set_xlim3d(-0.8, 0.8)
bx.set_ylim3d(-0.8, 0.8)
bx.set_zlim3d(-0.8, 0.8)
bx.set_axis_off()

fig3 = plt.figure(3)
cx = fig3.add_subplot(111, projection='3d')
cx.scatter(x, y, z, c=f, cmap=colormap, norm=normalize, s=5)
cx.set_xlim3d(-0.02, 0.02)
cx.set_ylim3d(-0.02, 0.02)
cx.set_zlim3d(-0.02, 0.02)
cx.set_axis_off()

fig4 = plt.figure(4)
dx = fig4.add_subplot(111)
dx.scatter(x, y, c=f, cmap=colormap, norm=normalize, s=5)
plt.xlim(-2, 2)
plt.ylim(-2, 2)
dx.set_axis_off()

plt.show()
