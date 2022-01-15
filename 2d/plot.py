import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

file1 = 'results/eq_results.txt'
file2 = 'results/results.txt'
file3 = 'results/forces.txt'
file4 = 'results/max_radius.txt'

data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)
data3 = np.loadtxt(file3)
max_r = np.loadtxt(file4)

y1 = np.asarray(data1)
x1 = np.linspace(0, max_r, 10000)

print (max_r)

x2 = np.asarray(data2[0])
y2 = np.asarray(data2[1])

x3 = np.asarray(data3[0])
y3 = np.asarray(data3[1])

plt.show()

plt.subplot(1, 3, 1)
plt.plot(x1, y1)

plt.subplot(1, 3, 2)
plt.plot(x2, y2)

plt.subplot(1, 3, 3)
plt.plot(x3, y3)

plt.show()
