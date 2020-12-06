import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

#Kerucut
t = 0.8
jaritop = 4.5
jaribase = 3.9

tbig = 2
jaribig = jaritop/t * tbig

length = 1000

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#ALPHA
theta = np.linspace(0,2*np.pi,90)
r = np.linspace(0,1,length)
T, R = np.meshgrid(theta, r)

# Then calculate X, Y, and Z
X0 = jaribig * R * np.sin(T)
Y0 = jaribig * R * np.cos(T)
Z0 = tbig * R

X1 = jaritop * R * np.sin(T)
Y1 = jaritop * R * np.cos(T)
Z1 = t * R + (tbig - t)

X2 = jaribig * R * np.sin(T)
Y2 = jaribig * R * np.cos(T)
Z2 = tbig * R

X0 = np.flip(X0, 0)
Y0 = np.flip(Y0, 0)
X1 = np.flip(X1, 0)
Y1 = np.flip(Y1, 0)
X2 = np.flip(X2, 0)
Y2 = np.flip(Y2, 0)


#Fault
dip = np.arccos(jaritop/((jaritop**2 + t**2)**0.5))*180/np.pi

Xf = np.linspace(-jaritop, jaritop, length)
Yf = np.linspace(-jaritop, jaritop, length)
Xf, Yf = np.meshgrid(Xf, Yf)
Zf = Xf*np.tan(dip*np.pi/180) + (tbig - t)

# SLice Surface
for i in range(len(X0)):
    X0[i][30:56] = None

ax.plot_surface(X2, Y2, Z2, color="purple", alpha=1)
ax.plot_surface(X1, Y1, Z1, color="yellow", alpha=1)
ax.plot_surface(X0, Y0, Z0, color="brown", alpha=1)
ax.plot_wireframe(Xf, Yf, Zf, alpha=0.5, color='black', label='fault')

ax.set_xlabel("EASTING(m)")
ax.set_ylabel("NORTHING(m)")
ax.set_zlabel("ELEVATION(m)")
plt.show()