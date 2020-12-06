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
Xbig = jaribig * R * np.sin(T)
Ybig = jaribig * R * np.cos(T)
Zbig = tbig * R

X1 = jaritop * R * np.sin(T)
Y1 = jaritop * R * np.cos(T)
Z1 = t * R + (tbig - t)

X2 = jaribase * R * np.sin(T)
Y2 = jaribase * R * np.cos(T)
Z2 = t * R

X1 = np.flip(X1, 0)
Y1 = np.flip(Y1, 0)
X2 = np.flip(X2, 0)
Y2 = np.flip(Y2, 0)
Xbig = np.flip(Xbig, 0)
Ybig = np.flip(Ybig, 0)


#Fault
dip = np.arccos(jaritop/((jaritop**2 + t**2)**0.5))*180/np.pi

Xf = np.linspace(-jaritop, jaritop, length)
Yf = np.linspace(-jaritop, jaritop, length)
Xf, Yf = np.meshgrid(Xf, Yf)

Zf = Xf*np.tan(dip*np.pi/180) + (tbig - t)

# Slice Big
now = -6.5
for i in range(len(Xbig)):
    for j in range(len(Xbig[i])):
        if(0 <= round(Xbig[i][j]) < len(Zf[0])):
            if(Xbig[i][j] > now and round(Zbig[i][0],1) >= round(Zf[0][round(Xbig[i][j])],1)):
                Xbig[i][j] = None
                Ybig[i][j] = None
            elif(Xbig[i][j] > now):
                Xbig[i][j] = None
                Ybig[i][j] = None
        elif(round(Xbig[i][j]) >= len(Zf[0]) or round(Xbig[i][j]) < 0):
            if(Xbig[i][j] > now):
                Xbig[i][j] = None
                Ybig[i][j] = None

    now += jaribig/length


# Slice Top
now = 0
for i in range(len(X1)):
    for j in range(len(X1[i])):
    # if(Zf[0][0] <= round(X1[i][j]) <= Zf[0][-1]):
        if(X1[i][j] > now and round(Z1[i][0],1) > round(Zf[0][round(X1[i][j])],1)):
            X1[i][j] = None
            Y1[i][j] = None
    now += jaritop/length

# Slice Base
now = 0
for i in range(len(X2)):
    for j in range(len(X2[i])):
    # if(Zf[0][0] <= round(X2[i][j]) <= Zf[0][-1]):
        if(X2[i][j] > now and round(Z2[i][0],1) > round(Zf[0][round(X2[i][j])],1)):
            X2[i][j] = None
            Y2[i][j] = None
    now += jaribase/length

ax.plot_wireframe(X1, Y1, Z1, color="blue")
ax.plot_wireframe(Xbig, Ybig, Zbig, color="brown", alpha=0.3)

# ax.plot_wireframe(X2, Y2, Z2, color="red", alpha=0.5)
ax.plot_wireframe(Xf, Yf, Zf, alpha=0.5, color='black', label='fault')

ax.set_xlabel("EASTING(m)")
ax.set_ylabel("NORTHING(m)")
ax.set_zlabel("ELEVATION(m)")
plt.show()