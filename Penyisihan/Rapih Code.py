import numpy as np
import matplotlib.pyplot as plt


#Alpha Parameters
rtop_alpha = 4500 #in m
rbase_alpha = 3900 #in m
htop_alpha = 800 #in m
hbase_alpha = 700 #in m
ntg_alpha = 0.35 #in decimal
por_alpha = 0.25 #in decimal
so_alpha = 0.8  #in decimal
re_alpha = 0.25 #in decimal
fvf_alpha = 0.9 #in decimal
oiltrapped_alpha = 0.62 #in decimal

#Beta Parameters
rtop_beta = 4000 #in m
rbase_beta = 3500 #in m
htop_beta = 800 #in m
hbase_beta = 700 #in m
ntg_beta = 0.3 #in decimal
por_beta = 0.22 #in decimal
so_beta = 0.8 #in decimal
re_beta = 0.2 #in decimal
fvf_beta = 0.9 #in decimal
oiltrapped_beta = 0.62 #in decimal

#GRV Calculation
def grv_calc (rtop, rbase, htop, hbase) :
    vol = 1/3*np.pi*(rtop**2*htop - rbase**2*hbase)
    return vol

# Slicing
def slice_arr(start, end, arr):
    for i in range(len(arr)):
        arr[i][start:end] = None
    return arr

def create_cone(R, T, r, H, hmax, hmin):
    X = r * R * np.sin(T)
    Y = r * R * np.cos(T)
    Z = H * R - (hmax - hmin)
    return X, Y, Z

def main():
    grv_alpha = grv_calc(rtop_alpha, rbase_alpha, htop_alpha, hbase_alpha) #meter cubic
    grv_beta =  grv_calc(rtop_beta, rbase_beta, htop_beta, hbase_beta) #meter cubic
    
    #VISUALIZATION
    theta = np.linspace(0,2*np.pi,90)
    r = np.linspace(0,1,10)
    T, R = np.meshgrid(theta, r)

    #Big Cone
    h_big_alpha = htop_alpha + 70/100*htop_alpha
    r_big_alpha = rtop_alpha/htop_alpha * h_big_alpha

    X_big_alpha, Y_big_alpha, Z_big_alpha = create_cone(R, T, r_big_alpha, h_big_alpha, h_big_alpha, htop_alpha)


    # Slice Open surface
    X_big_alpha = slice_arr(20, 60, X_big_alpha)
    Y_big_alpha = slice_arr(20, 60, Y_big_alpha)

    #Little Cone
    h_lit_alpha = hbase_alpha + 70/100*hbase_alpha
    r_lit_alpha = rbase_alpha/hbase_alpha * h_lit_alpha

    X_lit_alpha, Y_lit_alpha, Z_lit_alpha = create_cone(R, T, r_lit_alpha, h_lit_alpha, h_lit_alpha, hbase_alpha)


    #Alpha Cone(Oil)
    X_top_alpha, Y_top_alpha, Z_top_alpha = create_cone(R, T, rtop_alpha, htop_alpha, 0, 0)

    X_big_alpha = np.flip(X_big_alpha, 0)
    Y_big_alpha = np.flip(Y_big_alpha, 0)
    X_top_alpha = np.flip(X_top_alpha, 0)
    Y_top_alpha = np.flip(Y_top_alpha, 0)
    X_lit_alpha = np.flip(X_lit_alpha, 0)
    Y_lit_alpha = np.flip(Y_lit_alpha, 0)

    #Fault
    dip = np.arccos(rtop_alpha/((rtop_alpha**2 + htop_alpha**2)**0.5))*180/np.pi

    Xf = np.linspace(-r_big_alpha, r_big_alpha, 10)
    Yf = np.linspace(-r_big_alpha, r_big_alpha, 10)
    Xf, Yf = np.meshgrid(Xf, Yf)
    Zf = (Xf-rtop_alpha)*np.tan(dip*np.pi/180)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')


    ax.plot_surface(X_lit_alpha, Y_lit_alpha, Z_lit_alpha,
                    cmap='gist_heat')

    ax.plot_surface(Xf, Yf, Zf, antialiased=True, color="grey")

    ax.plot_surface(X_top_alpha, Y_top_alpha, Z_top_alpha,
                    cmap='summer')

    ax.plot_surface(X_big_alpha, Y_big_alpha, Z_big_alpha,
                    cmap='gist_heat')

    ax.view_init(47,-126)

    ax.set_xlabel("EASTING(m)")
    ax.set_ylabel("NORTHING(m)")
    ax.set_zlabel("ELEVATION(m)")
    plt.show()

if __name__ == "__main__":
    main()