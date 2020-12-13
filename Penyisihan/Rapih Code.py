import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as line
import pandas as pd
from pert import PERT
from scipy.special import gamma
import statistics as stats
from scipy.integrate import simps

reservoir_thickness = 100 #in m

# Alpha Parameters
wc_alpha = 2600 # in m
rtop_alpha = 4500  # in m
rbase_alpha = 3900  # in m
htop_alpha = 800  # in m
hbase_alpha = htop_alpha - reservoir_thickness  # in m
ntg_alpha = 0.35  # in decimal
por_alpha = 0.25  # in decimal
so_alpha = 0.8  # in decimal
re_alpha = 0.25  # in decimal
fvf_alpha = 0.9  # in decimal
chance_alpha = 0.62  # in decimal

# Beta Parameters
wc_beta = 1800 #in m
rtop_beta = 4000  # in m
rbase_beta = 3500  # in m
htop_beta = 800  # in m
hbase_beta = htop_beta-reservoir_thickness  # in m
ntg_beta = 0.3  # in decimal
por_beta = 0.22  # in decimal
so_beta = 0.8  # in decimal
re_beta = 0.2  # in decimal
fvf_beta = 0.9  # in decimal
chance_beta = 0.62  # in decimal

####ECONOMICS ANALYSIS VARIABLE####
oilprice = 50 # in USD per barrel
project_time = 20 #in years
capex = 900000000 #in USD, Pso(MMUSD 800) and drilling 10 well (@ MMUSD100/well)
opex = 300000000 #in USD, facilities, electricity, etc
tax = 20/100

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

#########################################################################################################################################################################

def no_zero(array):
    rbaru = []
    for i in array:
        if i < 0:
            rbaru.append(0)
        else:
            rbaru.append(i)
    return np.array(rbaru)

def simulation_pert (low,most,high):
    pert = PERT(low, most, high)
    r = pert.rvs(200000)
    count, bins, ignored = plt.hist(r, 30, density=True, color='indigo', rwidth=0.9)
    alpha = 1 + 4*((most - low)/(high-low))
    beta = 1 + 4*((high - most)/(high-low))
    B = (gamma(alpha)*gamma(beta))/(gamma(alpha + beta))
    pdf = (((bins-low)**(alpha-1))*((high-bins)**(beta-1)))/(B*((high-low)**(alpha+beta-1)))
    mode = most

    plt.plot(bins, pdf, linewidth=1, color='r')
    plt.plot(np.linspace(np.percentile(r, 10), np.percentile(r, 10), len(count)), count,
             label=('P10 =', round(np.percentile(r, 10), 2)))
    plt.plot(np.linspace(np.percentile(r, 50), np.percentile(r, 50), len(count)), count,
             label=('P50 =', round(np.percentile(r, 50), 2)))
    plt.plot(np.linspace(np.percentile(r, 90), np.percentile(r, 90), len(count)), count,
             label=('P90 =', round(np.percentile(r, 90), 2)))
    plt.plot(np.linspace(np.mean(r), np.mean(r), len(count)), count,
             label=('MEAN =', round(np.mean(r), 2)))
    plt.plot(bins,pdf, color='orchid')
    plt.legend()

    return    {
            "res": r,
            "val_min": round(min(r), 2),
            "val_max": round(max(r), 2),
            "mean": round(np.mean(r), 2),
            "mode": round(mode, 2),
            "p50": round(np.percentile(r, 50), 2),
            "p1": round(np.percentile(r, 1), 2),
            "p10": round(np.percentile(r, 10), 2),
            "p15": round(np.percentile(r, 15), 2),
            "p85": round(np.percentile(r, 85), 2),
            "p90": round(np.percentile(r, 90), 2),
            "p99": round(np.percentile(r, 99), 2),
            }

# GRV Calculation
def grv_calc(rtop, rbase, htop, hbase):
    vol = round((1 / 3 * np.pi * (rtop ** 2 * htop - rbase ** 2 * hbase))/1233.48)
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

def volumetric_deterministic(GRV, NtoG, por, so, re, fvf, chance) :
    vol = round(((7758*(GRV)*NtoG*por*so*re*chance)/((10**6)*fvf)),2) #GRV in acre feet
    return vol

def volumetric_probabilistic(area_top, area_base, htop, rc, ntg, por, so, re, chance, fvf) :
    recoverable_vol_risked = ((7758*((1 / 3 * (area_top["res"] * htop - area_base["res"] * (htop-rc["res"])))/1233.48)*ntg*por*so*re*chance)/((10**6)*fvf))
    count, bins, ignored = plt.hist(recoverable_vol_risked, 50, density=True, color='indigo', rwidth=0.9)
    plt.plot(np.linspace(np.percentile(recoverable_vol_risked,10),np.percentile(recoverable_vol_risked,10), len(count)), count,
            label=("P10", round(np.percentile(recoverable_vol_risked, 10),2)))
    plt.plot(np.linspace(np.percentile(recoverable_vol_risked,50),np.percentile(recoverable_vol_risked,50), len(count)), count,
             label=("P50", round(np.percentile(recoverable_vol_risked, 50),2)))
    plt.plot(np.linspace(np.percentile(recoverable_vol_risked,90),np.percentile(recoverable_vol_risked,90), len(count)), count,
             label=("P90", round(np.percentile(recoverable_vol_risked, 90),2)))
    plt.plot(np.linspace(np.mean(recoverable_vol_risked),np.mean(recoverable_vol_risked), len(count)), count,
             label=("MEAN", round(np.mean(recoverable_vol_risked),2)))
    plt.legend()
    plt.xlim(0,)
    plt.ylim(0,max(count+count*5/100))
    return {
        "res": recoverable_vol_risked,
        "val_min": round(min(recoverable_vol_risked), 2),
        "val_max": round(max(recoverable_vol_risked), 2),
        "mean": round(np.mean(recoverable_vol_risked), 2),
        "p50": round(np.percentile(recoverable_vol_risked, 50), 2),
        "p1": round(np.percentile(recoverable_vol_risked, 1), 2),
        "p10": round(np.percentile(recoverable_vol_risked, 10), 2),
        "p15": round(np.percentile(recoverable_vol_risked, 15), 2),
        "p85": round(np.percentile(recoverable_vol_risked, 85), 2),
        "p90": round(np.percentile(recoverable_vol_risked, 90), 2),
        "p99": round(np.percentile(recoverable_vol_risked, 99), 2),
        # "pdf": [bins, pdf],
    }

    #alpha_det
def break_even_analysis(vol, capex, opex, price, projecttime, tax):
    if isinstance(vol,dict):
        vol = vol["mean"]
    else:
        vol = vol
    
    total_income = vol*(10**6)*price #in USD
    income_per_year = total_income/projecttime
    income_per_year_aftertax = income_per_year - tax*income_per_year
    profit_per_year = income_per_year_aftertax - opex

    year = np.arange(0,projecttime,0.1)
    total_cost = capex + opex*year
    total_revenue = income_per_year_aftertax*year
    x_BEP = capex/(income_per_year_aftertax-opex)
    y_BEP = capex/(profit_per_year/income_per_year_aftertax)

    def tot_cost(yearr):
        return capex + opex*yearr
    def tot_revenue(yearr):
        return income_per_year_aftertax*yearr

    x = np.linspace(x_BEP, project_time, 10)
    tot_cos = simps(tot_cost(x), x)
    tot_rev = simps(tot_revenue(x), x)
    total_profit = tot_rev - tot_cos

    plt.plot(year, total_cost, label='Total Cost', linewidth=3)
    plt.plot(year, total_revenue, label='Total Revenue', linewidth=3)
    plt.scatter(x_BEP, y_BEP, label="Break Even Point", color='navy', linewidths=5)
    plt.plot((0, x_BEP), (y_BEP, y_BEP), color='black', linestyle='--')
    plt.plot((x_BEP, x_BEP), (0, y_BEP), color='black', linestyle='--')
    plt.text(x_BEP + 5 * x_BEP / 100, 0, round(x_BEP, 2))
    plt.fill_between(year, total_cost, total_revenue, where=(year < x_BEP), facecolor='pink', alpha=0.5, label='Loss')
    plt.fill_between(year, total_cost, total_revenue, where=(year > x_BEP), facecolor='green', alpha=0.5,
                     label='Profit')
    plt.xlabel('Year', fontsize=15)
    plt.ylabel('USD', fontsize=15)
    plt.legend(fontsize=20)
    plt.xlim(0, )
    plt.ylim(0, )
    return {
        "payout": x_BEP,
        "total profit": total_profit
    }
    
def sorting_sensitivity(valuesx, variablesx, lowx):
    datax = []
    for i in range(len(valuesx)):
        datax.append([valuesx[i], variablesx[i], lowx[i]])

    datax.sort(reverse=True)

    values_sorted = []
    variables_sorted = []
    low_sorted = []
    for i in datax:
        values_sorted.append(i[0])
        variables_sorted.append(i[1])
        low_sorted.append(i[2])
    return (values_sorted, variables_sorted, low_sorted)

def sensitivity_calculation(data, mean_resource):
    val_min = data["p10"]
    val_max = data["p90"]
    mean = data["mean"]
    median = data["p50"]

    sen_min = (mean_resource / mean * val_min) - mean_resource
    sen_max = (mean_resource / mean * val_max) - mean_resource
    value = np.absolute(sen_min) + np.absolute(sen_max)

    return {
        "sen_min": sen_min,
        "sen_max": sen_max,
        "value": value
    }

# Function to plot the bars, one by one
def plot_bars(ys, lows, values, base):
    for y, low, value in zip(ys, lows, values):
        # The width of the 'low' and 'high' pieces
        low_width = (base - low)
        high_width = (low + value - base)

        # Each bar is a "broken" horizontal bar chart
        plt.broken_barh(
            [(low, low_width), (base, high_width)],
            (y - 0.4, 0.8),
            facecolors=['blue', 'blue'],  # Try different colors if you like
            edgecolors=['black', 'black'],
            linewidth=1,
        )

        # Display the value as text. It should be positioned in the center of
        # the 'high' bar, except if there isn't any room there, then it should be
        # next to bar instead.
        x = base + high_width / 2
        if x <= base + 15:
            x = base + high_width + 7
        plt.text(x, y, str(value), va='center', ha='center')

#########################################################################################################################################################################
def monte_carlo_main():
    grv_alpha = grv_calc(rtop_alpha, rbase_alpha, htop_alpha, hbase_alpha)  # meter cubic
    grv_beta = grv_calc(rtop_beta, rbase_beta, htop_beta, hbase_beta)  # meter cubic

    volumetric_det_alpha = volumetric_deterministic(grv_alpha, ntg_alpha, por_alpha, so_alpha, re_alpha, fvf_alpha, chance_alpha)
    volumetric_det_beta = volumetric_deterministic(grv_beta, ntg_beta, por_beta, so_beta, re_beta, fvf_beta, chance_beta)

    # #Parameters Simulation for Probabilistic Analysis
    area_top_alpha = np.pi*(rtop_alpha**2) # in m^2
    area_base_alpha = np.pi*(rbase_alpha**2) # in m^2
    area_top_beta = np.pi*(rtop_beta**2) # in m^2
    area_base_beta = np.pi*(rbase_beta**2) # in m^2

    plt.subplot(121)
    dist_area_top_alpha = simulation_pert(area_top_alpha-area_top_alpha*20/100, area_top_alpha,
                                        area_top_alpha+area_top_alpha*30/100)
    plt.title('Area of Top Alpha Distribution')
    plt.xlabel('Area[m^2]')
    plt.ylabel('Probability Density Function')

    plt.subplot(122)
    dist_area_base_alpha = simulation_pert(area_base_alpha-area_base_alpha*20/100, area_base_alpha,
                                        area_base_alpha+area_base_alpha*50/100)
    plt.title('Area of Base Alpha Distribution')
    plt.xlabel('Area[m^2]')
    plt.ylabel('Probability Density Function')
    plt.suptitle('AREA OF ALPHA DISTRIBUTION', fontsize=20, fontweight='bold')
    plt.show()

    #########################################################################################################################################################################

    plt.subplot(121)
    dist_area_top_beta = simulation_pert(area_top_beta-area_top_beta*20/100, area_top_beta,
                                        area_top_beta+area_top_beta*30/100)
    plt.title('Area of Top Beta Distribution')
    plt.xlabel('Area[m^2]')
    plt.ylabel('Probability Density Function')

    plt.subplot(122)
    dist_area_base_beta = simulation_pert(area_base_beta-area_base_beta*20/100, area_base_beta,
                                        area_base_beta+area_base_beta*50/100)
    plt.title('Area of Top Beta Distribution')
    plt.xlabel('Area[m^2]')
    plt.ylabel('Probability Density Function')
    plt.suptitle('AREA OF BETA DISTRIBUTION', fontsize=20, fontweight='bold')

    plt.show()

    #########################################################################################################################################################################

    dist_reservoir_thickness = simulation_pert(reservoir_thickness-reservoir_thickness*20/100, reservoir_thickness,
                                            reservoir_thickness+reservoir_thickness*20/100)
    plt.title('RESERVOIR THICKNESS DISTRIBUTION', fontsize=20, fontweight='bold')
    plt.xlabel('Thickness [m]')
    plt.ylabel('Probability Density Function')
    plt.show()

    #########################################################################################################################################################################

    plt.subplot(121)
    volumetric_prob_alpha = volumetric_probabilistic(dist_area_top_alpha,dist_area_base_alpha,htop_alpha,dist_reservoir_thickness,ntg_alpha,
                                            por_alpha,so_alpha,re_alpha,chance_alpha,fvf_alpha)
    plt.title("ALPHA")
    plt.xlabel('Recoverable Oil With Risk (MMBO)')
    plt.ylabel('Probability Density Function')

    plt.subplot(122)
    volumetric_prob_beta = volumetric_probabilistic(dist_area_top_beta,dist_area_base_beta,htop_beta,dist_reservoir_thickness,ntg_beta,
                                            por_beta,so_beta,re_beta,chance_beta,fvf_beta)
    plt.title("BETA")
    plt.xlabel('Recoverable Oil With Risk (MMBO)')
    plt.ylabel('Probability Density Function')

    plt.suptitle('RECOVERABLE OIL WITH RISK DISTRIBUTION', fontsize=20, fontweight='bold')

    plt.show()

    #########################################################################################################################################################################

    ####ECONOMICS USING BREAK EVEN ANALYSIS####

    #Break Even Analysis of Alpha & Beta Deterministic
    plt.subplot(121)
    bea_alpha_det = break_even_analysis(volumetric_det_alpha,capex,opex,oilprice,project_time,tax)
    plt.title('Break-Even Analysis of Alpha', fontsize=15)

    plt.subplot(122)
    bea_beta_det = break_even_analysis(volumetric_det_beta,capex,opex,oilprice,project_time,tax)
    plt.title('Break-Even Analysis of Beta', fontsize=15)

    plt.suptitle('Break-Even Analysis of Deterministic Calculation', fontweight='bold', fontsize='20')
    plt.show()

    #########################################################################################################################################################################
    #Break Even Analysis of Alpha & Beta Probabilistic
    plt.subplot(121)
    bea_alpha_prob = break_even_analysis(volumetric_prob_alpha,capex,opex,oilprice,project_time,tax)
    plt.title('Break-Even Analysis of Alpha', fontsize=15)

    plt.subplot(122)
    bea_beta_prob = break_even_analysis(volumetric_prob_beta,capex,opex,oilprice,project_time,tax)
    plt.title('Break-Even Analysis of Beta', fontsize=15)

    plt.suptitle('Break-Even Analysis of Probabilistik Calculation', fontweight='bold', fontsize='20')
    plt.show()

    #########################################################################################################################################################################
    # SENSITIVITY DAN TORNADO
    mean_resource = round(np.mean(volumetric_det_alpha),2)

    sensi_area_top_alpha = sensitivity_calculation(dist_area_top_alpha, mean_resource)
    sensi_area_base_alpha = sensitivity_calculation(dist_area_base_alpha, mean_resource)

    sensi_area_top_beta = sensitivity_calculation(dist_area_top_beta, mean_resource)
    sensi_area_base_beta = sensitivity_calculation(dist_area_base_beta, mean_resource)

    sensi_rc = sensitivity_calculation(dist_reservoir_thickness, mean_resource)


    values_alpha = [round(sensi_area_top_alpha["value"], 2), round(sensi_area_base_alpha["value"], 2),
                    round(sensi_rc["value"], 2)]

    values_beta = [round(sensi_area_top_beta["value"], 2), round(sensi_area_base_beta["value"], 2),
                    round(sensi_rc["value"], 2)]

    lows_alpha = []
    lows_beta = []
    base = 0

    for i in range(len(values_alpha)):
        lows1 = base - values_alpha[i] / 2
        lows_alpha.append(lows1)

    for i in range(len(values_beta)):
        lows1 = base - values_beta[i] / 2
        lows_beta.append(lows1)

    # PLOTTING TORNADO
    variables_alpha = ['Area of Top Alpha',
                    'Area of Base Alpha',
                    'Reservoir Thickness']

    variables_beta = ['Area of Top Beta',
                    'Area of Base Beta',
                    'Reservoir Thickness']

    ###########################################################
    sorted_alpha = sorting_sensitivity(values_alpha, variables_alpha, lows_alpha)
    sorted_beta = sorting_sensitivity(values_beta, variables_beta, lows_beta)

    values_alpha = sorted_alpha[0]
    values_beta = sorted_beta[0]
    variables_alpha = sorted_alpha[1]
    variables_beta = sorted_beta[1]
    lows_alpha = sorted_alpha[2]
    lows_beta = sorted_beta[2]

    # The actual drawing part

    # The y position for each variable
    ys_alpha = range(len(values_alpha))[::-1]  # top to bottom
    ys_beta = range(len(values_beta))[::-1]  # top to bottom

    plt.subplot(121)

    # Plot the bars, one by one
    plot_bars(ys_alpha, lows_alpha, values_alpha, base)

    # Draw a vertical line down the middle
    plt.axvline(base, color='black')

    # Position the x-axis on the top, hide all the other spines (=axis lines)
    axes = plt.gca()  # (gca = get current axes)
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.xaxis.set_ticks_position('top')

    # Make the y-axis display the variables
    plt.yticks(ys_alpha, variables_alpha)

    # Set the portion of the x- and y-axes to show
    plt.xlim(base - 100, base + 100)
    plt.ylim(-1, len(variables_alpha))
    plt.title('ALPHA')

    plt.subplot(122)

    # Plot the bars, one by one
    plot_bars(ys_beta, lows_beta, values_beta, base)

    # Draw a vertical line down the middle
    plt.axvline(base, color='black')

    # Position the x-axis on the top, hide all the other spines (=axis lines)
    axes = plt.gca()  # (gca = get current axes)
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.xaxis.set_ticks_position('top')

    # Make the y-axis display the variables
    plt.yticks(ys_alpha, variables_beta)

    # Set the portion of the x- and y-axes to show
    plt.xlim(base - 100, base + 100)
    plt.ylim(-1, len(variables_beta))
    plt.title('BETA')
    plt.suptitle("SENSITIVITY CHART", fontsize=20, fontweight='bold')
    plt.show()

    #########################################################################################################################################################################

    gatherdata_parameters = {'':["GRV [acre.ft]", "Net to Gross", "Porosity", "Oil Saturation", "Recovery Effeciency",
                                "Formation Volume Factor", "Chance Oil Trapped"],
                            'ALPHA': [grv_alpha, ntg_alpha, por_alpha, so_alpha, re_alpha, fvf_alpha, chance_alpha],
                            'BETA': [grv_beta, ntg_beta, por_beta, so_beta, re_beta, fvf_beta, chance_beta]}

    gatherdata_volumetric = {'':["Volumetric Deterministic", "Volumetric Probabilistic", "Payout Period Deterministik",
                                "Total Profit Deterministic"],
                            'ALPHA': [volumetric_det_alpha, volumetric_prob_alpha["mean"], bea_alpha_det["payout"],
                                    bea_alpha_det["total profit"]],
                            'BETA': [volumetric_det_beta, volumetric_prob_beta["mean"], bea_beta_det["payout"],
                                    bea_beta_det["total profit"]]}

    df_par = pd.DataFrame(gatherdata_parameters)
    df_vol = pd.DataFrame(gatherdata_volumetric)

    print(df_par)
    print(df_vol)

def model_3d_main():
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

def main():
    monte_carlo_main()
    model_3d_main()

if __name__ == "__main__":
    main()
