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

grv_alpha = grv_calc(rtop_alpha, rbase_alpha, htop_alpha, hbase_alpha) #meter cubic
grv_beta =  grv_calc(rtop_beta, rbase_beta, htop_beta, hbase_beta) #meter cubic
