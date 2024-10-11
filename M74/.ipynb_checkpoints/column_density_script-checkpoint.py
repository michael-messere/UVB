import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
from scipy.interpolate import interp1d
from scipy.integrate import trapz

params = {'mathtext.default': 'regular' } 
plt.rcParams.update(params)
plt.rcParams['font.family'] = 'serif'
from matplotlib.gridspec import GridSpec

import pyCloudy as pc
pc.log_.level = 3

########################################################

def column_density(path):
    
    zheight_max = 4000
    zheight_delta = 250
    
    lower_interp_r = 0
    upper_interp_r = 95
    
    lower_interp_z = 0
    upper_interp_z = 3.99
    
    zheight = np.arange(0,zheight_max+zheight_delta,zheight_delta)
    model = []

    for i in range(len(zheight)):

        dir_ = path + '_' + str(zheight[i]) + '/'
        model_name = 'test'  
        full_model_name = '{0}{1}'.format(dir_, model_name)
        model.append(pc.CloudyModel(full_model_name))
        
    H_p0_interp = []
    H_p1_interp = []
    H_interp = []
    
    N = 1000 # how many points to interpolate over

    for i in range(len(model)):

        x = np.abs(model[i].radius[::-1]/3.06e21-model[i].radius[-1]/3.06e21)
        y = model[i].get_ionic('H', 0)[::-1]*model[i].nH[::-1]
        interp_linear = interp1d(x, y, kind='linear')
        x_new = np.linspace(lower_interp_r, upper_interp_r, N) 
        y_new_linear = interp_linear(x_new)
        H_p0_interp.append(y_new_linear)

        x = np.abs(model[i].radius[::-1]/3.06e21-model[i].radius[-1]/3.06e21)
        y = model[i].get_ionic('H', 1)[::-1]*model[i].nH[::-1]
        interp_linear = interp1d(x, y, kind='linear')
        x_new = np.linspace(lower_interp_r, upper_interp_r, N)  
        y_new_linear = interp_linear(x_new)
        H_p1_interp.append(y_new_linear)

        x = np.abs(model[i].radius[::-1]/3.06e21-model[i].radius[-1]/3.06e21)
        y = model[i].nH[::-1]
        interp_linear = interp1d(x, y, kind='linear')
        x_new = np.linspace(lower_interp_r, upper_interp_r, N)  
        y_new_linear = interp_linear(x_new)
        H_interp.append(y_new_linear)
        
    H_p0_interp = np.array(H_p0_interp)
    H_p1_interp = np.array(H_p1_interp)
    H_interp = np.array(H_interp)
    
    H_p0_interp_z = []
    H_p0_column_density = []
    H_p1_interp_z = []
    H_p1_column_density = []
    H_interp_z = []
    H_column_density = []

    for i in range(N):

        x = zheight/1000
        y = H_p0_interp[:,i]
        interp_linear = interp1d(x, y, kind='linear')
        x_new = np.linspace(lower_interp_z, upper_interp_z, N) 
        y_new_linear = interp_linear(x_new)
        H_p0_interp_z.append(y_new_linear)
        H_p0_column_density.append(2*trapz(y_new_linear, x_new)*3.086e+18*zheight_max)

        x = zheight/1000
        y = H_p1_interp[:,i]
        interp_linear = interp1d(x, y, kind='linear')
        x_new = np.linspace(lower_interp_z, upper_interp_z, N) 
        y_new_linear = interp_linear(x_new)
        H_p1_interp_z.append(y_new_linear)
        H_p1_column_density.append(2*trapz(y_new_linear, x_new)*3.086e+18*zheight_max)

        x = zheight/1000
        y = H_interp[:,i]
        interp_linear = interp1d(x, y, kind='linear')
        x_new = np.linspace(lower_interp_z, upper_interp_z, N) 
        y_new_linear = interp_linear(x_new)
        H_interp_z.append(y_new_linear)
        H_column_density.append(2*trapz(y_new_linear, x_new)*3.086e+18*zheight_max)
        
    H_alpha_interp = []

    for i in range(len(model)):

        x = np.abs(model[i].radius[::-1]/3.06e21-model[i].radius[-1]/3.06e21)
        y = model[i].get_emis('H__1_656281A')[::-1]
        interp_linear = interp1d(x, y, kind='linear')
        x_new = np.linspace(lower_interp_r, upper_interp_r, N) 
        y_new_linear = interp_linear(x_new)
        H_alpha_interp.append(y_new_linear)
        
    H_alpha_interp = np.array(H_alpha_interp)
    
    H_alpha_interp_z = []
    H_alpha_total = []

    for i in range(N):

        x = zheight/1000
        y = H_alpha_interp[:,i]
        interp_linear = interp1d(x, y, kind='linear')
        x_new = np.linspace(lower_interp_z, upper_interp_z, N) 
        y_new_linear = interp_linear(x_new)
        H_alpha_interp_z.append(y_new_linear)
        H_alpha_total.append(2*trapz(y_new_linear, x_new)*3.086e+18*zheight_max)
        
    #fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20,10), frameon=False)
    #plt.subplots_adjust(wspace=0.1, hspace=0.0)

    #ax.plot(np.linspace(lower_interp_r, upper_interp_r, N),H_p0_column_density,color='black')
    #ax.plot(np.linspace(lower_interp_r, upper_interp_r, N),H_p1_column_density,color='salmon')
    #ax.plot(np.linspace(lower_interp_r, upper_interp_r, N),H_column_density,color='black')

    #ax.set_yscale('log')
    #plt.show()
    
    return np.linspace(lower_interp_r, upper_interp_r, N), H_p0_column_density, H_p1_column_density, H_column_density, H_alpha_total