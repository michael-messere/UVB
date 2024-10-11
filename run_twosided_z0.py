# 

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable, Table, Column
from astropy import units as u
from astropy.io import ascii
import pandas as pd

import pyCloudy as pc
pc.log_.level = 3

pc.config.cloudy_exe = '/Users/messeremichael/Documents/Software/c23.01/source/cloudy.exe'

for i, zheight in enumerate(np.arange(0,5050,50)):

    dir_ = 'runs_twosided/power_law_density_30kpc_disk' + '_' + str(zheight) + '/'
    model_name = 'test'  
    full_model_name = '{0}{1}'.format(dir_, model_name)
    
    print(full_model_name)

    c_input = pc.CloudyInput(full_model_name)

    pc.log_.message('Running {0}'.format(model_name), calling = 'test1')
    pc.log_.timer('Starting Cloudy', quiet = True, calling = 'test1')

    c_input.run_cloudy()
    pc.log_.timer('Cloudy ended after seconds:', calling = 'test1')
