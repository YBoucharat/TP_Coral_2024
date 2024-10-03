import os
import sys
#sys.path.insert(0, os.path.abspath('./Library'))
sys.path.insert(0, os.path.abspath('./'))
import numpy as np
import xsimlab as xs
from reef_models import reef, reef_platform
from tools import nlines
from Dicts import Dicos
from Dict_models import DicoModels
#from tools_models import CheckSimu, ZarrName, shore
from datetime import datetime
#import dask
#import cProfile

# Set the number of CPU cores to be used
#dask.config.set(scheduler='threads', num_workers=4)

dt = 100
dico = Dicos()

#SL = 'Waelbroeck2002b'
#SL = 'Spratt2016-430k'
#SL = 'Bintanja2008-1500k'
#SL = 'Grant2014-100k'
# SL = 'Rohling2009-430k'
SL = 'Holocene-SL.dat'

tmax = nlines(dico.SL_files['Holocene'])*1e3

ds_in = xs.create_setup(
    model=reef,
    clocks={
        'time' : np.arange(0., tmax+dt, dt)
},
    master_clock = 'time',
    input_vars={
        ## Usefull parameters to test. Everything in meters and years.
        # vertical land motion rate. [-10 - 2]e-3 (huge values, you can test before with 0.something e-3)
        'vertical__u': 0.2e-3,
        'construct__hmax': 20,
        # maximum water height for reef growth. [10 - 50] 
        'grid__dmax': 20,
        # initial slope of the substrate. [2 - 8]
        'grid__slopi': 2e-2,
        # maximum reef growth rate. [2 - 15]
        'construct__Gm': 10e-3,
        # width of wave action (controls the width of the coral reef). [10 - 3000] 
        'hfactor__Dbar': 100,
        # Eroded volume. [50 - 1000]e-3
        'eros__Ev': 200e-3,
        # Water height for wave base (the waves will start to erode at this depth). [1 - 10]
        'eros__hwb': 3.,
        # # Elevation of antecedent terrace. As you want, depends on vertical__u.
        # 'init__zterr': -20,
        # # Length of antecedent terrace. Put it to 0 if no platform, 10000 is already big
        # 'init__lterr': 10000,
        ## No need to test 
        # Water height for open ocean
        'hfactor__how': 2,
        # uniform spacing
        'grid__spacing': 1,
        # filename for RSL reconstruction
        'SLstory__RSLin': SL,
        # Coefficient for erosion efficiency, sea-bed
        'eros__beta1': 0.1,
        # Coefficient for erosion efficiency, cliff retreat
        'eros__beta2': 1,
        # Height of notch for volume eroded during cliff retreat
        'eros__hnotch': 1,
        # ---
        'depot__repos': 15e-2,
    },
    output_vars={
        'init__x'       : None,
        'profile__z'    : None,
        'sealevel__asl' : None,
        'profile__xmin' : None,
        'profile__xmax' : None,
    }
)

ds_in.attrs['model_name'] = 'reef'
ds_in.attrs['store'] = 'LastProfile'

dm   = DicoModels()
tstart = datetime.now()
print(tstart)

with dm.models[ds_in.model_name]:
    ds_out = (ds_in
      .xsimlab.run()
              )

#x = ds_out.x[:].values
#y = ds_out.profile__z[:].values

#data = np.column_stack((x, y))
#np.savetxt('topo_holo.dat', data, delimiter='\t', fmt='%.1f')

#with dm.models[ds_in.model_name]:ds_out = cProfile.run('ds_in.xsimlab.run()')

print("duration :", datetime.now()-tstart)
