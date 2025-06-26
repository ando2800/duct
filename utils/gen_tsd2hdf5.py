import os
import numpy as np
import pandas as pd
from io import StringIO
import subprocess


print('[MESSAGE] Create time series data from historical log of duct simulation.')
#==== Constants and parameters (refer to simulation of duct & runparams2.dat)  === #
CMASSFLOW = 1.92783343620839 
h = 1.0
ub = CMASSFLOW / 4.0 
nstep = 140001
nhist = 10
dt = 2.0e-2
dt_nondim = dt * ub/h 

inputfile = "./data/raw/output.txt"
outputfile =  "./data/processed/statistics_tsd.h5"


#==== Functions to calculate new statistics ====#
def calc_I(sectors):
    Sn, Sw, Ss, Se = sectors
    total = Sn + Sw + Ss + Se
    if total == 0:
        return 0
    return (Sn + Ss - Sw - Se) / total

def calc_S(sectors):
    Sn, Sw, Ss, Se = sectors
    return (Sn + Sw + Ss + Se) * ub / h

# ==== 1: Fetch statistics from duct simulation log  ==== #
# grep + awk processe 
grep = subprocess.Popen(["grep", "ENER", inputfile], stdout=subprocess.PIPE, text=True)
awk = subprocess.Popen(["awk", "{print $14, $15, $16, $17}"],
                       stdin=grep.stdout, stdout=subprocess.PIPE, text=True)
output, _ = awk.communicate()
# read by pandas
df = pd.read_csv(StringIO(output), sep=r'\s+', header=None)
df.columns = ['Sn', 'Sw', 'Ss', 'Se'] 

# ==== 2: Calculate non-dimensional time #
df['time_nondim'] =  dt_nondim + 10*dt_nondim * np.arange(1 + (nstep-1)/nhist)

# ==== 3: Calculate new statistics ==== #
# dimentionless indicator function(Uhlmann et al. 2007)
df['I'] = df[['Sn', 'Sw', 'Ss', 'Se']].apply(lambda row: calc_I(row), axis=1)
# Sum of the enstrophy contained in the streamwise-averaged vorticity
df['Ssum'] = df[['Sn', 'Sw', 'Ss', 'Se']].apply(lambda row: calc_S(row), axis=1)


# ==== Save to hdf5 ====
with pd.HDFStore(outputfile, mode='w') as store:
    store.put('/statistics/Sn', df['Sn'])
    store.put('/statistics/Sw', df['Sw'])
    store.put('/statistics/Ss', df['Ss'])
    store.put('/statistics/Se', df['Se'])
    store.put('/new_statistics/Indicator', df['I'])
    store.put('/new_statistics/Ssum', df['Ssum'])
    store.put('/time_nondim', df['time_nondim'])
print("[MESSAGE] Time series data saved to HDF5 format successfully!")