# Data and code for "Birds breathe and sing at resonances of the biomechanics"

# F. Fainstein (1,2), F. Goller (3,4), and G. B. Mindlin (1,2)

# 1. Universidad de Buenos Aires, Facultad de Ciencias Exactas y Naturales, Departamento de Física, Ciudad Universitaria, 1428 Buenos Aires, Argentina.
# 2. CONICET - Universidad de Buenos Aires, Instituto de Física Interdisciplinaria y Aplicada (INFINA), Ciudad Universitaria, 1428 Buenos Aires, Argentina.
# 3. Institute of Integrative Cell Biology and Physiology, University of Münster, Münster 48143, Germany.
# 4. School of Biological Sciences, University of Utah, Salt Lake City, Utah 84112, USA.


#Code to integrate the dynamical equations for each bird, during quiet breathing.

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from scipy import signal
import os
import pandas as pd

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def get_spectrogram(data, sampling_rate, window=1024, overlap=1/1.1,
                    sigma=102.4, fmax=8000, drange=6):
    
    #Create spectrogram 
    
    fu, tu, Sxx = signal.spectrogram(data, sampling_rate, nperseg=window,
                                     noverlap=window*overlap,
                                     window=signal.get_window
                                     (('gaussian', sigma), window),
                                     scaling='spectrum')
    Sxx = np.clip(Sxx, a_min=np.amax(Sxx)/10**drange, a_max=np.amax(Sxx))
    return fu, tu, Sxx

def rk4(dxdt, x, t, dt, *args, **kwargs):   
    x = np.asarray(x)
    k1 = np.asarray(dxdt(x, t, *args, **kwargs))*dt
    k2 = np.asarray(dxdt(x + k1*0.5, t, *args, **kwargs))*dt
    k3 = np.asarray(dxdt(x + k2*0.5, t, *args, **kwargs))*dt
    k4 = np.asarray(dxdt(x + k3, t, *args, **kwargs))*dt
    return x + (k1 + 2*k2 + 2*k3 + k4)/6

def sigm(dp):
    #Heaviside function
    return (np.sign(dp)+1)/2

def f_model_singular(v, t, pars):
    #Model equations
    
    x = v[0]
    p = v[1]
    
    (taux, taup, alfa_in, alfa_out, F0, force) = [par for par in pars]
    
    alfa =   alfa_in - ( alfa_in - alfa_out ) *  sigm(p) 
    
    dxdt = (1/taux) * ( - (1 + 1 * x**2) * x - p + F0 * force )
    
    dpdt = (1/taup) * ( - (1 + 1 * x**2) * x - (1 + alfa) * p + F0 * force )
    
    return [dxdt, dpdt]

#%% Get data

bird_names = ["CaFF016-VioVio", "CaFF-NeVe", "CaFF-BlaVe", "CaFF073-RoVio", 
              "CaFF909-NaRo"]

# Change to analyze different birds
bird_number = 1   #bird numbers are [1, 2, 3, 4, 5]

#Get path to data (or add manually)
current_path = os.path.dirname(os.path.abspath(__file__))
current_path = current_path[:current_path.rfind("\\")]
path_data = current_path + "/" + "Data/bird" + str(bird_number) + "/" 

#Load data
df = pd.read_csv(path_data+"data_quiet_"+bird_names[bird_number-1]+".csv", sep=";")  
time = df["time"].to_numpy()
force = df["force"].to_numpy()
p_day = df["pressure"].to_numpy()

fs = 1 / np.mean(np.diff(time))

# Load parameters
df_pars = pd.read_csv(current_path + "/" + "Data/" + "pars_quiet.csv", sep=";")  
values_bird = df_pars[df_pars["birdID"] == bird_names[bird_number-1]].iloc[0]
pars = values_bird[["taux", "taup", "alpha_in", "alpha_out", "F0"]].to_numpy()
pars = np.concatenate((pars, [0]))

#Integrate the model
dt = 1 / fs
p = np.zeros_like(time)
x = np.zeros_like(time)

p[0] = p_day[0]
for ix, tt in enumerate(time[:-1]):
    
    pars[-1] = force[ix]

    x[ix+1], p[ix+1] = rk4(f_model_singular, [x[ix], p[ix]], tt, dt,
                                               pars)

#%% Plot results

fsize=6
lwi = .5

fig, ax = plt.subplots(2, 1, figsize=(3.37, 1), sharex=True)
ax[0].plot(time, p_day, color='k', lw=lwi)
ax[0].set_ylabel("Pressure \n(arb. u.)", fontsize=fsize)

ax[1].plot(time, p, lw=lwi)
ax[1].set_ylabel("Pressure \n(arb. u.)", fontsize=fsize)

for k in range(2):
    ax[k].spines["right"].set_visible(False)
    ax[k].spines["top"].set_visible(False)
    ax[k].tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax[k].spines["bottom"].set_visible(False)
    ax[k].tick_params(axis='both', which='both', labelsize=fsize)
    ax[k].get_yaxis().set_label_coords(-0.1,0.5)
    ax[k].set_yticks([-1, 0, 1])
fig.tight_layout()
ax[0].set_xlim([0.,time[-1]])

fig.subplots_adjust(bottom=0.05, top=0.93, hspace=0.25) 
plt.show()