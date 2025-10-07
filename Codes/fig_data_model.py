# Data and code for "Birds breathe and sing at resonances of the biomechanics"

# F. Fainstein (1,2), F. Goller (3,4), and G. B. Mindlin (1,2)

# 1. Universidad de Buenos Aires, Facultad de Ciencias Exactas y Naturales, Departamento de Física, Ciudad Universitaria, 1428 Buenos Aires, Argentina.
# 2. CONICET - Universidad de Buenos Aires, Instituto de Física Interdisciplinaria y Aplicada (INFINA), Ciudad Universitaria, 1428 Buenos Aires, Argentina.
# 3. Institute of Integrative Cell Biology and Physiology, University of Münster, Münster 48143, Germany.
# 4. School of Biological Sciences, University of Utah, Salt Lake City, Utah 84112, USA.


#Code to integrate the dynamical equations for each bird. 
#Make figures (Fig. 2 and Fig. S1)

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
df = pd.read_csv(path_data+"data_"+bird_names[bird_number-1]+".csv", sep=";")  
time = df["time"].to_numpy()
m_day = df["emg"].to_numpy()
force = df["force"].to_numpy()
p_day = df["pressure"].to_numpy()
s_day = df["sound"].to_numpy()

#Sampling rate
fs = 44150

#filter sound
order = 5
cutoff_1, cutoff_2 = 800, 15000
b, a = signal.butter(order, [cutoff_1, cutoff_2],btype='bandpass', fs=fs)
s_day = signal.filtfilt(b, a, s_day)

# Load parameters
df_pars = pd.read_csv(current_path + "/" + "Data/" + "pars_song.csv", sep=";")  
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
msize = 3
box_width = .5

fig, ax = plt.subplots(4,1,figsize=(6.7,2.5),sharex=True, height_ratios=[.75, .75, 1, 1])

ax[0].tick_params(axis='x',which='both', bottom=False, top=False, labelbottom=False)

fu, tu, Sxx = get_spectrogram(s_day, 44150)
Sxx = np.clip(Sxx, a_min=np.amax(Sxx)/4000, a_max = np.amax(Sxx))

ax[0].pcolormesh(tu, fu/1000, np.log(Sxx), rasterized=True, shading='auto', 
                  cmap='binary')
ax[0].set_ylim([0,15])
ax[0].set_yticks([0, 15])
ax[1].plot(time, abs(m_day),'k', lw=lwi*.25)

ax[2].plot(time, p_day, '-', color='k', lw=lwi, zorder=-1, label="Data")
ax[3].plot(time[int(0.06*fs):], p[int(0.06*fs):], label="Model", lw=lwi*1.25)
ax[2].axhline(0, color='k', lw=lwi*.25, alpha=1, zorder=-2)
ax[3].axhline(0, color='k', lw=lwi*.25, alpha=1, zorder=-2)

ax[0].set_ylabel("Frequency \n(kHz)", fontsize=fsize)
ax[1].set_ylabel("EMG \n(V)", fontsize=fsize)
ax[2].set_ylabel("Pressure \n(arb. units)", fontsize=fsize)
ax[3].set_ylabel("Pressure \n(arb. units)", fontsize=fsize)

ax[1].set_yticks([0, 0.2])
ax[2].set_yticks([-10, 30])
ax[3].set_yticks([-10, 30])

for k in range(4):
    if k>0:
        ax[k].spines["right"].set_visible(False)
        ax[k].spines["top"].set_visible(False)
        ax[k].tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        ax[k].spines["bottom"].set_visible(False)
    ax[k].tick_params(axis='both', which='both', labelsize=fsize, width=box_width)

    ax[k].get_yaxis().set_label_coords(-0.05,0.5)
    
    for axis in ['top','bottom','left','right']:
        ax[k].spines[axis].set_linewidth(box_width)
        
fig.tight_layout()
ax[0].set_xlim([0.,time[-1]])

fig.subplots_adjust(bottom=0.05, top=0.93, hspace=0.25) # or whatever
plt.show()