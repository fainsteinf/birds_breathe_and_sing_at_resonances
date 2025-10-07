# Data and code for "Birds breathe and sing at resonances of the biomechanics"

# F. Fainstein (1,2), F. Goller (3,4), and G. B. Mindlin (1,2)

# 1. Universidad de Buenos Aires, Facultad de Ciencias Exactas y Naturales, Departamento de Física, Ciudad Universitaria, 1428 Buenos Aires, Argentina.
# 2. CONICET - Universidad de Buenos Aires, Instituto de Física Interdisciplinaria y Aplicada (INFINA), Ciudad Universitaria, 1428 Buenos Aires, Argentina.
# 3. Institute of Integrative Cell Biology and Physiology, University of Münster, Münster 48143, Germany.
# 4. School of Biological Sciences, University of Utah, Salt Lake City, Utah 84112, USA.


#Magnification curves and rate distributions figures (Fig. 3 and Fig. S4)

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from scipy.stats import gaussian_kde
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
import os
import pandas as pd

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

bird_names = ['CaFF016-VioVio', 'CaFF-NeVe', 'CaFF-BlaVe', 'CaFF073-RoVio',  
              'CaFF909-NaRo']

# Change to analyze different birds
bird_number = 1   #bird numbers are [1, 2, 3, 4, 5]

####################### Get results ##################### 
#Get path to data  (or add manually)
current_path = os.path.dirname(os.path.abspath(__file__))
current_path = current_path[:current_path.rfind("\\")]
path_data = current_path + "/" + "Data/bird" + str(bird_number) + "/" 

# Get song magnification
df_aux = pd.read_csv(path_data+"magnification_song_"+bird_names[bird_number-1]+".csv", sep=";")  
frequency_s = df_aux["frequency"].to_numpy()
magnification_s = df_aux["magnification"].to_numpy()

# Get quiet breathing magnification
df_aux = pd.read_csv(path_data+"magnification_quiet_"+bird_names[bird_number-1]+".csv", sep=";")  
frequency_q = df_aux["frequency"].to_numpy()
magnification_q = df_aux["magnification"].to_numpy()

# Get song breathing rate and labels of syllable clusters
df_aux = pd.read_csv(path_data+"song_rate_labels_"+bird_names[bird_number-1]+".csv", sep=";")  
period_s = df_aux["period"].to_numpy()
labels = df_aux["label"].to_numpy()

# Get quiet breathing rate
df_aux = pd.read_csv(path_data+"quiet_rate_"+bird_names[bird_number-1]+".csv", sep=";")  
period_q = df_aux["period"].to_numpy()

#%% Computations prior to figure making

#Compute number of song syllable clusters
clustered = (labels >= 0)
if -1 in np.unique(labels):    
    n_clusters = len(np.unique(labels))-1
else:
    n_clusters = len(np.unique(labels))

#Interpolation of the magnification curves. To plot mode syllabic rates.  
spline_s = UnivariateSpline(frequency_s, magnification_s, s=0, k=3)
spline_q = UnivariateSpline(frequency_q, magnification_q, s=0, k=3)

#Gaussian KDE of the distribution of syllabic rates given a syllable type
x_vals_clusters = []
y_vals_clusters = []
for i in range(n_clusters):
    indices = labels == i
    data_k = 1 / period_s[indices]
    kde = gaussian_kde(data_k)
    med = np.mean(data_k)
    std = np.std(data_k)
    x_vals = np.linspace(med-4*std, med+4*std, 500)
    y_vals = kde(x_vals)
    y_vals_clusters.append(y_vals)
    x_vals_clusters.append(x_vals)

#%% Plot figure

lwi = 1
fsize = 6
box_width = .5
msize = 5

cmap = matplotlib.colormaps['tab20']
colors = cmap(np.linspace(0, 1, n_clusters))

fig, axs = plt.subplots(3, 1, figsize=(3.37 / 2, 2.5), 
                        height_ratios=[.75, .25, .25], sharex=True)

#Plot magnification curves
axs[0].plot(frequency_s, magnification_s, color='k', 
               lw = lwi)
axs[0].plot(frequency_q, magnification_q, 'Grey', 
               lw = lwi, linestyle='dashed', alpha=.7)

axs[0].text(0.87, 0.9, "bird "+str(bird_number), horizontalalignment='center',
     verticalalignment='center', transform=axs[0].transAxes, fontsize=fsize)

#Plot histograms 
for i in range(n_clusters):
    axs[1].plot(x_vals_clusters[i], y_vals_clusters[i],color=colors[i], lw=lwi, 
                zorder=i)
    
    #Graficar histograma tambien
    indices = labels == i
    data_k = 1 / period_s
    nbins = int(np.sqrt(np.sum(indices)))
    counts, bins = np.histogram(1 / period_s[indices], bins=nbins, 
                                density=True)
    bin_centers = .5 * (bins[:-1] + bins[1:])
    bin_width = bins[1] - bins[0]

    axs[1].bar(bin_centers, counts, width=bin_width, bottom=-i*.0, color=colors[i],
            alpha=0.8, zorder=i)
    
    #Plot modes of distributions
    peaks, _ = find_peaks(y_vals_clusters[i], height=0.1, prominence=0.05)
    for peak in peaks:
        maximum = x_vals_clusters[i][peak]
        axs[0].scatter(maximum, spline_s(maximum), color=colors[i], 
                  s=msize, zorder=10)
        
#Histogram of song and quiet breathing rates
ax_aux = axs[2].twinx()
ax_aux.hist(1/period_q, 100, edgecolor='k', alpha=.75,histtype='stepfilled',
            density=True, zorder=0, facecolor='Grey', lw=lwi/4)
ax_aux.set_yticks([0, 2])
ax_aux.tick_params(axis='both', which='both', labelsize=fsize, width=0.5)

histsong = axs[2].hist(1/period_s, int(np.sqrt(len(period_s))*4), edgecolor='k', alpha=1,histtype='stepfilled', 
            facecolor='k', density=True, lw=lwi/4)

#Set axis limits
if bird_number in [3, 5]:
    axs[0].set_xlim([0, 30])
    axs[0].set_xticks([0, 10, 20, 30])
else:
    axs[0].set_xlim([0, 25])
    axs[0].set_xticks(np.arange(0, 30, 5))
axs[0].set_ylim([.4, 1.05])
axs[1].set_yticks([0, np.round(np.max(y_vals_clusters)*1.4, 1)])
axs[2].set_yticks([0, np.round(np.max(histsong[0])*1.5, 1)])

axs[0].set_ylabel("Pressure response", fontsize=fsize)
axs[1].set_ylabel("Density", fontsize=fsize)
axs[2].set_ylabel("Density", fontsize=fsize)
axs[2].set_xlabel("Frequency (Hz)", fontsize=fsize, labelpad=2)

for i in range(3):    
    axs[i].yaxis.set_label_coords(-0.175,0.5)
    axs[i].tick_params(axis='both', which='both', labelsize=fsize, 
                       width=box_width)
    for axis in ['top','bottom','left','right']:
        axs[i].spines[axis].set_linewidth(box_width)
for axis in ['top','bottom','left','right']:
    ax_aux.spines[axis].set_linewidth(0)

fig.tight_layout()
fig.subplots_adjust(bottom=0.125, top=0.93, hspace=0.25, 
                    left=0.2, right=0.9) 
plt.show()

