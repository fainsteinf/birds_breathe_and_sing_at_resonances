# Data and code for "Birds breathe and sing at resonances of the biomechanics"

# F. Fainstein (1,2), F. Goller (3,4), and G. B. Mindlin (1,2)

# 1. Universidad de Buenos Aires, Facultad de Ciencias Exactas y Naturales, Departamento de Física, Ciudad Universitaria, 1428 Buenos Aires, Argentina.
# 2. CONICET - Universidad de Buenos Aires, Instituto de Física Interdisciplinaria y Aplicada (INFINA), Ciudad Universitaria, 1428 Buenos Aires, Argentina.
# 3. Institute of Integrative Cell Biology and Physiology, University of Münster, Münster 48143, Germany.
# 4. School of Biological Sciences, University of Utah, Salt Lake City, Utah 84112, USA.


#Grid search fitting procedure of song data. 


import numpy as np
from scipy.interpolate import UnivariateSpline
import time as time_clock
from multiprocessing import Pool, cpu_count
import os
import pandas as pd
from tqdm.notebook import tqdm

#integrador
def rk4(dxdt, x, t, dt, *args, **kwargs):
    x = np.asarray(x)
    k1 = np.asarray(dxdt(x, t, *args, **kwargs))*dt
    k2 = np.asarray(dxdt(x + k1*0.5, t, *args, **kwargs))*dt
    k3 = np.asarray(dxdt(x + k2*0.5, t, *args, **kwargs))*dt
    k4 = np.asarray(dxdt(x + k3, t, *args, **kwargs))*dt
    return x + (k1 + 2*k2 + 2*k3 + k4)/6

def sigm(dp):
    return (np.sign(dp)+1)/2

def f_model_singular(v, t, pars):
    
    x = v[0]
    p = v[1]
    
    (taux, taup, alfa_in, alfa_rel, F0, force) = [par for par in pars]
    
    alfa = alfa_in * ( 1 - ( 1 - alfa_rel ) *  sigm(p) ) 
    
    dxdt = (1/taux) * ( - (1 + 1 * x**2) * x - p + F0 * force )
    
    dpdt = (1/taup) * ( - (1 + 1 * x**2) * x - (1 + alfa) * p + F0 * force )
    
    return [dxdt, dpdt]

######## Load data ########
bird_names = ["CaFF016-VioVio", "CaFF-NeVe", "CaFF-BlaVe", "CaFF073-RoVio", 
              "CaFF909-NaRo"]

# Change to analyze different birds
bird_number = 5   #bird numbers are [1, 2, 3, 4, 5]

#Get path to data (or add manually)
current_path = os.path.dirname(os.path.abspath(__file__))
current_path = current_path[:current_path.rfind("\\")]
path_data = current_path + "/" + "Data/bird" + str(bird_number) + "/" 

#Load data
df = pd.read_csv(path_data+"data_"+bird_names[bird_number-1]+".csv", sep=";")  
time = df["time"].to_numpy()
force = df["force"].to_numpy()
p_day = df["pressure"].to_numpy()

######## Resample data ########
dt = 3 * 10**-4
time_int = np.arange(time[0], time[-1], dt)
force_int = UnivariateSpline(time, force, s=0, k=1)(time_int)
p_day_int = UnivariateSpline(time, p_day, s=0, k=1)(time_int)

# === Parameter grid (to be changed) ===
tauxs = np.linspace(0.025, 0.1, 2)
taups = np.linspace(0.005, 0.1, 2)
alfa_ins = np.linspace(.1, 3, 2)
alfa_rels = np.linspace(0.025, 0.125, 2)
F0s = np.linspace(20, 100, 2)

param_grid = [
    (taux, taup, alfa_in, alfa_rel, F0)
    for taux in tauxs
    for taup in taups
    for alfa_in in alfa_ins
    for alfa_rel in alfa_rels
    for F0 in F0s
]
#%%
#### Fitting function to paralellize ### 
def evaluate_combination(args):
    taux, taup, alfa_in, alfa_rel, F0 = args
    pars = [taux, taup, alfa_in, alfa_rel, F0, force[0]]

    p = np.zeros_like(time_int)
    x = np.zeros_like(time_int)

    for ix, tt in enumerate(time_int[:-1]):
        pars[-1] = force_int[ix]
        x[ix+1], p[ix+1] = rk4(f_model_singular, [x[ix], p[ix]], tt, dt, pars)

    chisqr_red = np.sum(np.power(p_day_int - p, 2)) / (len(p) - 5)
    return [taux, taup, alfa_in, alfa_rel, F0, chisqr_red]

#%%# Paralellized fitting algorithm

if __name__ == "__main__":
    start_time = time_clock.time()
    
    parameters_chisqr = []
    save_every = 1000  # Save every 1000 iterations
    save_folder = current_path+"/"
    base_name = "fit_song_bird"+str(bird_number)
    final_path = os.path.join(save_folder, base_name + ".npy")

    with Pool(cpu_count()) as pool:
        for i, result in enumerate(tqdm(pool.imap_unordered(evaluate_combination, param_grid), total=len(param_grid))):
            parameters_chisqr.append(result)

            if (i + 1) % save_every == 0 or (i + 1) == len(param_grid):
                partial_path = os.path.join(save_folder, f"{base_name}_parcial_{i+1}.npy")
                np.save(partial_path, np.array(parameters_chisqr))
                print(f"Partial save: {partial_path}")
    
    # Save final result
    np.save(final_path, np.array(parameters_chisqr))
    print("--- Total time: %.2f seconds ---" % (time_clock.time() - start_time))
