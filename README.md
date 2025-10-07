# Readme

The repository contains data and codes that support the analysis in a paper under consideration for publication.

-----------------------------------------------------------------------------------

Data and code for "Birds breathe and sing at resonances of the biomechanics"

F. Fainstein (1,2), F. Goller (3,4), and G. B. Mindlin (1,2)

1. Universidad de Buenos Aires, Facultad de Ciencias Exactas y Naturales, Departamento de Física, Ciudad Universitaria, 1428 Buenos Aires, Argentina.
2. CONICET - Universidad de Buenos Aires, Instituto de Física Interdisciplinaria y Aplicada (INFINA), Ciudad Universitaria, 1428 Buenos Aires, Argentina.
3. Institute of Integrative Cell Biology and Physiology, University of Münster, Münster 48143, Germany.
4. School of Biological Sciences, University of Utah, Salt Lake City, Utah 84112, USA.

Data was collected at the Dynamical Systems Lab, Physics Department, University of Buenos Aires, following 
the procedures described in the article and Supplemental Material. 

This repository contains the data and code supporting the results presented in the study. Below is a description of the contents of the Data and Codes directories.

📁 Data
Summary Files
pars_quiet.csv
Fitted model parameters during quiet breathing for each bird.

pars_song.csv
Fitted model parameters during song production for each bird.

Bird-Specific Data
Each bird has a corresponding folder inside the Data directory. Bird identifiers are listed in Table 2 of the Supplemental Material. Within each bird's folder:

data_<bird_name>.csv
Columns: time, expiratory EMG, total force, air sac pressure, and sound during song production.

data_quiet_<bird_name>.csv
Columns: time, force, and air sac pressure during quiet breathing.

magnification_quiet_<bird_name>.csv
Columns: frequencies and magnification values for quiet breathing parameters.

magnification_song_<bird_name>.csv
Columns: frequencies and magnification values for song production parameters.

quiet_rate_<bird_name>.csv
Column: measured breathing periods during quiet breathing.

song_rate_labels_<bird_name>.csv
Columns: breathing periods during song and corresponding labels based on pressure pattern clustering.

📁 Codes
All scripts are written in Python. Summary Files.

fig_data_model.py
Loads measured data during song, integrates the model, and generates Fig. 2 (or subfigures in Fig. S1).

fig_resonances.py
Loads processed results and generates Fig. 3 (or subfigures in Fig. S4).

fit_song.py
Fits model parameters to reproduce song data using a parallelized grid search.

fig_quiet_breathing.py
Loads quiet breathing data, integrates the model, and visualizes results.
