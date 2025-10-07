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

See at: https://doi.org/10.5281/zenodo.15636149

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
