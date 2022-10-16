# Simulation and evaluation of imperfect X-ray diffraction data

This repository contains everything needed to reproduce the results of my master's thesis:

**Linking data and model quality in X-ray Crystallography**

Successfully running the scripts in this repository will require at least 2 GB of free space, as well as much processing power and / or time.

For details on the purpose of this project, please look at *master_thesis_adler.pdf*

To reproduce my results, run *new_model.py* in the *models* directory.

Once the new directory exists, move there and enter "sh bg_simulate_all &", then enter "disown", in order to simulate everything in the background;
if you do this, you can close the console window while you wait; the simulations will likely take several hours (or days);
bg_simulate_all will execute the following scripts:
create_ensemble.sh: this will create 20 models and "shake" the coordinates for each one, then combine them into an "ensemble";
create_reference.sh: this will create a reference file to ensure correct indexing by XDS;
adjust_max.py: this will determine maximum values for each parameter based on CC1/2;
simulate.py: this will simulate frames, then reduce and refine the data with an adjustable number of values (default: 20) for each parameter, 
and save relevant data statistics in the "datafiles" directory;
plot_all.py: this will create plots for the contents of the "datafiles" directory;

You can check progress by looking at the logs in the "logfiles" directory.
If you want to skip / change any of the steps, simply alter the relevant lines in bg_simulate_all, or run each script individually.

To simulate data for a specific parameter value, use "fix.py";
