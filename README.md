#to create all relevant scripts and folders for a new protein, enter "python3.6 new_model.py", then enter the model name (i.e. "2bn3.pdb");
#if your model file contains a ligand, you can enter "python3.6 new_model.py -l" instead (not necessary, the program will also ask whether there is a ligand);

#once the model directory exists, you can move there and enter "sh bg_simulate_all &", then enter "disown", in order to simulate everything in the background;
#if you do this, you can close the console window while you wait; the simulations will likely take several hours (or days);
#bg_simulate_all will execute the following scripts:
#create_ensemble.sh: this will create 20 models and "shake" the coordinates for each one, then combine them into an "ensemble";
#create_reference.sh: this will create a reference file to ensure correct indexing by XDS;
#adjust_max.py: this will determine maximum values for each parameter based on CC1/2;
#simulate.py: this will simulate frames, then reduce and refine the data with an adjustable number of values (default: 20) for each parameter, 
#and save relevant data statistics in the "datafiles" directory;
#plot_all.py: this will create plots for the contents of the "datafiles" directory;

#you can check progress by looking at the logs in the "logfiles" directory;
#if you want to skip / change any of the steps, simply alter the relevant lines in bg_simulate_all, or run each script individually;

#if you want to simulate data for a specific parameter value, use "fix.py";
