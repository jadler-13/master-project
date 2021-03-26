#python3.6 simulate.py
import os
import math
import numpy as np
pnumber = int(input("please enter parameter number:\n0: wavelength\n1: beam\n2: cell\n3: orientation\n>"))
newvalue = input("please enter parameter value:\n")
parameter = ["WAVELENGTH_STDDEV","BEAM_STDDEV","CELL_STDDEV", "ORIENTATION_STDDEV"]
os.system("rm -r FIX")
os.system("cp -r new FIX")
pathin = "new/SIMULATE.INP"
pathout = "FIX/SIMULATE.INP"


#read starting values for chosen parameter from input file and store them in a list;
with open(pathin, 'r') as inp:
    for line in inp:
        if parameter[pnumber] in line:
            print("line contents: {}".format(line.split())) #for debugging;
            for element in line.split():
                try:
                    value = float(element)
                    break
                except ValueError:
                    pass
value = str(value)

for m in range(1):
    with open(pathin, 'r') as fin:
        with open(pathout, 'w') as fout:
            for line in fin:
                if parameter[pnumber] in line:
                    newline = line.replace("!(don't)", "!")
                    fout.write(newline.replace(value, newvalue))
                elif "NAME_TEMPLATE_OF_DATA_FRAMES=" in line:
                    fout.write(line.replace("new", parameter[pnumber]))
                else:
                    fout.write(line)
    os.system("rm temp.INP")
    os.system("cp {} temp.INP".format(pathout))
    os.chdir(parameter[pnumber])

    #run sim_mx for different resolution ranges and record # of overlapping + contributing pixels;
    ranges = ["50 30", "29.9 10", "9.9 5", "4.9 3", "2.9 1.79"]
    for i in range(5):
        #print("\ni = {}\n".format(i))           #DEBUG
        with open ("../temp.INP", "r") as fin:
            with open ("SIMULATE.INP", "w") as fout:
                for line in fin:
                    if "INCLUDE_RESOLUTION_RANGE=" in line:
                        fout.write("INCLUDE_RESOLUTION_RANGE= {}\n".format(ranges[i]))
                    else:
                        fout.write(line)
        os.system("sh run_sim.sh")
        with open ("simlog", "r") as fin:
            with open ("../datafiles/fix_reso_pixels{}".format(str(pnumber)), "a") as fout:
                if i == 0:
                    fout.write("{} = {}:\n".format(parameter[pnumber], newvalue))
                fout.write("resolution range {}".format(ranges[i]))
                simcount = 0
                for line in fin:
                    lc = line.split()
                    if simcount == 2:
                        break
                    if lc[0] == "average":
                        fout.write(line)
                        simcount += 1

    os.system("rm SIMULATE.INP\ncp ../temp.INP SIMULATE.INP")
    os.system("sh run.sh") #simulate frames with new parameter values and analyse them;
    #os.system("cp xds/FRAME.cbf ../frames/{pn}_{val}_FRAME.cbf".format(pn = str(pnumber), val = str(values[0])[:6]))
    
    #write starting and final R-factors for phenix.refine to output file 1;
    with open ("refine/rs", "r") as fin:
        with open ("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            fout.write("{} = {}:\n".format(parameter[pnumber], newvalue))
            for line in fin:
                fout.write(line)
    
    #write beam divergence e.s.d and reflecting range e.s.d determined by xds to output file 1;
    with open ("xds/esd", "r") as fin:
        with open ("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            for line in fin:
                fout.write(line)

    #write cc1/2 for highest resolution bin to output file 1;
    with open ("xds/cc12", "r") as fin:
        with open ("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            for line in fin:
                lc = line.split()
                fout.write("CC(1/2) in resolution bin {x} is: {y}\n".format(x = lc[0], y = lc[10]))

    #write r.m.s.d. to output file 1;
    with open ("refine/rmsd", "r") as fin:
        with open ("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            for line in fin:
                if "r.m.s.d:" in line:
                    fout.write(line)
                    break
    
    #write ISa (determined by XDS) to output file 1;
    ISacount = 0
    with open ("xds/xdslog", "r") as fin:
        with open ("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            for line in fin:
                if ISacount == 2:
                    objects = line.split()
                    fout.write("ISa= {}\n".format(objects[2]))
                    ISacount = 0
                if "ISa" in line:
                    ISacount += 1

    #write protein mean B-factor and resolution-dependent R-factors determined by phenix.refine to output file 1;
    rescount = 0
    bcount = 0
    with open ("refine/reflog", "r") as fin:
        with open ("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            for line in fin:
                if "Resolution" in line:
                    #print(line) #for debugging;
                    lc = line.split()
                    if lc[0]=="|" and lc[3]=="Compl.":
                        rescount += 1
                if rescount==2:
                    lc = line.split()
                    if lc[1]=="10:":
                        fout.write(line)
                        break
                    else:
                        fout.write(line)
                if "Individual ADP refinement" in line:
                    bcount += 1
                if bcount==6:
                    lc = line.split()
                    if lc!=[] and lc[0]=="Protein:":
                        fout.write("Protein mean B-factor = {}\n".format(str(lc[3])))
                        continue
    
    #write Wilson B-factors for F-obs and F-model to output file 1;
    with open("refine/wilsonlog1", "r") as fin:
        with open("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            for line in fin:
                if "Least squares" in line:
                    lc = line.split()
                    fout.write("Wilson B-factor for F-obs: {}".format(str(lc[7][:4])))
                    break
    with open("refine/wilsonlog2", "r") as fin:
        with open("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            for line in fin:
                if "Least squares" in line:
                    lc = line.split()
                    fout.write("\nWilson B-factor for F-model: {}\n".format(str(lc[7][:4])))
                    break
    
    #write "observed" (simulated) and model amplitudes, as well as Resolution and S**2 values to output file 2;
    with open("wilsonrefs", "r") as fin:
        with open("../datafiles/fix_wilsonvalues{}".format(str(pnumber)), "a") as fout:
            fout.write("\n {} = {}:\nReso   S**2   F-obs   F-model".format(parameter[pnumber], newvalue))
            lowcount = 0
            reslimit = 3
            start = False
            for line in fin:
                lc = line.split()
                if len(lc) < 6:
                    continue
                if lc[0] == "User:":
                    start = True
                    continue
                if not start:
                    continue
                if lc[0] == "List":
                    break
                if lc[4] == "?" or lc[5] == "?":
                    continue
                Ssq = float(lc[3])
                reso = 1/math.sqrt(Ssq)
                if reso > reslimit:
                    lowcount += 1
                    continue
                fout.write("\n{x1}   {x2}   {x3}   {x4}".format(x1 = str(reso)[:4], x2 = str(Ssq)[:6], x3 = str(lc[4]), x4 = str(lc[5]))) 
            print("{x1} reflections had resolution values >{x2}".format(x1 = lowcount, x2 = reslimit))

    #write geometry quality indicators to output file 3;
    with open("refine/molprobity.out", "r") as fin:
        with open ("../datafiles/fix_geometry{}".format(str(pnumber)), "a") as fout:
            fout.write("{} = {}:\n".format(parameter[pnumber], newvalue))
            for line in fin:
                lc = line.split()
                if not lc:
                    continue
                elif lc[0] == "whole:":
                    fout.write("Rama-Z "+line)
                elif lc[0] == "Clashscore" or lc[0] == "RMS(bonds)" or lc[0] == "RMS(angles)":
                    fout.write(line)
                elif lc[0] == "MolProbity":
                    fout.write(line)
                    break

    #write number of indexed reflections to output file 1;
    indexcount = 0
    with open("xds/IDXREF.LP", "r") as fin:
        with open("../datafiles/fix_rvalues{}".format(str(pnumber)), "a") as fout:
            for line in fin:
                if indexcount == 1:
                    fout.write(line)
                    break
                elif "INDEXING OF OBSERVED" in line:
                    fout.write(line)
                    indexcount += 1

