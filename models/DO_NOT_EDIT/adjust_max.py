#python3.6 simulate.py
import os
import numpy as np
parameter = ["WAVELENGTH_STDDEV","BEAM_STDDEV","CELL_STDDEV", "ORIENTATION_STDDEV"]
values = []
old = []
new = []
check = 10

def get_diff():
    with open("xds/cc12", "r") as fin:
        for line in fin:
            lc = line.split()
            print("\nCC(1/2) in resolution bin {x} is: {y}\n".format(x = lc[0], y = lc[10]))
            return float(lc[10][:2]) - 40

for i in range(len(parameter)):
    directory = "{}_temp".format(parameter[i])
    os.system("rm -r {x}\ncp -r new {x}".format(x = directory))
    #activate current parameter and store starting values in a list;
    values.append([])
    with open("new/SIMULATE.INP", "r") as fin:
        with open("{}/SIMULATE.INP".format(directory), "w") as fout:
            for line in fin:
                if parameter[i] in line:
                    #print("line contents: {}".format(line.split())) #for debugging;
                    for value in line.split():
                        try:
                            values[i].append(float(value))
                        except ValueError:
                            pass
                    fout.write(line.replace("(don't)", ""))
                elif "NAME_TEMPLATE_OF_DATA_FRAMES=" in line:
                    fout.write(line.replace("new", directory))
                else:
                    fout.write(line)
    values[i] = np.array(values[i])
    oldval = list(map(str, values[i]))
    oldstr = ""
    for j in range(len(values[i])):
        oldstr = oldstr + oldval[j][:6] + " "
    old.append(oldstr)
    os.chdir(directory)
    os.system("sh run_xds.sh")
    diff = get_diff()
    while not (-3) <= diff <= 4:
        if diff == check != 100:
            inc = (-0.001)
        else:
            if diff > 35:
                mod = 0.01
            elif diff > 30:
                mod = 0.007
            elif diff > 25:
                mod = 0.005
            elif diff > 15:
                mod = 0.002
            else:
                mod = 0.001
            inc = diff * mod * values[i][0]
        newvalues = values[i] + inc
        val = list(map(str, values[i])) #convert values to strings to allow replacement later;
        newval = list(map(str, newvalues))
        valstr = ""
        newstr = ""
        for n in range(len(values[i])):
            valstr = valstr + val[n][:6] + " "
            newstr = newstr + newval[n][:6] + " "
        print("\nvalues: {}\nnew values: {}".format(valstr, newstr))
        with open("SIMULATE.INP", "r") as fin:
            with open("temp.INP", "w") as fout:
                for line in fin:
                    if parameter[i] in line:
                        fout.write(line.replace(valstr, newstr))
                    else:
                        fout.write(line)
        values[i] = newvalues
        os.system("mv temp.INP SIMULATE.INP\nsh run_xds.sh")
        diff = get_diff()
    new.append(newstr)
    os.chdir("..")

with open("{}/SIMULATE.INP".format(directory), "r") as fin:
    with open("new/SIMULATE.INP", "w") as fout:
        for line in fin:
            if directory in line:
                fout.write(line.replace(directory, "new"))
                continue
            else:
                stop = 0
                for i in range(len(parameter)):
                    if parameter[i] in line:
                        fout.write(line.replace(old[i], new[i]))
                        stop = 1
                if stop == 1:
                    continue
                else:
                    fout.write(line)

for i in range(len(parameter)):
    os.system("rm -r {}_temp".format(parameter[i]))

