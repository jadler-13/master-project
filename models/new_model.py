#python3.6 new_model.py (-l)
import sys
import os

model = input("Please enter model name!\ninput example: 2bn3.pdb\n> ")
if not os.path.isfile(model):
    sys.exit("{} not found. Please make sure to move the model to this directory!".format(model))
name = model[:-4]
if os.path.isdir(name):
    sys.exit("Directory already exists. Please delete before proceeding.")
else:
    os.system("cp -r DO_NOT_EDIT {}".format(name))
os.system("ln -s {x}/new/SIMULATE.INP {x}/new/xds/XDS.INP".format(x = name))
try:
    ligand = sys.argv[1]
    if ligand != "-l":
        sys.exit("\ninvalid argument: {}".format(ligand))
except IndexError:
    query = input("Does your model file contain a ligand? (y/n)\n> ")
    if query[0] == "y":
        ligand = "-l"
    else:
        ligand = 0
newname = "{}_modified".format(name)            #default output for phenix.pdbtools is [input file name]_modified.pdb
newmodel = "{}.pdb".format(newname)
path = "{}/{}".format(os.getcwd(), name)
os.system("phenix.pdbtools {} set_b_iso=20 remove_alt_confs=True".format(model))
if ligand == "-l":
    os.system("phenix.ready_set {x}\nrm TLA* {x} {y}.eff\nmv {y}.ligands.cif {z}/{y}.ligands.cif".format(x = newmodel, y = newname, z = name))
    os.system("mv {}.updated.pdb {}".format(newname, newmodel))
os.system("grep -v ANISOU {x} > {y}/{x}".format(x = newmodel, y = name))
os.system("rm {}".format(newmodel))
os.system("phenix.pdbtools {} remove_alt_confs=True".format(model))
os.system("grep -v ANISOU {} > {}/{}_b_unchanged.pdb".format(newmodel, name, newname))

#prepare script for easy creation of ensemble model;
with open("{}/create_ensemble.sh".format(path), "r") as fin:
    with open("{}/temp_ce.sh".format(path), "w") as fout:
        for line in fin:
            if "phenix.refine" in line and ligand == "-l":
                fout.write(line.replace("refine.mtz", "refine.mtz ../{}.ligands.cif".format(newname)))
            elif "replace_this" in line:
                fout.write(line.replace("replace_this", newname))
            else:
                fout.write(line)
os.system("rm {x}/create_ensemble.sh\nmv {x}/temp_ce.sh {x}/create_ensemble.sh".format(x = path))

#determine unit cell constants;
with open(newmodel, "r") as fin:
    for line in fin:
        lc = line.split()
        if not lc:
            continue
        elif lc[0] == "CRYST1":
            cell_constants = lc
            cell_constants_str = ""
            for i in range(1, 7):
                cell_constants_str += "{} ".format(lc[i])
            break
        else:
            continue
print("\nunit cell constants: ",cell_constants_str)
os.system("rm {}".format(newmodel))

#set absolute path for image directory (important, as SIMULATE.INP will be called from 2 different directories);
with open("{}/new/SIMULATE.INP".format(path), "r") as fin:
    with open("{}/new/TEMP.INP".format(path), "w") as fout:
        for line in fin:
            lc = line.split()
            if not lc:
                fout.write(line)
                continue
            elif lc[0] == "!simulate" and lc[1] == "UNIT_CELL_A-AXIS=":
                newline = line.replace(lc[2], cell_constants[1])
                fout.write(newline)
            elif lc[0] == "!simulate" and lc[1] == "UNIT_CELL_B-AXIS=":
                newline = line.replace(lc[3], cell_constants[2])
                fout.write(newline)
            elif lc[0] == "!simulate" and lc[1] == "UNIT_CELL_C-AXIS=":
                newline = line.replace(lc[4], cell_constants[3])
                fout.write(newline)
            elif lc[0] == "UNIT_CELL_CONSTANTS=":
                fout.write("{} {}\n".format(lc[0], cell_constants_str))
            elif lc[0] == "NAME_TEMPLATE_OF_DATA_FRAMES=":
                newline = line.replace(lc[1], "{}/new/xds/image/insu_???.cbf".format(path))
                fout.write(newline)
            else:
                fout.write(line)
os.system("rm {x}/new/SIMULATE.INP\nmv {x}/new/TEMP.INP {x}/new/SIMULATE.INP".format(x = path))

#create SIMULATE.INP for ref directory;
with open("{}/new/SIMULATE.INP".format(path), "r") as fin:
    with open("{}/ref/SIMULATE.INP".format(path), "w") as fout:
        for line in fin:
            lc = line.split()
            if not lc:
                fout.write(line)
                continue
            elif lc[0] == "NAME_TEMPLATE_OF_DATA_FRAMES=":
                newline = line.replace(lc[1], "{}/ref/xds/image/insu_???.cbf".format(path))
                fout.write(newline)
            else:
                fout.write(line)
os.system("ln -s {x}/ref/SIMULATE.INP {x}/ref/xds/XDS.INP".format(x = path))

#update run.sh to use current model;
with open("{}/new/run.sh".format(path), "r") as fin:
    with open("{}/new/temp_run.sh".format(path), "w") as fout:
        for line in fin:
            if "phenix.refine" in line and ligand == "-l":
                templine = line.replace("reflections_1.mtz", "reflections_1.mtz ../../{}.ligands.cif".format(newname))
                fout.write(templine.replace("replace_this", newname))
            elif "replace_this" in line:
                fout.write(line.replace("replace_this", newname))
            else:
                fout.write(line)
os.system("rm {x}/new/run.sh\nmv {x}/new/temp_run.sh {x}/new/run.sh".format(x = path))

#create run_sim.sh and run_xds.sh;
with open("{}/new/run.sh".format(path), "r") as fin:
    with open("{}/new/run_sim.sh".format(path), "w") as fout:
        for line in fin:
            if "#break1" in line:
                break
            else:
                fout.write(line)
with open("{}/new/run.sh".format(path), "r") as fin:
    with open("{}/new/run_xds.sh".format(path), "w") as fout:
        for line in fin:
            if "#break2" in line:
                break
            else:
                fout.write(line)

#update run1.sh and run2.sh in ref directory to use current model;
with open("{}/ref/run1.sh".format(path), "r") as fin:
    with open("{}/ref/temp_run.sh".format(path), "w") as fout:
        for line in fin:
            if "replace_this" in line:
                fout.write(line.replace("replace_this", newname))
            else:
                fout.write(line)
os.system("rm {x}/ref/run1.sh\nmv {x}/ref/temp_run.sh {x}/ref/run1.sh".format(x = path))
with open("{}/ref/run2.sh".format(path), "r") as fin:
    with open("{}/ref/temp_run.sh".format(path), "w") as fout:
        for line in fin:
            if "phenix.refine" in line and ligand == "-l":
                templine = line.replace("reflections_1.mtz", "reflections_1.mtz ../../{}.ligands.cif".format(newname))
                fout.write(templine.replace("replace_this", newname))
            elif "replace_this" in line:
                fout.write(line.replace("replace_this", newname))
            else:
                fout.write(line)
os.system("rm {x}/ref/run2.sh\nmv {x}/ref/temp_run.sh {x}/ref/run2.sh".format(x = path))

#prepare script for easy creation of reference model;
with open("{}/create_reference.sh".format(path), "r") as fin:
    with open("{}/temp_cr.sh".format(path), "w") as fout:
        for line in fin:
            if "replace_this" in line:
                fout.write(line.replace("replace_this", newname))
            else:
                fout.write(line)
os.system("rm {x}/create_reference.sh\nmv {x}/temp_cr.sh {x}/create_reference.sh".format(x = path))

#update plot_all.py to use current model;
with open("{}/plot_all.py".format(path), "r") as fin:
    with open("{}/temp_plot.py".format(path), "w") as fout:
        for line in fin:
            if "replace_this" in line:
                fout.write(line.replace("replace_this", newname))
            else:
                fout.write(line)
os.system("mv {x}/temp_plot.py {x}/plot_all.py".format(x = path))

print("\n{x} directory created. Next, create the ensemble and then the reference model, using the respective scripts in {x}.".format(x = name))
