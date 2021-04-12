#python3.6 plot_all.py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
matplotlib.use('Agg')

parameters = ["WAVELENGTH_STDDEV", "BEAM_STDDEV", "CELL_STDDEV", "ORIENTATION_STDDEV"]

#this function unites the legends for two separate axes plotted in the same figure;
def unite_legends(axes):
    h, l = [], []
    for ax in axes:
        tmp = ax.get_legend_handles_labels()
        h.extend(tmp[0])
        l.extend(tmp[1])
    return h, l

#these functions calculate mean for all F_obs or F_model values in one resolution shell and add mean to respective list;
def add_to_FP(x, y, z):
    FPi = np.array(x)
    mn = np.mean(FPi, dtype=np.float64)
    FP[z][y].append(float(str(mn)[:4]))
def add_to_FM(x, y, z):
    FMi = np.array(x)
    mn = np.mean(FMi, dtype=np.float64)
    FM[z][y].append(float(str(mn)[:4]))

#define empty lists for all values to be plotted;
parameter = []
rmsd = []
parval = []
resrange = []
B_origin = []
B_origin_het = []
B_ref = []
B_ref_het = []
overlap = []
contrib = []
cc12 = []
beam_esd = []
ref_esd = []
startrwork = []
startrfree = []
finalrwork = []
finalrfree = []
ISa = []
B = []
Wilson_B_obs = []
Wilson_B_model = []
ramaz = []
clashscore = []
rms_bonds = []
rms_angles = []
mpscore = []
rwork = []
rfree = []
parline = []
FP = []
FM = []
FPi = []
FMi = []
reso = []
ssq = []

print("reading values...")

#read atomic B-factors from original model and add them to respective lists;
with open("replace_this_b_unchanged.pdb", "r") as fin:
    for line in fin:
        col = line.split()
        if not col:
            continue
        elif col[0] == "ATOM":
            if len(col) < 12:
                if len(col[2]) < 5 and len(col[9]) > 5:
                    B_origin.append(float(col[9][4:]))
                else:
                    B_origin.append(float(col[9]))
            else:
                B_origin.append(float(col[10]))
        elif col[0] == "HETATM":
            if len(col) < 11:
                if len(col[8]) > 4:
                    B_origin_het.append(float(col[8][4:]))
                else:
                    print("new error when reading atomic B-factor from line:\n",col)
            elif len(col[9]) > 5:
                B_origin_het.append(float(col[9][4:]))
            elif len(col) < 12:
                B_origin_het.append(float(col[9]))
            else:
                B_origin_het.append(float(col[10]))
        elif col[0] == "HETATM":
            if len(col) < 11:
                if len(col[8]) > 4:
                    B_origin_het.append(float(col[8][4:]))
                else:
                    print("new error when reading atomic B-factor from line:\n",col)
            elif len(col[9]) > 5:
                B_origin_het.append(float(col[9][4:]))
            elif len(col) < 12:
                B_origin_het.append(float(col[9]))
            else:
                B_origin_het.append(float(col[10]))

#read atomic B-factors from refined simulation models and add them to respective lists;
for i in range(4):
    with open("{}/refine/replace_this_refine_001.pdb".format(parameters[i]), "r") as fin:
        B_ref.append([])
        B_ref_het.append([])
        for line in fin:
            col = line.split()
            if not col:
                continue
            elif col[0] == "ATOM":
                if len(col) < 12:
                    if len(col[2]) < 5 and len(col[9]) > 5:
                        B_ref[i].append(float(col[9][4:]))
                    else:
                        B_ref[i].append(float(col[9]))
                else:
                    B_ref[i].append(float(col[10]))
            elif col[0] == "HETATM":
                if len(col) < 11:
                    if len(col[8]) > 4:
                        B_ref_het[i].append(float(col[8][4:]))
                    else:
                        print("new error when reading atomic B-factor from line:\n",col)
                elif len(col[9]) > 5:
                    B_ref_het[i].append(float(col[9][4:]))
                elif len(col) < 12:
                    B_ref_het[i].append(float(col[9]))
                else:
                    B_ref_het[i].append(float(col[10]))

#read overlapping and contributing pixels per reflection in different resolution ranges
for i in range(4):
    cyclecount = 0
    overlap.append([])
    contrib.append([])
    with open("datafiles/reso_pixels{}".format(i), 'r') as fin:
        for line in fin:
            lc = line.split()   #lc = line content
            if lc[1] == "=" and cyclecount == 0:
                cyclecount = 1
                rangecount = 0
                allranges = 0
            elif lc[1] == "=":
                rangecount = 0
                allranges = 1
            elif "resolution range" in line:
                if not allranges:
                    resrange.append("{} - {} Å".format(lc[2], lc[3]))
                    overlap[i].append([])
                    contrib[i].append([])
                overlap[i][rangecount].append(float(lc[-1:][0]))
            elif "contributing" in line:
                contrib[i][rangecount].append(float(lc[-1:][0]))
                rangecount +=1


#read various indicators from datafiles/rvalues{0-3} and add them to respective lists;
for i in range(4):
    parameter.append(0)
    parval.append([])
    startrwork.append([])
    startrfree.append([])
    finalrwork.append([])
    finalrfree.append([])
    ISa.append([])
    B.append([])
    Wilson_B_obs.append([])
    Wilson_B_model.append([])
    rmsd.append([])
    beam_esd.append([])
    ref_esd.append([])
    cc12.append([])
    rwork.append([])
    rfree.append([])
    parline.append([])
    parcount = 0
    bincount = 0
    with open("datafiles/rvalues{}".format(i), 'r') as fin:
        for line in fin:
            lc = line.split()   #lc = line content
            if lc[1] == "=" and parameter[i] == 0:
                parameter[i] = lc[0]
                parval[i].append(float(lc[2]))
                parline[i].append(lc)
                rwork[i].append([])
                rfree[i].append([])
                cyclecount = 1
            elif lc[1] == "=":
                parval[i].append(float(lc[2]))
                parline[i].append(lc)
                rwork[i].append([])
                rfree[i].append([])
                #print("i = {}\nparameter[i] =  {}\nline = {}".format(i, parameter[i], line))
                cyclecount += 1
            elif lc[0] == "|" and lc[1] == "1:":
                #print("parameter: {}\ncycle: {}\n".format(parameter[i], cyclecount))
                rwork[i][parcount].append(float(lc[8]))
                rfree[i][parcount].append(float(lc[9]))
                bincount = 1
            elif lc[0] == "|" and 0 < bincount < 9:
                #print(line)
                rwork[i][parcount].append(float(lc[8]))
                rfree[i][parcount].append(float(lc[9]))
                bincount += 1
            elif lc[0] == "|" and bincount == 9:
                #print(line)
                rwork[i][parcount].append(float(lc[8]))
                rfree[i][parcount].append(float(lc[9]))
                bincount = 0
                parcount += 1
            elif lc[0] == "Start":
                startrwork[i].append(float(lc[3].strip(",")))
                startrfree[i].append(float(lc[6]))
            elif lc[0] == "Final":
                finalrwork[i].append(float(lc[3].strip(",")))
                finalrfree[i].append(float(lc[6]))
            elif lc[0] == "ISa=":
                ISa[i].append(float(lc[1]))
            elif lc[0] == "Protein":
                B[i].append(float(lc[4]))
            elif lc[0] == "Wilson" and lc[3] == "F-obs:":
                Wilson_B_obs[i].append(float(lc[4]))
            elif lc[0] == "Wilson" and lc[3] == "F-model:":
                Wilson_B_model[i].append(float(lc[4]))
            elif lc[0] == "r.m.s.d:":
                rmsd[i].append(float(lc[1]))
            elif lc[0] == "BEAM_DIVERGENCE=":
                beam_esd[i].append(float(lc[3]))
            elif lc[0] == "REFLECTING_RANGE=":
                ref_esd[i].append(float(lc[3]))
            elif lc[0] == "CC(1/2)":
                cc12[i].append(float(lc[6].strip("*")))
            elif lc[0] == "average":
                if lc[3] == "overlapping":
                    overlap[i].append(float(lc[8]))
                elif lc[3] == "contributing":
                    contrib[i].append(float(lc[8]))
            
#read F_obs and F_model intensities from datafiles/wilsonvalues{0-3} and add them to lists; also add resolution shell and s**2 values;
for i in range(4):
    FP.append([])
    FM.append([])
    FPi.append([])
    FMi.append([])
    reso.append([])
    ssq.append([])
    wilsoncount = -1
    with open("datafiles/wilsonvalues{}".format(i), 'r') as fin:
        for line in fin:
            lc  = line.split()
            if not lc:
                continue
            if lc[1] == "=":
                if FPi[i]:
                    add_to_FP(FPi[i], wilsoncount, i)
                    add_to_FM(FMi[i], wilsoncount, i)
                    FPi[i] = []
                    FMi[i] = []
                wilsoncount += 1
                FP[i].append([])
                FM[i].append([])
                reso[i].append([])
                ssq[i].append([])
                new = True
            elif lc[0] == "Reso":
                continue
            else:
                if reso[i][wilsoncount] and reso[i][wilsoncount][-1] != (math.floor(float(lc[0])*10+0.5)/10):
                    new = True
                if new:
                    if FPi[i]:
                        add_to_FP(FPi[i], wilsoncount, i)
                        add_to_FM(FMi[i], wilsoncount, i)
                    reso[i][wilsoncount].append(math.floor(float(lc[0])*10+0.5)/10)
                    ssq[i][wilsoncount].append(float(lc[1]))
                    FPi[i] = [float(lc[2])]
                    FMi[i] = [float(lc[3])]
                    new = False
                elif not new:
                    FPi[i].append(float(lc[2]))
                    FMi[i].append(float(lc[3]))
        else:
            add_to_FP(FPi[i], wilsoncount, i)
            add_to_FM(FMi[i], wilsoncount, i)

#add NaN values to fill gaps in resolution bins;
reso_bins = len(reso[0][0])
for i in range(4):
    for j in range(len(parval[0])):
        try:
            indextest = reso[i][j][reso_bins-1]
        except IndexError:
            reso[i][j].append(1.8)
            ssq[i][j].append(0.2922)
            FP[i][j].append(np.nan)
            FM[i][j].append(np.nan)

#check array shape;
reso_bins = len(reso[0][0])
if len(FP) != 4 or len(FM) != 4:
    sys.exit("unexpected number of parameters! (expected: 4)")
if not len(parval[0]) == len(parval[1]) == len(parval[2]) == len(parval[3]):
    sys.exit("number of simulated values differs between parameters!")
for i in range(4):
    if len(FP[i]) != len(parval[i]) or len(FM[i]) != len(parval[i]):
        sys.exit("number of intensity arrays does not match number of simulated values!")
    for j in range(20):
        if len(reso[i][j]) != reso_bins:
            sys.exit("number of resolution bins differs between parameters! (check whether resolutions of <1.85 were recorded at max STDDEV-values)")

'''
#FOR DEBUGGING:
for i in range(4):
    print(parameters[i],":\n")
    for j in range(20):
        print("reso[{}][{}]: ".format(i, j), reso[i][j])
        print("ssq[{}][{}]: ".format(i, j), ssq[i][j])
        print("FP[{}][{}]: ".format(i, j), FP[i][j])
        print("FM[{}][{}]: ".format(i, j), FM[i][j],"\n")
'''

#convert amplitudes F into natural logarithm of intensities I;
FP = np.log(np.square(FP))
FM = np.log(np.square(FM))

#read various geometry indicators from datafiles/geometry{0-3} and add them to respective lists;
for i in range(4):
    ramaz.append([])
    clashscore.append([])
    rms_bonds.append([])
    rms_angles.append([])
    mpscore.append([])
    with open("datafiles/geometry{}".format(i), 'r') as fin:
        for line in fin:
            lc = line.split()
            if lc[0] == "Rama-Z":
                ramaz[i].append(float(lc[2]))
            elif lc[0] == "Clashscore":
                clashscore[i].append(float(lc[2]))
            elif lc[0] == "RMS(bonds)":
                rms_bonds[i].append(float(lc[2]))
            elif lc[0] == "RMS(angles)":
                rms_angles[i].append(float(lc[2]))
            elif lc[0] == "MolProbity":
                mpscore[i].append(float(lc[3]))

valnum = len(parval[0])

#create plots for number of overlapping / contributing pixels per reflection;
print("plotting overlapping / contributing pixels per reflection for different resolution ranges...")
res_added = np.add(overlap, contrib)
for i in range(4):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(14, 6))
    x = parval[i]
    mark = ["o", "s", "v", "^", "<", ">"]
    for j in range(rangecount):
        ax1.plot(x, np.array(overlap[i][j]), label=resrange[j], marker=mark[j])
        ax2.plot(x, np.array(contrib[i][j]), label=resrange[j], marker=mark[j])
    ax1.set_xlabel(parameter[i])
    ax1.set_ylabel("average # of overlapping pixels per reflection")
    ax2.set_xlabel(parameter[i])
    ax2.set_ylabel("average # of contributing pixels per reflection")
    ax1.set_title("Reflection overlap as a function of\n{} in different resolution ranges".format(parameter[i]), fontweight="bold")
    ax2.set_title("Number of contributing pixels per reflection as a\nfunction of {} in different resolution ranges".format(parameter[i]), fontweight="bold")
    ax1.legend()
    ax1.set_ylim(0, None)
    ax2.set_ylim(0, None)
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(11)
    plt.savefig("plots/reso_pixels{}.png".format(i), bbox_inches="tight")
    plt.close()
for i in [0, 2]:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    mark = ["o", "s", "v", "^", "<", ">"]
    for j in range(rangecount):
        ax1.plot(parval[i], np.array(res_added[i][j]), label=resrange[j], marker=mark[j])
        ax2.plot(parval[i+1], np.array(res_added[i+1][j]), label=resrange[j], marker=mark[j])
    ax1.set_xlabel(parameter[i])
    ax2.set_xlabel(parameter[i+1])
    ax1.set_ylabel("overlapping + contributing pixels per reflection")
    ax1.set_title("Overlapping + contributing pixels as a function\nof {} in different resolution ranges".format(parameter[i]), fontweight="bold")
    ax2.set_title("Overlapping + contributing pixels as a function\nof {} in different resolution ranges".format(parameter[i+1]), fontweight="bold")
    ax1.legend()
    ax1.set_ylim(0, None)
    ax2.set_ylim(0, None)
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(11)
    plt.savefig("plots/sum_reso_pixels_{}_{}.png".format(i, i+1), bbox_inches="tight")
    plt.close()

#create plots for R-factors, mean protein B-factors, Wilson B-factors, and ISa;
print("plotting R-factors, B-factors, and ISa...")
for i in range(4):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(14, 6))
    ax1.plot(parval[i], np.array(startrwork[i]), label=r'start-R$_{work}$', marker="o")
    ax1.plot(parval[i], np.array(startrfree[i]), label=r'start-R$_{free}$', marker="s")
    ax1.plot(parval[i], np.array(finalrwork[i]), label=r'final-R$_{work}$', marker="^")
    ax1.plot(parval[i], np.array(finalrfree[i]), label=r'final-R$_{free}$', marker="v")
    ax2.plot(parval[i], np.array(ISa[i]), label='ISa', color = 'red', marker="o")
    ax1.set_xlabel(parameter[i])
    ax2.set_xlabel(parameter[i])
    ax3 = ax2.twinx()
    ax3.plot(parval[i], np.array(B[i]), label='Protein mean B-factor', marker="s")
    ax3.plot(parval[i], np.array(Wilson_B_obs[i]), label=r'Wilson B-factor F$_{sim}$', marker="^")
    ax3.plot(parval[i], np.array(Wilson_B_model[i]), label=r'Wilson B-factor F$_{model}$', marker="v")
    ax1.set_ylabel('R-factors')
    ax2.set_ylabel('ISa')
    ax3.set_ylabel('B-factors')
    ax1.set_title('R-factors as functions of {}'.format(parameter[i]), fontweight="bold")
    ax2.set_title('ISa and B-factors as functions of {}'.format(parameter[i]), fontweight="bold")
    ax1.legend()
    handles, labels = unite_legends([ax2, ax3])
    ax3.legend(handles, labels, loc="upper center")
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax3.xaxis.label, ax3.yaxis.label] + ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(11)
    plt.savefig("plots/main{}.png".format(i), bbox_inches="tight")
    plt.close()

#create plots with shared x and y axes for resolution dependent r-values;
print("plotting resolution dependent R-factors...")
#DEBUG
#for i in range(4):
#    for j in range(len(parline[i])):
#        print("len(rwork[{}][{}] = {}\n".format(i, j, len(rwork[i][j])))
#        print("len(rfree[{}][{}] = {}\n".format(i, j, len(rfree[i][j])))
rworkmin = np.amin(rwork)
rfreemax = np.amax(rfree)
x = np.linspace(1, 10, 10)
for i in range(4):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(14, 6))
    for j in range(len(parline[i])):
        ix = float(j) / valnum
        if j == 0:
            ax1.plot(x, np.array(rwork[i][j]), label=parline[i][j][2], color=(ix, 0.2, 0.2), marker='o')
            ax2.plot(x, np.array(rfree[i][j]), label=parline[i][j][2], color=(ix, 0.2, 0.2), marker='o')
        elif j == len(parline[i])-1:
            ax1.plot(x, np.array(rwork[i][j]), label=parline[i][j][2], color=(ix, 0.2, 0.2), marker='s')
            ax2.plot(x, np.array(rfree[i][j]), label=parline[i][j][2], color=(ix, 0.2, 0.2), marker='s')
        else:
            ax1.plot(x, np.array(rwork[i][j]), label=parline[i][j][2], color=(ix, 0.2, 0.2))
            ax2.plot(x, np.array(rfree[i][j]), label=parline[i][j][2], color=(ix, 0.2, 0.2))
    ax1.set_xlabel('Resolution Bin Number')
    ax2.set_xlabel('Resolution Bin Number')
    ax1.set_ylabel(r'R$_{work}$')
    ax2.set_ylabel(r'R$_{free}$')
    ax1.set_title(r'R$_{work}$ by resolution and '+parline[i][0][0], fontweight="bold")
    ax2.set_title(r'R$_{free}$ by resolution and '+parline[i][0][0], fontweight="bold")
    ax1.set_ylim([rworkmin, rfreemax])
    ax2.set_ylim([rworkmin, rfreemax])
    ax1.legend(ncol=2)
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(11)
    plt.savefig("plots/res{}.png".format(i), bbox_inches="tight")
    plt.close()

#optional: display 2 wilson plots in one figure;
print("creating Wilson plots...")
reso = np.array(reso)
ssq = np.array(ssq)
for i in [0, 2]:
    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(14, 6))
    ax2 = ax1.twinx()
    ax4 = ax3.twinx()
    for j in range(wilsoncount+1):
        ix = float(j) / valnum
        ax1.plot(ssq[i][j], FP[i][j], label='F-obs '+str(parline[i][j][2]), color=(ix, 0.1, 0.5), marker='o')
        ax2.plot(ssq[i][j], FM[i][j], label='F-model '+str(parline[i][j][2]), color=(ix, 0.5, 0.1), marker='s')
        ax3.plot(ssq[i+1][j], FP[i+1][j], label='F-obs '+str(parline[i+1][j][2]), color=(ix, 0.1, 0.5), marker='o')
        ax4.plot(ssq[i+1][j], FM[i+1][j], label='F-model '+str(parline[i+1][j][2]), color=(ix, 0.5, 0.1), marker='s')
    ax1.set_xlabel('Resolution [Å]')
    ax1.set_ylabel(r'$ln(mn(F_{sim}^2))$', color='purple')
    ax1.set_title('Wilson plot for '+parline[i][0][0], fontweight="bold")
    ax3.set_xlabel('Resolution [Å]')
    ax4.set_ylabel(r'$ln(mn(F_{model}^2))$', color='green')
    ax3.set_title('Wilson plot for '+parline[i+1][0][0], fontweight="bold")
    plt.sca(ax1)
    plt.xticks(ssq[i][0], reso[i][0])
    plt.sca(ax3)
    plt.xticks(ssq[i+1][0], reso[i+1][0])
    every_nth = 2
    for n, label in enumerate(ax1.xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)
    for n, label in enumerate(ax3.xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)
    ax1_ymin, ax1_ymax = ax1.get_ylim()
    ax2_ymin, ax2_ymax = ax2.get_ylim()
    if ax1_ymin < ax2_ymin:
        y_min1 = ax1_ymin
    else:
        y_min1 = ax2_ymin
    if ax1_ymax > ax2_ymax:
        y_max1 = ax1_ymax
    else:
        y_max1 = ax2_ymax
    ax3_ymin, ax3_ymax = ax3.get_ylim()
    ax4_ymin, ax4_ymax = ax4.get_ylim()
    if ax3_ymin < ax4_ymin:
        y_min2 = ax3_ymin
    else:
        y_min2 = ax4_ymin
    if ax3_ymax > ax4_ymax:
        y_max2 = ax3_ymax
    else:
        y_max2 = ax4_ymax
    ax1.set_ylim(y_min1, y_max1)
    ax2.set_ylim(y_min1, y_max1)
    ax3.set_ylim(y_min2, y_max2)
    ax4.set_ylim(y_min2, y_max2)
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] + ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax4.xaxis.label, ax4.yaxis.label] + ax4.get_xticklabels() + ax4.get_yticklabels()):
        item.set_fontsize(11)
    plt.savefig("plots/wilson_{}_{}.png".format(i, i+1), bbox_inches="tight")
    plt.close()
'''
#create Wilson plot; x-axis is s**2, but corresponding resolution shell is used for x-labels;
print("creating Wilson plots...")
reso = np.array(reso)
ssq = np.array(ssq)
for i in range(4):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax2 = ax.twinx()
    for j in range(wilsoncount+1):
        ix = float(j) / valnum
        ax.plot(ssq[i][j], FP[i][j], label='F-obs '+str(parline[i][j][2]), color=(ix, 0.1, 0.5), marker='o')
        ax2.plot(ssq[i][j], FM[i][j], label='F-model '+str(parline[i][j][2]), color=(ix, 0.5, 0.1), marker='s')
    ax.set_xlabel('Resolution [Å]')
    ax.set_ylabel(r'$ln(mn(F_{sim}^2))$', color='purple')
    ax2.set_ylabel(r'$ln(mn(F_{model}^2))$', color='green')
    ax.set_title('Wilson plot for '+parline[i][0][0])
    plt.xticks(ssq[i][0], reso[i][0])
    every_nth = 2
    for n, label in enumerate(ax.xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)
    ax_ymin, ax_ymax = ax.get_ylim()
    ax2_ymin, ax2_ymax = ax2.get_ylim()
    if ax_ymin < ax2_ymin:
        y_min = ax_ymin
    else:
        y_min = ax2_ymin
    if ax_ymax > ax2_ymax:
        y_max = ax_ymax
    else:
        y_max = ax2_ymax
    ax.set_ylim(y_min, y_max)
    ax2.set_ylim(y_min, y_max)
    plt.savefig("plots/wilson{}.png".format(i), bbox_inches="tight")
    plt.close()
'''

#create plot for rmsd between original model and refined simulation models;
print("plotting r.m.s.d. between original model and refined simulation models...")
fig, ax = plt.subplots(1, 1, figsize=(7, 6))
x = np.linspace(0,valnum,valnum)
for i in range(4):
    ax.plot(x, np.array(rmsd[i]), label=parameter[i], marker=i+5)
ax.set_xlabel('STDDEV parameter level')
ax.set_ylabel('r.m.s.d. (Å) between model with and without added error')
ax.set_title('r.m.s.d. (Å) as a function of STDDEV parameters', fontweight="bold")
ax.legend()
plt.xticks([0, valnum], ["0", "max"])
ax.set_ylim(0, None)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(11)
plt.savefig("plots/rmsd.png", bbox_inches="tight")
plt.close()

#create plots for beam e.s.d. / reflecting range e.s.d. with  a) wavelength and cell stddev   b) beam and orientation stddev;
print("plotting e.s.d. for beam divergence and reflecting range...")
mark = ["o", "o", "s", "s"]
for par in [[0, 2], [1, 3]]:
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(14, 6))
    ax1.plot(parval[par[0]], np.array(beam_esd[par[0]]), label=parameter[par[0]], marker=mark[par[0]], color="orange")
    ax3 = ax1.twiny()
    ax3.plot(parval[par[1]], np.array(beam_esd[par[1]]), label=parameter[par[1]], marker=mark[par[1]], color="blue")
    ax2.plot(parval[par[0]], np.array(ref_esd[par[0]]), label=parameter[par[0]], marker=mark[par[0]], color="orange")
    ax4 = ax2.twiny()
    ax4.plot(parval[par[1]], np.array(ref_esd[par[1]]), label=parameter[par[1]], marker=mark[par[1]], color="blue")
    ax1.set_xlabel(parameter[par[0]])
    ax1.set_ylabel("Beam divergence e.s.d.")
    ax2.set_xlabel(parameter[par[0]])
    ax2.set_ylabel("Reflecting range e.s.d.")
    ax3.set_xlabel(parameter[par[1]])
    ax4.set_xlabel(parameter[par[1]])
    ax1.set_title("Beam divergence e.s.d. as a function of\n{x} and {y}".format(x=parameter[par[0]], y=parameter[par[1]]), fontweight="bold")
    ax2.set_title("Reflecting range e.s.d. as a function of\n{x} and {y}".format(x=parameter[par[0]], y=parameter[par[1]]), fontweight="bold")
    handles, labels = unite_legends([ax1, ax3])
    ax3.legend(handles, labels)
    ax1.set_ylim(0, None)
    ax2.set_ylim(0, None)
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax3.xaxis.label, ax3.yaxis.label] + ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax4.xaxis.label, ax4.yaxis.label] + ax4.get_xticklabels() + ax4.get_yticklabels()):
        item.set_fontsize(11)
    plt.savefig("plots/esd_{x}_{y}.png".format(x=str(par[0]), y=str(par[1])), bbox_inches="tight")
    plt.close()

#create plots for atomic B-factors of protein and water atoms;
print("plotting atomic B-factors...")
x = np.array(B_origin)
x2 = np.array(B_origin_het)
for i in range(4):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    y = np.array(B_ref[i])
    y2 = np.array(B_ref_het[i])
    coef = np.polyfit(x, y, 1)
    coef2 = np.polyfit(x2, y2, 1)
    regression = np.poly1d(coef)
    regression2 = np.poly1d(coef2)
    ic = str(coef[1]).split(".")    #ic = intercept
    ic2 = str(coef2[1]).split(".")
    ax1.plot(x, regression(x), "-k", label="slope: "+str(coef[0])[:3]+"\nintercept: "+ic[0]+"."+ic[1][:1])
    ax2.plot(x2, regression2(x2), "-k", label="slope: "+str(coef2[0])[:3]+"\nintercept: "+ic2[0]+"."+ic2[1][:1])
    ax1.plot(x, y, "bo", markersize=3)
    ax2.plot(x2, y2, "bo", markersize=3)
    ax1.set_xlabel("B-factor in ideal model")
    ax2.set_xlabel("B-factor in ideal model")
    ax1.set_ylabel("B-factor in refined simulation model")
    ax2.set_ylabel("B-factor in refined simulation model")
    ax1.set_title("Protein B-factors at max {} vs ideal model".format(parameters[i]), fontweight="bold")
    ax2.set_title("Heteroatom B-factors at max {} vs ideal model".format(parameters[i]), fontweight="bold")
    ax1.set_xlim(0, None)
    ax2.set_xlim(0, None)
    ax1.set_ylim(0, None)
    ax2.set_ylim(0, None)
    ax1.legend()
    ax2.legend()
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(11)
    plt.savefig("plots/atomic_b_{}.png".format(i), bbox_inches="tight")
    plt.close()

#create plot for CC(1/2) in highest resolution bin;
print("plotting CC(1/2) in highest resolution bin...")
fig, ax = plt.subplots(1, 1, figsize=(7, 6))
x = np.linspace(0, valnum, valnum)
mark = ["o", "s", "v", "^"]
for i in range(4):
    ax.plot(x, np.array(cc12[i]), label=parameter[i], marker=mark[i])
ax.set_xlabel("STDDEV parameter level")
ax.set_ylabel(r"CC$_{1/2}$")
ax.set_title(r"CC$_{1/2}$" + " in highest resolution bin\nas a function of STDDEV parameters", fontweight="bold")
ax.legend()
ax.xaxis.set_ticks([])
ax.set_ylim(0, None)
plt.xticks([0, valnum], ["0", "max"])
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(11)
plt.savefig("plots/cc12.png", bbox_inches="tight")
plt.close()

#create plots for geometry indicators;
print("plotting geometry indicators...")
x = np.linspace(0, valnum, valnum)
mark = ["o", "s", "^", "v"]
for j in range(2):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(14, 6))
    for i in range(4):
        if j == 0:
            ax1.plot(x, np.array(ramaz[i]), label="Rama-Z ({})".format(parameter[i]), marker=mark[i])
            ax2.plot(x, np.array(clashscore[i]), label="Clashscore ({})".format(parameter[i]), marker=mark[i])
            ax1.plot(x, np.array(mpscore[i]), label="MolProbity score ({})".format(parameter[i]), marker=str(i+1))
        elif j == 1:
            ax1.plot(x, np.array(rms_bonds[i]), label="r.m.s.d.(bonds) ({})".format(parameter[i]), marker=mark[i])
            ax2.plot(x, np.array(rms_angles[i]), label="r.m.s.d.(angles) ({})".format(parameter[i]), marker=mark[i])
    ax1.set_xlabel("STDDEV parameter level")
    ax2.set_xlabel("STDDEV parameter level")
    ax1.set_ylabel("Geometry quality indicators")
    if j == 0:
        ax1.set_title("Rama-Z and MolProbity score\nas functions of STDDEV parameters", fontweight="bold")
        ax2.set_title("Clashscore as a function\nof STDDEV parameters", fontweight="bold")
        #ax1.set_ylim(None, 0)
        ax2.set_ylim(0, None)
    elif j == 1:
        ax1.set_title("r.m.s.d.(bonds) as a function of STDDEV parameters", fontweight="bold")
        ax2.set_title("r.m.s.d.(angles) as a function of STDDEV parameters", fontweight="bold")
        ax1.set_ylim(0, 0.05)
        ax2.set_ylim(0, 2)
    ax1.legend()
    ax2.legend()
    ax.xaxis.set_ticks([])
    plt.xticks([0, valnum], ["0", "max"])
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(11)
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(11)
    plt.savefig("plots/geo{}.png".format(j), bbox_inches="tight")
    plt.close()
print("all tasks complete.")

