#!/bin/bash

#calculate structure factors from electron density;
phenix.fmodel ../replace_this_ensemble.pdb high_resolution=1.79 k_sol=0.362 b_sol=72.5 > fmodellog

#expand reflections to spacegroup P1;
phenix.reflection_file_converter --expand_to_p1 replace_this_ensemble.pdb.mtz --mtz=temp-p1.pdb.mtz > rfclog

#write reflection data into a log file that can easily be edited with awk;
mtzdump hklin temp-p1.pdb.mtz<<EOF > log 
nref -1
go
EOF

#extract lines with reflection info from log;
awk '/LIST OF/,/FONT COLOR/ {print $0}' log | tail -n +4 | grep -v FONT > temp.hkl
awk '{print $1, $2, $3, $4/20, $5}' temp.hkl > tempsmall.hkl

#fort.1 is the expected input file for fullsphere from P1, hence this link;
ln -sfn tempsmall.hkl fort.1

#write all Friedel mates into the file because SIM_MX has no built-in space-group library and can't perform symmetry expansion;
echo 4 4 | fullsphere_from_P1

#fort.2 is the default output name for fullsphere_from_P1, intensities.hkl is a default input for sim_mx;
mv fort.2 intensities.hkl

#cleanup;
rm fort.1 log temp.hkl temp-p1.pdb.mtz

#simulate diffraction experiment; -n specifies number of simulated 'mosaic blocks', -r determines seed for random number for photon count and background noise;
rm -r xds/image
mkdir xds/image
sim_mx.64 -n 100000 -R -r $(date +%N) > simlog

#break1

#analyze simulated frames, calculate ISa, write beam divergence and reflecting range to temporary file esd;
cd xds || exit
rm XDS_ASCII.HKL
xds_par > xdslog
grep -i -A 1 "beam_divergence_e.s.d.= " INTEGRATE.LP > esd
grep -i -A 11 "cc(1/2)" CORRECT.LP | tail -1 > cc12

#break2

#convert reflection data into mtz format to make it readable for uniqueify;
echo "INPUT_FILE= XDS_ASCII.HKL" > XDSCONV.INP
echo "OUTPUT_FILE= temp.hkl CCP4" >> XDSCONV.INP
echo "FRIEDEL'S_LAW= TRUE" >> XDSCONV.INP
xdsconv > xdsconvlog
f2mtz HKLOUT temp.mtz<F2MTZ.INP > f2mtzlog
cad HKLIN1 temp.mtz HKLOUT XDS_ASCII.mtz<<EOF > cadlog
LABIN FILE 1 ALL
END
EOF
rm XDSCONV.INP temp.hkl temp.mtz

#create unique R free flags if not yet available and extends them to cover all reflections;
cd ..
if ! [ -f uniquerfree.mtz ]
then
uniqueify xds/XDS_ASCII.mtz uniquetemp.mtz > uniqlog
sftools <<EOF > sftoolslog
read uniquetemp.mtz
write uniquerfree.mtz col FreeR_flag
Y
stop
EOF
rm uniquetemp.mtz uniquetemp.log
fi
cad hklin1 uniquerfree.mtz hklin2 xds/XDS_ASCII.mtz hklout HKLOUT.mtz<<EOF > cadlog2
LABIN FILE 1 ALL
LABIN FILE 2 ALL
EOF
rm reflections_*.mtz uniquerfree.mtz reflections_*.log
uniqueify -f FreeR_flag HKLOUT.mtz reflections_1.mtz > extendlog

#create and refine a new model and electron density, calculate R-factors, mean protein B-factor, and r.m.s.d. between original and refined model;
cd refine || exit
phenix.refine ../../replace_this.pdb ../reflections_1.mtz  sites.shake=0.2 strategy=individual_sites+individual_adp main.number_of_macro_cycles=6 xray_data.high_resolution=1.79 refinement.output.n_resolution_bins=10 --overwrite > reflog
tail -2 reflog > rs
superpose ../../replace_this_reference.pdb replace_this_refine_001.pdb > rmsd

#save geometry quality indicators;
phenix.molprobity replace_this_refine_001.pdb > molprobity.out

#calculate wilson B factor for F-obs and F-model;
wilson hklin replace_this_refine_001.mtz <<EOF >wilsonlog1
CONTENTS H 413 C 249 N 64 O 71 S 6
NRESIDUE 56
LABIN FP=F-obs SIGFP=SIGF-obs
RSCALE 3 1.79
END
EOF
wilson hklin replace_this_refine_001.mtz <<EOF >wilsonlog2
CONTENTS H 413 C 249 N 64 O 71 S 6
NRESIDUE 56
LABIN FP=F-model SIGFP=PHIF-model
RSCALE 3 1.79
END
EOF

#create reflection file with S**2 for wilson plot;
mtzutils hklin replace_this_refine_001.mtz hklout wilson.mtz <<EOF > utilslog
INCLUDE F-obs F-model
RUN
EOF
mtzdump hklin wilson.mtz <<EOF > wrefs
RESO 1.79 50
LRESO
NREF -1
EOF
sort -nk4 wrefs > ../wilsonrefs
