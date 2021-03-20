#!/bin/bash

#calculate structure factors from electron density;
phenix.fmodel ../replace_this_b_unchanged.pdb high_resolution=1.6 k_sol=0.362 b_sol=72.5 > fmodellog
cp replace_this_b_unchanged.pdb.mtz ../

#determine space group number and change it in SIMULATE.INP;
mtzdmp replace_this_b_unchanged.pdb.mtz > spacelog
python3.6 spacegroup.py
mv temp_SIMULATE.INP SIMULATE.INP
mv ../new/temp_SIMULATE.INP ../new/SIMULATE.INP

#expand reflections to spacegroup P1;
phenix.reflection_file_converter --expand_to_p1 replace_this_b_unchanged.pdb.mtz --mtz=temp-p1.pdb.mtz > rfclog

#write reflection data into a log file that can easily be edited with awk;
mtzdump hklin temp-p1.pdb.mtz<<EOF > log 
nref -1
go
EOF

#extract lines with reflection info from log;
awk '/LIST OF/,/FONT COLOR/ {print $0}' log | tail -n +4 | grep -v FONT > temp.hkl
awk '{print $1, $2, $3, $4, $5}' temp.hkl > tempsmall.hkl

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
sim_mx.64 -n 10000 -R -r $(date +%N) > simlog

#analyze simulated frames, calculate ISa, write beam divergence and reflecting range to temporary file esd;
cd xds || exit
rm XDS_ASCII.HKL
xds_par > xdslog
grep -i -A 1 "beam_divergence_e.s.d.= " INTEGRATE.LP > esd
grep -i -A 11 "cc(1/2)" CORRECT.LP | tail -1 > cc12

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

#create reference data set for correct indexing;
pointless XYZIN ../../replace_this_b_unchanged.pdb HKLOUT temp.mtz XDS_ASCII.mtz > pointlesslog
#square model amplitudes:
sftools <<EOF
read temp.mtz
select col FP = present
calc col I-model = col FP col FP *
write I-model.mtz col I-model
quit
EOF
#dump to ASCII format:
mtzdump hklin I-model.mtz > temp.hkl <<EOF
nref -1
end
EOF
#prepare I-model.hkl:
echo \!FORMAT=XDS_ASCII > I-model.hkl
echo \!END_OF_HEADER   >> I-model.hkl
echo pick reflection info from temp.hkl:
awk '/LIST OF REFLECTIONS/,/<B><FONT COLOR=/' temp.hkl | tail -n +4 | head -n -1 | awk '{print $0,1}' >> I-model.hkl
echo \!END_OF_DATA   >> I-model.hkl
#I-model.hkl is now ready to be used as REFERENCE_DATA_SET
rm temp.hkl I-model.mtz
mv I-model.hkl ../../indexing_reference.hkl
