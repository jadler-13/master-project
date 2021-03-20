#!/bin/bash

###phenix.ready_set replace_this_b_unchanged.pdb
###rm TLA* replace_this_b_unchanged.eff

echo generating ensemble model, this may take a few hours...
cd ref || exit
sh run1.sh
cd ..

#change column labels in refine.mtz to make it usable for phenix.refine;
rm refine.mtz
sftools <<EOF
read replace_this_b_unchanged.pdb.mtz
calc F col FOBS = col FMODEL
calc Q col SIGFOBS = 1
write refine.mtz
quit
EOF
rm replace_this_b_unchanged.pdb.mtz

rm -r ensemble_files
mkdir ensemble_files || exit
cd ensemble_files || exit

for f in {0..19}
do
    phenix.pdbtools ../replace_this_b_unchanged.pdb output.file_name="$f.pdb" modify.sites.shake=1.0 || exit 
    phenix.refine "$f.pdb" ../refine.mtz allow_polymer_cross_special_position=True pdb_interpretation.max_reasonable_bond_distance=200 refinement.input.xray_data.r_free_flags.generate=True refinement.input.xray_data.r_free_flags.ignore_r_free_flags=True simulated_annealing=True refinement.main.number_of_macro_cycles=20 ordered_solvent=True --overwrite || exit
done

iotbx.pdb.join_fragment_files *_refine_001.pdb > ../replace_this_ensemble.pdb || exit
cd ..

echo ensemble-model created.
