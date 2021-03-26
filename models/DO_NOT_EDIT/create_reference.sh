#!/bin/bash

echo generating reference model, this may take a few minutes...
cd ref
sh run2.sh
mv refine/replace_this_refine_001.pdb ../replace_this_reference.pdb
cd ..
