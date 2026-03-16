#!/bin/bash
set -e 

#cd VFscan
ind=./test_files/Rick22
oud=./test_oud
refd=./database/VFDB_setB_nt # setB
bash batch_blastn_setB.sh -i $ind -o $oud -r $refd

