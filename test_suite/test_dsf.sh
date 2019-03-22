#!/bin/bash

for fg in 1 2 4 8 16 32; do
    #dynamic_structure_factor traj.xyz --max 0 --fg $fg
    #mv dsf_00000.dat bf_dsf_${fg}.dat
    mv bf_dsf_${fg}.dat dsf_${fg}.dat
done

