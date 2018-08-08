#!/bin/bash                                                                                                                                                                          

# Command to run ARGweaver
arg-sample -s mytest2/mytest2.sites -N 10000 -r .075e-10 -m 1e-9 --ntimes 40 --maxtime 200e3 -c 1 -n 10000 -o mytest2/mytest2.sample2/out --overwrite

# Commands to run Arbores (note that for the second line you need to create the initialisation file)
./Arbores arb_data4.txt 200000 1e-8 0.001 1 ms_test_8_v1
./Arbores arb_data4.txt 200000 1e-8 0.001 10 ms_test_8_v2 init_8_v2.txt
