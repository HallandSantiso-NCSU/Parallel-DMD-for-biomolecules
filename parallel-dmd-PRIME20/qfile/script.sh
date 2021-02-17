#!/bin/bash

source /usr/local/bin/setupics

cd ../code

mpif90 -Os -xAVX -no-prec-div -r8 -arch host -align dcommons -g -traceback -Dchaptype=1 -Dnop1=672 -Dnop2=672 -Dchnln1=7 -Dchnln2=7 -Dnumbeads1=28 -Dnumbeads2=28 -Dnumbin=2000 -Dn_wrap=2 -Ddebugging -Dcanon -Drunr -o ../dmd main.F90

cd ..

temps='050 045 040 035 030 028 026 024 022' 
for i in $temps
do ./dmd < temp_$i > out_$i wait
done

for i in {1..20}
do ./dmd < temp_018 > out_018_$i wait
done
