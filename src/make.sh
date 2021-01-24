!#/bin/bash

cp ./NtCForces.cpp /github/MMB/src/

cd /github/MMB/build/

make clean
make -j4

cp MMB /usr/bin/MMB2

cd
