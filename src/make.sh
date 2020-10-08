#!/bin/bash

cp ./NtCForces.cpp /github/MMB/src/

cd /github/MMB/build/

make -j4

cp MMB /usr/bin/

cd
