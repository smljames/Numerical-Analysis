#!/bin/bash
for i in 1 2 3
do
    g++ hw04_"$i"_norm.cpp MAT_"$i"_norm.cpp VEC.cpp
    ./a.out 40 | tee -a time_40.txt
done
