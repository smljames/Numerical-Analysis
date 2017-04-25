#!/bin/bash
rm ../output/output.txt
for i in 2 4 10 20 40
do
    g++ -o hw06_$i hw06_$i.cpp MAT.cpp VEC.cpp
    ./hw06_$i $i >> ../output/output.txt
done
