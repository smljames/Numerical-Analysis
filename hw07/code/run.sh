#!/bin/bash
g++ hw07.cpp MAT.cpp VEC.cpp
rm ../outputs/time.txt

for i in {3..8}
do
    echo m"$i".dat | tee -a ../outputs/time.txt
    ./a.out < m"$i".dat | tee -a ../outputs/time.txt
done
