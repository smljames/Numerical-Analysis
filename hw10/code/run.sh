#!/bin/bash

g++ hw10.cpp MAT.cpp VEC.cpp

for i in 12 24 48 96 192 384 768 1536
do
    for j in 1 2 3 4 5 6
    do
        ./a.out 0 2 $i $j
    done
done
