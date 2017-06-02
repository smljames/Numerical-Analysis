#!/bin/bash
g++ proj.cpp VEC.cpp Complex.cpp
for i in f
do
    for j in 1 2 3 4 5
    do
        ./a.out ../data/t$i$j.dat
    done
done
