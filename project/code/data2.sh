#!/bin/bash
g++ proj.cpp VEC.cpp Complex.cpp
for i in 1
do
    for j in a b c d e f 
    do
        ./a.out ../data/t$j$i.dat
    done
done
