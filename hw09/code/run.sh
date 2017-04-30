#!/bin/bash
g++ hw09.cpp MAT.cpp VEC.cpp

for i in 3 5 7 13 21
do
    ./a.out f"$i".dat > ../output/f"$i".txt
done
