#!/bin/bash
g++ hw05.cpp MAT.cpp VEC.cpp
for i in 20 40 60 80 100
do
    time ./a.out $i | tee -a time.txt
done
