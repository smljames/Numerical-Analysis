#!/bin/bash
for((i=3; i<11; i++))
do 
    ./a.out m$i.dat |tee -a run_time.txt
done
