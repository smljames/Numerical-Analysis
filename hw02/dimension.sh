#!/bin/bash
for((i=3; i<11; i++))
do
    cat m$i.dat | head -n 1
done
