#!/bin/bash
./a.out 2 |tee -a run_time.txt
./a.out 4 |tee -a run_time.txt
./a.out 10 |tee -a run_time.txt
./a.out 20 |tee -a run_time.txt
./a.out 40 |tee -a run_time.txt
./a.out 50 |tee -a run_time.txt
