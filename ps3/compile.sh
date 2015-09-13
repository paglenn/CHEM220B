#!/bin/bash 

g++ liq3-pairs.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o pairs -O3
#./pairs 
