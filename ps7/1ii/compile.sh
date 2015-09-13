#!/bin/bash 

g++  simulation.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o dynamics -O3
#./pairs 
