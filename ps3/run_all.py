#!/usr/bin/env python
import os

for p in [0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    os.system("./pairs %f"%p)
os.system("mv density_hist-0.0.dat density_hist-0.02.dat")
os.system("./plot.py")
