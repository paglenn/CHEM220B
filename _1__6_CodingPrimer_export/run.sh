#!/bin/sh
NAME=cpp_prog
g++ hello.cc -o $NAME -lgsl
./$NAME 10

