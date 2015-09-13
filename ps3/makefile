# makefile for prog with gsl 
field: 
	g++ liq3-pairs.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o field -O3
clean: 
	./clean.sh 
	make field 
