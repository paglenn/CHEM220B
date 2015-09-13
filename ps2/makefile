# makefile for prog with gsl 
field: 
	g++ liq3-field.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o field -O3
hist: 
	g++ liq3-field-hist.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o hist -O3
clean: 
	./clean.sh 
	make field 
	make hist 
