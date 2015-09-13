# makefile for prog with gsl 
probe:  
	g++ liq3-probe.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o probe -O3

mfp: 
	g++ liq3-MeanFreePath.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o mfp -O3

ecmc: 
	g++ liq3-ecmc.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o ecmc -O3

field: 
	g++ liq3-field.cpp -I /usr/local/Cellar/gsl/1.16/include -L /usr/local/Cellar/gsl/1.16/lib -lgsl -lm -o field -O3
clean: 
	rm 
	make 
