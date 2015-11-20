graficar.py: data.dat
	python graficar.py
data.dat: mcmc_solar.x
	./mcmc_solar.x 150000 20000 > data.dat
mcmc_solar.x: mcmc_solar.c
	 cc mcmc_solar.c -lgsl -lgslcblas -lm -o mcmc_solar.x


clean:
	rm data.dat 
