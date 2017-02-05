scattering:	scattering.C
	g++ -std=c++11 -O3 -o scattering scattering.C `root-config --cflags --glibs` `psrchive --cflags --libs`

scattering8:	scattering_nchn8.C
	g++ -Wall -std=c++11 -O3 -o scattering8 scattering_nchn8.C `root-config --cflags --glibs` `psrchive --cflags --libs`

scatternest:	scatternest.C
	mpiicpc -c Parameters.C utils.C scatter_likelihood.C `pkg-config --cflags glib-2.0` `pkg-config --libs glib-2.0`
	mpiicpc -Wall -std=c++11 -O3 -o scatternest scatternest.C Parameters.o utils.o scatter_likelihood.o `psrchive --cflags --libs` -L$(MULTINEST_DIR) -lnest3 `pkg-config --libs glib-2.0`

scatternestX:	scatternestX.C
	mpiicpc -c Parameters.C utils.C scatter_likelihood2.C `pkg-config --cflags glib-2.0` `pkg-config --libs glib-2.0`
	mpiicpc -Wall -std=c++11 -O3 -o scatternestX scatternestX.C Parameters.o utils.o scatter_likelihood2.o `psrchive --cflags --libs` -L$(MULTINEST_DIR) -lnest3 `pkg-config --libs glib-2.0`
