scattering:	scattering.C
	g++ -std=c++11 -O3 -o scattering scattering.C `root-config --cflags --glibs` `psrchive --cflags --libs`

scattering8:	scattering_nchn8.C
	g++ -Wall -std=c++11 -O3 -o scattering8 scattering_nchn8.C `root-config --cflags --glibs` `psrchive --cflags --libs`

scatternest:    scatternest.C
	mpiicpc -c Parameters.C utils.C scatter_likelihood_MN.C get_scatter_chi.C read_results.C scatter_likelihood_PC.C c_interface.cpp `pkg-config --cflags glib-2.0` `pkg-config --libs glib-2.0`
	mpiicpc -DHAVE_POLYCHORD -Wall -std=c++11 -O3 -o scatternest scatternest.C Parameters.o utils.o scatter_likelihood_MN.o get_scatter_chi.o read_results.o scatter_likelihood_PC.o c_interface.o `psrchive --cflags --libs` -L$(MULTINEST_DIR) -lnest3 `pkg-config --libs glib-2.0` -L/u/gdesvign/src/PolyChord_v1.9/lib -lchord
