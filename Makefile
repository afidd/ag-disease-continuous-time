# Requires:
#   BOOST libraries: boost.org
#   Gnu scientific libraries (GSL): https://www.gnu.org/software/gsl
#   HDF5 library: http://www.hdfgroup.org/HDF5/
#   Semi-Markov library: https://github.com/afidd/Semi-Markov
#

BOOST=/home/ajd27/Documents/boost_1_57_0
# Different Boost installations have different suffixes.
# If there is no suffix, use "BOOSTVARIANT=".
BOOSTVARIANT=
SEMIMARKOV=/usr/local/include/semimarkov-0.1
HDF5=/usr/local/hdf5-1.8.11

CXX=g++
# -DSMVHIDELOG -pg
OPT=-g -O2 -DSMVHIDELOG
INCLUDES=-I$(SEMIMARKOV) -I. -I$(BOOST)/include -I$(HDF5)/include 
LIBS=-L$(BOOST)/lib -L$(HDF5)/lib  \
    -lboost_unit_test_framework$(BOOSTVARIANT) \
	-lboost_log_setup$(BOOSTVARIANT) -lboost_log$(BOOSTVARIANT) \
	-lboost_chrono$(BOOSTVARIANT) -lboost_thread$(BOOSTVARIANT) \
	-lboost_date_time$(BOOSTVARIANT) -lboost_filesystem$(BOOSTVARIANT) \
	-lboost_program_options$(BOOSTVARIANT) -lboost_random$(BOOSTVARIANT) \
	-lboost_system$(BOOSTVARIANT) -lhdf5 -lhdf5_hl -lgsl -lgslcblas \
	-lpthread


contact: sir_exp.o main.o hdf_file.o
	g++ $(OPT) -fPIC -o sirexp sir_exp.o main.o hdf_file.o $(LIBS)

sir_exp.o: sir_exp.cpp sir_exp.hpp
	g++ sir_exp.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o sir_exp.o

hdf_file.o: hdf_file.cpp hdf_file.hpp sir_exp.hpp
	g++ hdf_file.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o hdf_file.o

main.o: main.cpp contact_version.hpp sir_exp.hpp
	g++ main.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o main.o

contact_version.hpp: Makefile
	python getgit.py contact_version.hpp

clean:
	rm -f *.o contact contact_version.hpp
