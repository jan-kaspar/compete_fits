all: distributions

distributions: Model.h distributions.cc
	g++ --std=c++11 -Wall `root-config --libs` -lMinuit `root-config --cflags` \
		distributions.cc -o distributions

distributions: Model_RPdPL2_20.h
distributions: Model_RPdPL2u_17.h
distributions: Model_RPdPL2u_19.h
distributions: Model_RPdPqcL2u_16.h

distributions: Model_RqcRcLqc_12.h
distributions: Model_RqcRcL2qc_12.h
