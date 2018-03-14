all: distributions uncertainty_RRPnfL2u_21

distributions: Model.h all_models.h distributions.cc
	g++ --std=c++11 -Wall `root-config --libs` -lMinuit `root-config --cflags` \
		distributions.cc -o distributions

uncertainty_RRPnfL2u_21: Model.h all_models.h uncertainty_RRPnfL2u_21.cc stat.h
	g++ --std=c++11 -Wall `root-config --libs` -lMinuit `root-config --cflags` \
		uncertainty_RRPnfL2u_21.cc -o uncertainty_RRPnfL2u_21

distributions: Model_RRdPL2_20.h
distributions: Model_RRdPL2u_17.h
distributions: Model_RRdPL2u_19.h
distributions: Model_RRdPqcL2u_16.h
distributions: Model_RqcRcL2qc_12.h
distributions: Model_RqcRcLqc_12.h
distributions: Model_RqcRLqc_14.h
distributions: Model_RRcdPL2u_15.h
distributions: Model_RRcdPqcL2u_14.h
distributions: Model_RRcL2qc_15.h
distributions: Model_RRcLqc_15.h
distributions: Model_RRcPL_19.h
distributions: Model_RRL_18.h
distributions: Model_RRLnf_19.h
distributions: Model_RRL2_18.h
distributions: Model_RRL2qc_17.h
distributions: Model_RRLqc_17.h
distributions: Model_RRPEu_19.h
distributions: Model_RRPL_21.h
distributions: Model_RRPL2_20.h
distributions: Model_RRPL2qc_18.h
distributions: Model_RRPL2u_19.h
distributions: Model_RRPnfL2u_21.h

uncertainty_RRPnfL2u_21 : Model_RRPnfL2u_21.h
