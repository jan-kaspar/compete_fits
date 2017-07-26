#include "Model.h"

class Model_RRLqc_17 : public Model
{
	public:

	/*
	double Z_pp = 37.276058;
	double B_pp = 0.33436876;
	double s0 = 55.768319;
	double Y_1_pp = 42.915501;
	double Y_2_pp = 31.143403;
	double eta = 0.53323572;
	*/

	Model_RRLqc_17()
	{
		name = "Model_RRLqc_17";
	}

	double si_p_p(double s) const override
	{
		//return Z_pp + B_pp * pow(log(s/s0), 2) + (Y_1_pp - Y_2_pp) * pow(s, -eta);
		return 0;
	}

	double si_p_ap(double s) const override
	{
		//return Z_pp + B_pp * pow(log(s/s0), 2) + (Y_1_pp + Y_2_pp) * pow(s, -eta);
		return 0;
	}

	double rho_p_p(double s) const override
	{
		//double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2*M_PI) - Y_2_pp * pow(s, -eta) * tan((1.-eta)/2*M_PI);
		//return rho_si / si_p_p(s);
		return 0;
	}

	double rho_p_ap(double s) const override
	{
		//double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2*M_PI) + Y_2_pp * pow(s, -eta) * tan((1.-eta)/2*M_PI);
		//return rho_si / si_p_ap(s);
		return 0;
	}
};
