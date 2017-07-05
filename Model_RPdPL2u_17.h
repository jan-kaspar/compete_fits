#include "Model.h"

class Model_RPdPL2u_17 : public Model
{
	public:

	double Z_pp = 36.863204;
	double B_pp = 0.32009024;
	double s0 = 43.88116;
	double Y_1_pp = 45.42289;
	double Y_2_pp = 30.53598;
	double eta = 0.52949847;

	Model_RPdPL2u_17()
	{
		name = "Model_RPdPL2u_17";
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B_pp * pow(log(s/s0), 2) + (Y_1_pp - Y_2_pp) * pow(s, -eta);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B_pp * pow(log(s/s0), 2) + (Y_1_pp + Y_2_pp) * pow(s, -eta);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2*M_PI) - Y_2_pp * pow(s, -eta) * tan((1.-eta)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2*M_PI) + Y_2_pp * pow(s, -eta) * tan((1.-eta)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
