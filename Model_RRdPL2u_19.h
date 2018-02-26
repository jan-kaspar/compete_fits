#include "Model.h"

class Model_RRdPL2u_19 : public Model
{
	public:

	double Z_pp = 37.046665;
	double B_pp = 0.32765691;
	double s0 = 49.056427;
	double Y_1_pp = 44.317753;
	double Y_2_pp = 30.819292;
	double eta = 0.53076334;

	Model_RRdPL2u_19()
	{
		name = "Model_RRdPL2u_19";
		label = "(RR)^d PL2_u\\ (19)";
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
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) - Y_2_pp * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) + Y_2_pp * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
