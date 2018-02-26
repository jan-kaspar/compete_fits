#include "Model.h"

class Model_RRL_18 : public Model
{
	public:

	double B_pp = 6.6282714;
	double s0 = 79.371146;
	double Y_1_pp = 104.68182;
	double Y_2_pp = 33.233187;
	double eta_1 = 0.21134153;
	double eta_2 = 0.54379057;

	Model_RRL_18()
	{
		name = "Model_RRL_18";
		label = "RRL\\ (18)";
	}

	double si_p_p(double s) const override
	{
		return B_pp * log(s/s0) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return B_pp * log(s/s0) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
