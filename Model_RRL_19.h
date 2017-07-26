#include "Model.h"

class Model_RRL_19 : public Model
{
	public:

	double A = -30.265138;
	double B = 6.7106141;
	double Y_1_pp = 105.82114;
	double Y_2_pp = 33.358907;
	double eta_1 = 0.20882981;
	double eta_2 = 0.54453128;

	Model_RRL_19()
	{
		name = "Model_RRL_19";
	}

	double si_p_p(double s) const override
	{
		return A + B * log(s) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return A + B * log(s) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI*B/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI*B/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
