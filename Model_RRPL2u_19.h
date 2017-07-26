#include "Model.h"

class Model_RRPL2u_19 : public Model
{
	public:

	double Z_pp = 35.865711;
	double B = 0.31573;
	double s0 = 34.409806;
	double Y_1_pp = 42.069653;
	double Y_2_pp = 32.155544;
	double eta_1 = 0.46822531;
	double eta_2 = 0.53956628;

	Model_RRPL2u_19()
	{
		name = "Model_RRPL2u_19";
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
