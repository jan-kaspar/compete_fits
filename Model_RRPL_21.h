#include "Model.h"

class Model_RRPL_21 : public Model
{
	public:

	double Z_pp = -28.963728;
	double B_pp = 6.5817218;
	double Y_1_pp = 103.61068;
	double Y_2_pp = 32.036043;
	double eta_1 = 0.20804983;
	double eta_2 = 0.53637805;

	Model_RRPL_21()
	{
		name = "Model_RRPL_21";
		label = "RR PL\\ (21)";
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B_pp * log(s) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B_pp * log(s) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI*B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI*B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
