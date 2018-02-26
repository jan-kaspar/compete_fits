#include "Model.h"

class Model_RRcPL_19 : public Model
{
	public:

	double Z_pp = -27.389097;
	double B_pp = 6.4776681;
	double Y_1_pp = 101.95344;
	double Y_2_pip = 6.1769888;
	double eta_1 = 0.21077046;
	double eta_2 = 0.53045536;

	Model_RRcPL_19()
	{
		name = "Model_RRcPL_19";
		label = "RR_c PL\\ (19)";
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B_pp * log(s) + Y_1_pp * pow(s, -eta_1) - 5.*Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B_pp * log(s) + Y_1_pp * pow(s, -eta_1) + 5.*Y_2_pip * pow(s, -eta_2);
		return 0;
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
