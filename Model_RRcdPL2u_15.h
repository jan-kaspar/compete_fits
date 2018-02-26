#include "Model.h"

class Model_RRcdPL2u_15 : public Model
{
	public:

	double Z_pp = 36.784336;
	double B = 0.31175883;
	double s0 = 39.731955;
	double Y_1_pp = 48.60371;
	double Y_2_pip = 6.5969507;
	double eta = 0.54348946;

	Model_RRcdPL2u_15()
	{
		name = "Model_RRcdPL2u_15";
		label = "(RR_{c})^{d} PL2_{u}\\ (15)";
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B * pow(log(s/s0), 2) + (Y_1_pp - 5.*Y_2_pip) * pow(s, -eta);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B * pow(log(s/s0), 2) + (Y_1_pp + 5.*Y_2_pip) * pow(s, -eta);
		return 0;
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) - 5.*Y_2_pip * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) + 5.*Y_2_pip * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_ap(s);
		return 0;
	}
};
