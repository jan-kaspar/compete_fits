#include "Model.h"

class Model_RRL2qc_17 : public Model
{
	public:

	double B = 0.016076625;
	double s0 = 0.00044830086;
	double Y_1_pp = 67.01912;
	double Y_2_pp = 35.528122;
	double eta_1 = 0.27261187;
	double eta_2 = 0.5552468;

	Model_RRL2qc_17()
	{
		name = "Model_RRL2qc_17";
		label = "RR L2^{qc}\\ (17)";
	}

	double si_p_p(double s) const override
	{
		return 9.*B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI*B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI*B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
