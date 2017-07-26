#include "Model.h"

class Model_RRL2_18 : public Model
{
	public:

	double B_pp = 0.14468157;
	double s0 = 0.00044828241;
	double Y_1_pp = 67.016128;
	double Y_2_pp = 35.55018;
	double eta_1 = 0.27258771;
	double eta_2 = 0.55536716;

	Model_RRL2_18()
	{
		name = "Model_RRL2_18";
	}

	double si_p_p(double s) const override
	{
		return B_pp * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return B_pp * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
