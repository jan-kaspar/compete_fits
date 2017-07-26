#include "Model.h"

class Model_RRcdPqcL2u_14 : public Model
{
	public:

	double Z = 4.1303189;
	double B = 0.31752323;
	double s0 = 47.118111;
	double Y_1_pp = 46.318568;
	double Y_2_pip = 6.8915267;
	double eta = 0.55122458;

	Model_RRcdPqcL2u_14()
	{
		name = "Model_RRcdPqcL2u_14";
	}

	double si_p_p(double s) const override
	{
		return 9.*Z + B * pow(log(s/s0), 2) + (Y_1_pp - 5.*Y_2_pip) * pow(s, -eta);
	}

	double si_p_ap(double s) const override
	{
		return 9.*Z + B * pow(log(s/s0), 2) + (Y_1_pp + 5.*Y_2_pip) * pow(s, -eta);
		return 0;
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2*M_PI) - 5.*Y_2_pip * pow(s, -eta) * tan((1.-eta)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2*M_PI) + 5.*Y_2_pip * pow(s, -eta) * tan((1.-eta)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
