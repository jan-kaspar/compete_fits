#include "Model.h"

class Model_RPdPqcL2u_16 : public Model
{
	public:

	double Z = 4.1266114;
	double B = 0.31640911;
	double s0 = 46.272116;
	double Y_1_pp = 46.412235;
	double Y_2_pp = 34.157737;
	double eta = 0.55062892;

	Model_RPdPqcL2u_16()
	{
		name = "Model_RPdPqcL2u_16";
	}

	double si_p_p(double s) const override
	{
		return 9.*Z + B * pow(log(s/s0), 2) + (Y_1_pp - Y_2_pp) * pow(s, -eta);
	}

	double si_p_ap(double s) const override
	{
		return 9.*Z + B * pow(log(s/s0), 2) + (Y_1_pp + Y_2_pp) * pow(s, -eta);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2*M_PI) - Y_2_pp * pow(s, -eta) * tan((1.-eta)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2*M_PI) + Y_2_pp * pow(s, -eta) * tan((1.-eta)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
