#include "Model.h"

class Model_RRLqc_17 : public Model
{
	public:

	double B = 0.73635067;
	double s0 = 79.027935;
	double Y_1_pp = 104.66237;
	double Y_2_pp = 33.205188;
	double eta_1 = 0.21146055;
	double eta_2 = 0.54360032;

	Model_RRLqc_17()
	{
		name = "Model_RRLqc_17";
	}

	double si_p_p(double s) const override
	{
		return 9.*B*log(s/s0) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*B*log(s/s0) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI*B/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI*B/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
