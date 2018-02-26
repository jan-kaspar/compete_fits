#include "Model.h"

class Model_RqcRLqc_14 : public Model
{
	public:

	double B = 0.75972486;
	double s0 = 119.34367;
	double Y_1_pp = 11.907349;
	double Y_2_pp = 35.454203;
	double eta_1 = 0.20193472;
	double eta_2 = 0.55543486;

	Model_RqcRLqc_14()
	{
		name = "Model_RqcRLqc_14";
		label = "R^{qc} R L^{qc}\\ (14)";
	}

	double si_p_p(double s) const override
	{
		return 9.*B * log(s/s0) + 9.*Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*B * log(s/s0) + 9.*Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI*B/2. - 9.*Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI*B/2. - 9.*Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
