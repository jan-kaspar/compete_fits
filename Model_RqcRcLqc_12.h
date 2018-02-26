#include "Model.h"

class Model_RqcRcLqc_12 : public Model
{
	public:

	double B = 0.75953082;
	double s0 = 118.96503;
	double Y_1_pp = 11.904844;
	double Y_2_pip = 7.0913207;
	double eta_1 = 0.20200261;
	double eta_2 = 0.55546195;

	Model_RqcRcLqc_12()
	{
		name = "Model_RqcRcLqc_12";
		label = "R^{qc} R_{c} L^{qc}\\ (12)";
	}

	double si_p_p(double s) const override
	{
		return 9. * B * log(s/s0) + 9. * Y_1_pp * pow(s, -eta_1) - 5. * Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9. * B * log(s/s0) + 9. * Y_1_pp * pow(s, -eta_1) + 5. * Y_2_pip * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9. * M_PI * B / 2. - 9. * Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - 5. * Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9. * M_PI * B / 2. - 9. * Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + 5. * Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
