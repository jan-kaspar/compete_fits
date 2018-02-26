#include "Model.h"

class Model_RqcRcL2qc_12 : public Model
{
	public:

	double B = 0.016080647;
	double s0 = 0.00045002716;
	double Y_1_pp = 7.4402016;
	double Y_2_pip = 7.1370767;
	double eta_1 = 0.2722584;
	double eta_2 = 0.55494024;

	Model_RqcRcL2qc_12()
	{
		name = "Model_RqcRcL2qc_12";
		label = "R^{qc} R_{c} L2^{qc}\\ (12)";
	}

	double si_p_p(double s) const override
	{
		return 9. * B * pow(log(s/s0), 2.) + 9. * Y_1_pp * pow(s, -eta_1) - 5. * Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9. * B * pow(log(s/s0), 2.) + 9. * Y_1_pp * pow(s, -eta_1) + 5. * Y_2_pip * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9. * M_PI * B * log(s/s0) - 9. * Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - 5. * Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9. * M_PI * B * log(s/s0) - 9. * Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + 5. * Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
