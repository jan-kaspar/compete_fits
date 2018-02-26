#include "Model.h"

class Model_RRcL2qc_15 : public Model
{
	public:

	double B = 0.016178984;
	double s0 = 0.00047862371;
	double Y_1_pp = 67.017688;
	double Y_2_pip = 7.1656375;
	double eta_1 = 0.2714367;
	double eta_2 = 0.55573821;

	Model_RRcL2qc_15()
	{
		name = "Model_RRcL2qc_15";
		label = "RR_{\\rm c} L2^{\\rm qc}\\ (15)";
	}

	double si_p_p(double s) const override
	{
		return 9.*B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) - 5.*Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) + 5.*Y_2_pip * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
