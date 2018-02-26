#include "Model.h"

class Model_RRcLqc_15 : public Model
{
	public:

	double B = 0.74989124;
	double s0 = 101.65274;
	double Y_1_pp = 106.03196;
	double Y_2_pip = 7.0873701;
	double eta_1 = 0.20535026;
	double eta_2 = 0.5554147;

	Model_RRcLqc_15()
	{
		name = "Model_RRcLqc_15";
		label = "RR_c L^{qc}\\ (15)";
	}

	double si_p_p(double s) const override
	{
		return 9.*B * log(s/s0) + Y_1_pp * pow(s, -eta_1) - 5.*Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*B * log(s/s0) + Y_1_pp * pow(s, -eta_1) + 5.*Y_2_pip * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI*B/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI*B/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
