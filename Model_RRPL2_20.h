#include "Model.h"

class Model_RRPL2_20 : public Model
{
	public:

	double A = -38.50015;
	double B = 7.5860589;
	double C = -0.02800196;
	double Y_1_pp = 113.76647;
	double Y_2_pp = 33.168757;
	double eta_1 = 0.20020869;
	double eta_2 = 0.54334494;

	Model_RRPL2_20()
	{
		name = "Model_RRPL2_20";
	}

	double si_p_p(double s) const override
	{
		return A + B * log(s) + C * pow(log(s), 2) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return A + B * log(s) + C * pow(log(s), 2) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * (B/2. + C * log(s)) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * (B/2. + C * log(s)) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
		return 0;
	}
};
