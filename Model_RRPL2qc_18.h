#include "Model.h"

class Model_RRPL2qc_18 : public Model
{
	public:

	double A = -3.7258004;
	double B = 0.79164674;
	double C = -0.0018160896;
	double Y_1_pp = 109.05955;
	double Y_2_pp = 33.077245;
	double eta_1 = 0.20641782;
	double eta_2 = 0.5428034;

	Model_RRPL2qc_18()
	{
		name = "Model_RRPL2qc_18";
	}

	double si_p_p(double s) const override
	{
		return 9.*(A + B*log(s) + C*pow(log(s), 2)) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*(A + B*log(s) + C*pow(log(s), 2)) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI * (B/2. + C*log(s)) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
		return 0;
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI * (B/2. + C*log(s)) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
