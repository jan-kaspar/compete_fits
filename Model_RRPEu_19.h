#include "Model.h"

class Model_RRPEu_19 : public Model
{
	public:

	double Z_pp = 4.1824033;
	double X = 15.685953;
	double ep = 0.10134777;
	double Y_1_pp = 57.450783;
	double Y_2_pp = 33.468752;
	double eta_1 = 0.34309186;
	double eta_2 = 0.54487001;

	Model_RRPEu_19()
	{
		name = "Model_RRPEu_19";
	}

	double si_p_p(double s) const override
	{
		return Z_pp + X * pow(s, ep) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + X * pow(s, ep) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = - X*pow(s, ep) / tan((1.+ep)/2*M_PI) - Y_1_pp*pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - Y_2_pp*pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = - X*pow(s, ep) / tan((1.+ep)/2*M_PI) - Y_1_pp*pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + Y_2_pp*pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
