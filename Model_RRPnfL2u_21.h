#include "Model.h"

class Model_RRPnfL2u_21 : public Model
{
	public:

	double Z_pp;
	double B;
	double s0;
	double Y_1_pp;
	double Y_2_pp;
	double eta_1;
	double eta_2;

	virtual void SetDefaultParameterValues()
	{
		Z_pp = 35.496563;
		B = 0.30763124;
		s0 = 29.203742;
		Y_1_pp = 42.592712;
		Y_2_pp = 33.363267;
		eta_1 = 0.46004023;
		eta_2 = 0.54541822;
	}

	Model_RRPnfL2u_21()
	{
		name = "Model_RRPnfL2u_21";
		label = "RR P_{nf} L2_u\\ (21)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(7);	
		par_unc(0) = 0.47062778;
		par_unc(1) = 0.0098010848;
		par_unc(2) = 5.3760124;
		par_unc(3) = 1.3539773;
		par_unc(4) = 1.0388395;
		par_unc(5) = 0.016521274;
		par_unc(6) = 0.0067507812;

		double corr_data[] = {
			100, 74.1, 94.1, -11.4, 2.54, 86, 4.3,
			74.1, 100, 91.4, -56.6, 5.03, 39.1, 4.99,
			94.1, 91.4, 100, -35.9, 4.84, 66.8, 5.84,
			-11.4, -51.6, -35.9, 100, 25.1, 40.1, 25.8,
			2.54, 5.03, 4.84, 25.1, 100, 11.1, 97.5,
			86, 39.1, 66.8, 40.1, 11.1, 100, 13.8,
			4.3, 4.99, 5.84, 25.8, 97.5, 13.8, 100
		};
		par_unc_corr.ResizeTo(par_unc.GetNrows(), par_unc.GetNrows());
		par_unc_corr.SetMatrixArray(corr_data);

		PrepareParameterErrorGeneratorMatrix();
	}

	virtual bool ApplyParameterChange(const TVectorD &de) override
	{
		Z_pp += de(0);
		B += de(1);
		s0 += de(2);
		Y_1_pp += de(3);
		Y_2_pp += de(4);
		eta_1 += de(5);
		eta_2 += de(6);

		return true;
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
