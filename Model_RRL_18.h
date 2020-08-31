#include "Model.h"

class Model_RRL_18 : public Model
{
	public:

	double B_pp;
	double s0;
	double Y_1_pp;
	double Y_2_pp;
	double eta_1;
	double eta_2; 

	virtual void SetDefaultParameterValues()
	{
		B_pp = 6.6282714;
		s0 = 79.371146;
		Y_1_pp = 104.68182;
		Y_2_pp = 33.233187;
		eta_1 = 0.21134153;
		eta_2 = 0.54379057;
	}

	Model_RRL_18()
	{
		name = "Model_RRL_18";
		label = "RRL\\ (18)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0.2193893;	// B_pp
		par_unc(1) = 30.364068;	// s0
		par_unc(2) = 2.8145916;	// Y_1_pp
		par_unc(3) = 0.955194;	// Y_2_pp
		par_unc(4) = 0.0079899165;	// eta_1
		par_unc(5) = 0.0063134496;	// eta_2 

		double corr_data[] = {
			100, 99.3, 99.6, -20.2, -97.2, -19.5,
			99.3, 100, 98.6, -21.8, -99.2, -20.8,
			99.6, 98.6, 100, -14.4, -95.7, -14.4,
			-20.2, -21.8, -14.4, 100, 26.6, 97.4,
			-97.2, -99.2, -95.7, 26.6, 100, 25,
			-19.5, -20.8, -14.4, 97.4, 25, 100
		};
		par_unc_corr.ResizeTo(par_unc.GetNrows(), par_unc.GetNrows());
		par_unc_corr.SetMatrixArray(corr_data);

		PrepareParameterErrorGeneratorMatrix();
	}

	virtual bool ApplyParameterChange(const TVectorD &de) override
	{
		B_pp   += de(0);
		s0     += de(1);
		Y_1_pp += de(2);
		Y_2_pp += de(3);
		eta_1  += de(4);
		eta_2  += de(5);

		// TODO: published s0 uncertainty too large
		return (fabs(de(1)) < 2.6 * par_unc(1));
	}

	double si_p_p(double s) const override
	{
		return B_pp * log(s/s0) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return B_pp * log(s/s0) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
