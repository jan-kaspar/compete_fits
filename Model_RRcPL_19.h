#include "Model.h"

class Model_RRcPL_19 : public Model
{
	public:

	double Z_pp;
	double B_pp;
	double Y_1_pp;
	double Y_2_pip;
	double eta_1;
	double eta_2;  

	virtual void SetDefaultParameterValues()
	{
		Z_pp = -27.389097;
		B_pp = 6.4776681;
		Y_1_pp = 101.95344;
		Y_2_pip = 6.1769888;
		eta_1 = 0.21077046;
		eta_2 = 0.53045536;
	}

	Model_RRcPL_19()
	{
		name = "Model_RRcPL_19";
		label = "RR_c PL\\ (19)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 3.8940536;	// Z_pp
		par_unc(1) = 0.24508299;	// B_pp
		par_unc(2) = 3.1559926;	// Y_1_pp
		par_unc(3) = 0.19714977;	// Y_2_pip
		par_unc(4) = 0.0091272007;	// eta_1
		par_unc(5) = 0.0069818493;	// eta_2  

		double corr_data[] = {
			100, -99.6, -98.7, 14.6, 98.3, 14.1,
			-99.6, 100, 99.6, -12.7, -96.3, -12.4,
			-98.7, 99.6, 100, -6.94, -94.2, -6.97,
			14.6, -12.7, -6.94, 100, 20.8, 98.5,
			98.3, -96.3, -94.2, 20.8, 100, 20,
			14.1, -12.4, -6.97, 98.5, 20, 100
		};
		par_unc_corr.ResizeTo(par_unc.GetNrows(), par_unc.GetNrows());
		par_unc_corr.SetMatrixArray(corr_data);

		PrepareParameterErrorGeneratorMatrix();
	}

	virtual bool ApplyParameterChange(const TVectorD &de) override
	{
		Z_pp    += de(0);
		B_pp    += de(1);
		Y_1_pp  += de(2);
		Y_2_pip += de(3);
		eta_1   += de(4);
		eta_2   += de(5);

		return true;
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B_pp * log(s) + Y_1_pp * pow(s, -eta_1) - 5.*Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B_pp * log(s) + Y_1_pp * pow(s, -eta_1) + 5.*Y_2_pip * pow(s, -eta_2);
		return 0;
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
