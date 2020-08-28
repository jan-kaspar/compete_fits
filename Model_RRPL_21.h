#include "Model.h"

class Model_RRPL_21 : public Model
{
	public:

	double Z_pp;
	double B_pp;
	double Y_1_pp;
	double Y_2_pp;
	double eta_1;
	double eta_2;

	virtual void SetDefaultParameterValues()
	{
		Z_pp = -28.963728;
		B_pp = 6.5817218;
		Y_1_pp = 103.61068;
		Y_2_pp = 32.036043;
		eta_1 = 0.20804983;
		eta_2 = 0.53637805;
	}

	Model_RRPL_21()
	{
		name = "Model_RRPL_21";
		label = "RR PL\\ (21)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0;	// Z_pp
		par_unc(1) = 0;	// B_pp
		par_unc(2) = 0;	// Y_1_pp
		par_unc(3) = 0;	// Y_2_pp
		par_unc(4) = 0;	// eta_1
		par_unc(5) = 0;	// eta_2

		double corr_data[] = {
			1, 0, 0, 0, 0, 0,
			0, 1, 0, 0, 0, 0,
			0, 0, 1, 0, 0, 0,
			0, 0, 0, 1, 0, 0,
			0, 0, 0, 0, 1, 0,
			0, 0, 0, 0, 0, 1
		};
		par_unc_corr.ResizeTo(par_unc.GetNrows(), par_unc.GetNrows());
		par_unc_corr.SetMatrixArray(corr_data);

		PrepareParameterErrorGeneratorMatrix();
	}

	virtual bool ApplyParameterChange(const TVectorD &de) override
	{
		Z_pp   += de(0);
		B_pp   += de(1);
		Y_1_pp += de(2);
		Y_2_pp += de(3);
		eta_1  += de(4);
		eta_2  += de(5);

		return true;
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B_pp * log(s) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B_pp * log(s) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI*B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI*B_pp/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
