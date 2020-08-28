#include "Model.h"

class Model_RRPL2u_19 : public Model
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
		Z_pp = 35.865711;
		B = 0.31573;
		s0 = 34.409806;
		Y_1_pp = 42.069653;
		Y_2_pp = 32.155544;
		eta_1 = 0.46822531;
		eta_2 = 0.53956628;
	}

	Model_RRPL2u_19()
	{
		name = "Model_RRPL2u_19";
		label = "RR PL2_u\\ (19)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(7);
 
		par_unc(0) = 0;	// Z_pp
		par_unc(1) = 0;	// B
		par_unc(2) = 0;	// s0
		par_unc(3) = 0;	// Y_1_pp
		par_unc(4) = 0;	// Y_2_pp
		par_unc(5) = 0;	// eta_1
		par_unc(6) = 0;	// eta_2

		double corr_data[] = {
			1, 0, 0, 0, 0, 0, 0,
			0, 1, 0, 0, 0, 0, 0,
			0, 0, 1, 0, 0, 0, 0,
			0, 0, 0, 1, 0, 0, 0,
			0, 0, 0, 0, 1, 0, 0,
			0, 0, 0, 0, 0, 1, 0,
			0, 0, 0, 0, 0, 0, 1
		};
		par_unc_corr.ResizeTo(par_unc.GetNrows(), par_unc.GetNrows());
		par_unc_corr.SetMatrixArray(corr_data);

		PrepareParameterErrorGeneratorMatrix();
	}

	virtual bool ApplyParameterChange(const TVectorD &de) override
	{
		Z_pp   += de(0);
		B      += de(1);
		s0     += de(2);
		Y_1_pp += de(3);
		Y_2_pp += de(4);
		eta_1  += de(5);
		eta_2  += de(6);

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
