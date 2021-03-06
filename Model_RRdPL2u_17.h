#include "Model.h"

class Model_RRdPL2u_17 : public Model
{
	public:

	double Z_pp;
	double B_pp;
	double s0;
	double Y_1_pp;
	double Y_2_pp;
	double eta;   

	virtual void SetDefaultParameterValues()
	{
		Z_pp = 36.863204;
		B_pp = 0.32009024;
		s0 = 43.88116;
		Y_1_pp = 45.42289;
		Y_2_pp = 30.53598;
		eta = 0.52949847;
	}

	Model_RRdPL2u_17()
	{
		name = "Model_RRdPL2u_17";
		label = "(RR)^d PL2_u\\ (17)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0;	// Z_pp
		par_unc(1) = 0;	// B_pp
		par_unc(2) = 0;	// s0
		par_unc(3) = 0;	// Y_1_pp
		par_unc(4) = 0;	// Y_2_pp
		par_unc(5) = 0;	// eta   

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
		s0     += de(2);
		Y_1_pp += de(3);
		Y_2_pp += de(4);
		eta    += de(5);

		return true;
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B_pp * pow(log(s/s0), 2) + (Y_1_pp - Y_2_pp) * pow(s, -eta);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B_pp * pow(log(s/s0), 2) + (Y_1_pp + Y_2_pp) * pow(s, -eta);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) - Y_2_pp * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) + Y_2_pp * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
