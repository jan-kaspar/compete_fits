#include "Model.h"

class Model_RRcdPL2u_15 : public Model
{
	public:

	double Z_pp;
	double B;
	double s0;
	double Y_1_pp;
	double Y_2_pip;
	double eta;    

	virtual void SetDefaultParameterValues()
	{
		Z_pp = 36.784336;
		B = 0.31175883;
		s0 = 39.731955;
		Y_1_pp = 48.60371;
		Y_2_pip = 6.5969507;
		eta = 0.54348946;
	}

	Model_RRcdPL2u_15()
	{
		name = "Model_RRcdPL2u_15";
		label = "(RR_{c})^{d} PL2_{u}\\ (15)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0;	// Z_pp
		par_unc(1) = 0;	// B
		par_unc(2) = 0;	// s0
		par_unc(3) = 0;	// Y_1_pp
		par_unc(4) = 0;	// Y_2_pip
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
		Z_pp    += de(0);
		B       += de(1);
		s0      += de(2);
		Y_1_pp  += de(3);
		Y_2_pip += de(4);
		eta     += de(5);

		return true;
	}

	double si_p_p(double s) const override
	{
		return Z_pp + B * pow(log(s/s0), 2) + (Y_1_pp - 5.*Y_2_pip) * pow(s, -eta);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + B * pow(log(s/s0), 2) + (Y_1_pp + 5.*Y_2_pip) * pow(s, -eta);
		return 0;
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) - 5.*Y_2_pip * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) + 5.*Y_2_pip * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_ap(s);
		return 0;
	}
};
