#include "Model.h"

class Model_RqcRcLqc_12 : public Model
{
	public:

	double B;
	double s0;
	double Y_1_pp;
	double Y_2_pip;
	double eta_1;
	double eta_2;

	virtual void SetDefaultParameterValues()
	{
		B = 0.75953082;
		s0 = 118.96503;
		Y_1_pp = 11.904844;
		Y_2_pip = 7.0913207;
		eta_1 = 0.20200261;
		eta_2 = 0.55546195;
	}

	Model_RqcRcLqc_12()
	{
		name = "Model_RqcRcLqc_12";
		label = "R^{qc} R_{c} L^{qc}\\ (12)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0;	// B
		par_unc(1) = 0;	// s0
		par_unc(2) = 0;	// Y_1_pp
		par_unc(3) = 0;	// Y_2_pip
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
		B       += de(0);
		s0      += de(1);
		Y_1_pp  += de(2);
		Y_2_pip += de(3);
		eta_1   += de(4);
		eta_2   += de(5);

		return true;
	}

	double si_p_p(double s) const override
	{
		return 9. * B * log(s/s0) + 9. * Y_1_pp * pow(s, -eta_1) - 5. * Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9. * B * log(s/s0) + 9. * Y_1_pp * pow(s, -eta_1) + 5. * Y_2_pip * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9. * M_PI * B / 2. - 9. * Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - 5. * Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9. * M_PI * B / 2. - 9. * Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + 5. * Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
