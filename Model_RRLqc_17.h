#include "Model.h"

class Model_RRLqc_17 : public Model
{
	public:

	double B;
	double s0;
	double Y_1_pp;
	double Y_2_pp;
	double eta_1;
	double eta_2;

	virtual void SetDefaultParameterValues()
	{
		B = 0.73635067;
		s0 = 79.027935;
		Y_1_pp = 104.66237;
		Y_2_pp = 33.205188;
		eta_1 = 0.21146055;
		eta_2 = 0.54360032;
	}

	Model_RRLqc_17()
	{
		name = "Model_RRLqc_17";
		label = "RR L^{qc}\\ (17)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0;	// B
		par_unc(1) = 0;	// s0
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
		B      += de(0);
		s0     += de(1);
		Y_1_pp += de(2);
		Y_2_pp += de(3);
		eta_1  += de(4);
		eta_2  += de(5);

		return true;
	}

	double si_p_p(double s) const override
	{
		return 9.*B*log(s/s0) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*B*log(s/s0) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI*B/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI*B/2. - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
