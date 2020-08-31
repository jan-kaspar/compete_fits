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
 
		par_unc(0) = 0.024365587;	// B
		par_unc(1) = 30.21895;	// s0
		par_unc(2) = 2.8112491;	// Y_1_pp
		par_unc(3) = 0.9540322;	// Y_2_pp
		par_unc(4) = 0.0079921224;	// eta_1
		par_unc(5) = 0.0063104218;	// eta_2 

		double corr_data[] = {
			100, 99.3, 99.6, -20.3, -97.2, -19.6,
			99.3, 100, 98.6, -21.9, -99.2, -20.9,
			99.6, 98.6, 100, -14.4, -95.7, -14.4,
			-20.3, -21.9, -14.4, 100, 26.6, 97.4,
			-27.2, -99.2, -95.7, 26.6, 100, 25.1,
			-19.6, -20.9, -14.4, 97.4, 25.1, 100
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

		// TODO: published s0 uncertainty too large
		return (fabs(de(1)) < 2.6 * par_unc(1));
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
