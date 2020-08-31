#include "Model.h"

class Model_RRL2qc_17 : public Model
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
		eta_1 = 0.27261187;
		eta_2 = 0.5552468;
		B = 0.016076625;
		s0 = 0.00044830086;
		Y_1_pp = 67.01912;
		Y_2_pp = 35.528122;
	}

	Model_RRL2qc_17()
	{
		name = "Model_RRL2qc_17";
		label = "RR L2^{qc}\\ (17)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
		par_unc(0) = 0.01145315;	// eta_1
		par_unc(1) = 0.011697771;	// eta_2
		par_unc(2) = 0.00084766039;	// B
		par_unc(3) = 0.00024455797;	// s0
		par_unc(4) = 1.110379;		// Y_1_pp
		par_unc(5) = 2.1696148;		// Y_2_pp

		double corr_data[] = {
			100., 23.2, -90.8, -93.8, 74.6, 21.8,
			23.2, 100., -9.92, -11.1, 46.7, 98.5,
			-90.8, -9.92, 100., 99.6, -41.8, -8.46,
			-93.8, -11.1, 99.6, 100., -48., -9.65,
			74.6, 46.7, -41.8, -48., 100., 46.6,
			21.8, 98.5, -8.46, -9.65, 46.6, 100.
		};
		par_unc_corr.ResizeTo(par_unc.GetNrows(), par_unc.GetNrows());
		par_unc_corr.SetMatrixArray(corr_data);

		PrepareParameterErrorGeneratorMatrix();
	}

	virtual bool ApplyParameterChange(const TVectorD &de) override
	{
		eta_1 += de(0);
		eta_2 += de(1);
		B += de(2);
		s0 += de(3);
		Y_1_pp += de(4);
		Y_2_pp += de(5);

		// TODO: published s0 uncertainty too large
		return (fabs(de(3)) < 1.8 * par_unc(3));
	}

	double si_p_p(double s) const override
	{
		return 9.*B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI*B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI*B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
