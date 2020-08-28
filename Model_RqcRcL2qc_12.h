#include "Model.h"

class Model_RqcRcL2qc_12 : public Model
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
		B = 0.016080647;
		s0 = 0.00045002716;
		Y_1_pp = 7.4402016;
		Y_2_pip = 7.1370767;
		eta_1 = 0.2722584;
		eta_2 = 0.55494024;
	}

	Model_RqcRcL2qc_12()
	{
		name = "Model_RqcRcL2qc_12";
		label = "R^{qc} R_{c} L2^{qc}\\ (12)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0.011413449;	// eta_1
		par_unc(1) = 0.011591606;	// eta_2
		par_unc(2) = 0.0008444799;	// B
		par_unc(3) = 0.00024455049;	// s0
		par_unc(4) = 0.12274107;	// Y_1_pp
		par_unc(5) = 0.42528539;	// Y_2_pip

		// TODO: published uncertainty too large
		//par_unc(3) = 0.;

		double corr_data[] = {
			100., 24., -90.8, -93.8, 75., 23.4,
			24., 100., -10.9, -12.1, 46.6, 99.,
			-90.8, -10.9, 100., 99.6, -42.2, -10.2,
			-93.8, -12.1, 99.6, 100., -48.4, -11.4,
			75., 46.6, -42.2, -48.4, 100., 46.7,
			23.4, 99., -10.2, -11.4, 46.7, 100.
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
		Y_2_pip += de(5);

		return (s0 > 1E-4);
	}

	double si_p_p(double s) const override
	{
		return 9. * B * pow(log(s/s0), 2.) + 9. * Y_1_pp * pow(s, -eta_1) - 5. * Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9. * B * pow(log(s/s0), 2.) + 9. * Y_1_pp * pow(s, -eta_1) + 5. * Y_2_pip * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9. * M_PI * B * log(s/s0) - 9. * Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - 5. * Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9. * M_PI * B * log(s/s0) - 9. * Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + 5. * Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
