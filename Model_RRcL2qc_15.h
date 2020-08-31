#include "Model.h"

class Model_RRcL2qc_15 : public Model
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
		eta_1 = 0.2714367;
		eta_2 = 0.55573821;
		B = 0.016178984;
		s0 = 0.00047862371;
		Y_1_pp = 67.017688;
		Y_2_pip = 7.1656375;
	}

	Model_RRcL2qc_15()
	{
		name = "Model_RRcL2qc_15";
		label = "RR_{\\rm c} L2^{\\rm qc}\\ (15)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
		par_unc(0) = 0.011373472;	// eta_1
		par_unc(1) = 0.01162258;	// eta_2
		par_unc(2) = 0.00084914046;	// B
		par_unc(3) = 0.0002593784;	// s0
		par_unc(4) = 1.0998903;		// Y_1_pp
		par_unc(5) = 0.38665906;	// Y_2_pip

		double corr_data[] = {
			100., 23.7, -90.8, -93.8, 74.7, 23.1,
			23.7, 100., -10.5, -11.7, 46.6, 99.,
			-90.8, -10.5, 100., 99.6, -41.8, -9.85,
			-93.8, -11.7, 99.6, 100., -48.1, -11.1,
			74.7, 46.6, -41.8, -48.1, 100., 46.7,
			23.1, 99., -9.85, -11.1, 46.7, 100.
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

		// TODO: published s0 uncertainty too large
		return (fabs(de(3)) < 1.8 * par_unc(3));
	}

	double si_p_p(double s) const override
	{
		return 9.*B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) - 5.*Y_2_pip * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*B * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) + 5.*Y_2_pip * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + 5.*Y_2_pip * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
