#include "Model.h"

class Model_RRL2_18 : public Model
{
	public:

	double B_pp;
	double s0;
	double Y_1_pp;
	double Y_2_pp;
	double eta_1;
	double eta_2;

	virtual void SetDefaultParameterValues()
	{
		s0 = 0.00044828241;
		eta_1 = 0.27258771;
		eta_2 = 0.55536716;
		B_pp = 0.14468157;
		Y_1_pp = 67.016128;
		Y_2_pp = 35.55018;
	}

	Model_RRL2_18()
	{
		name = "Model_RRL2_18";
		label = "RR L2\\ (18)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0.00024460879;	// s0
		par_unc(1) = 0.011454336;	// eta_1
		par_unc(2) = 0.011704097;	// eta_2
		par_unc(3) = 0.0076301918;	// B_pp
		par_unc(4) = 1.1102999;	// Y_1_pp
		par_unc(5) = 2.1720112;	// Y_2_pp

		// TODO: published uncertainty too large
		//par_unc(0) = 0.;

		double corr_data[] = {
			100., -93.8, -11.1, 99.6, -48., -9.64,
			-93.8, 100., 23.2, -90.8, 74.6, 21.8,
			-11.1, 23.2, 100., -9.91, 46.7, 98.5,
			99.6, -90.8, -9.91, 100., -41.8, -8.45,
			-48., 74.6, 46.7, -41.8, 100., 46.6,
			-9.64, 21.8, 98.5, -8.45, 46.6, 100.
		};
		par_unc_corr.ResizeTo(par_unc.GetNrows(), par_unc.GetNrows());
		par_unc_corr.SetMatrixArray(corr_data);

		PrepareParameterErrorGeneratorMatrix();
	}

	virtual bool ApplyParameterChange(const TVectorD &de) override
	{
		s0 += de(0);
		eta_1 += de(1);
		eta_2 += de(2);
		B_pp += de(3);
		Y_1_pp += de(4);
		Y_2_pp += de(5);

		return (s0 > 1E-4);
	}

	double si_p_p(double s) const override
	{
		return B_pp * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return B_pp * pow(log(s/s0), 2) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B_pp * log(s/s0) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
