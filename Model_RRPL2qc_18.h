#include "Model.h"

class Model_RRPL2qc_18 : public Model
{
	public:

	double A;
	double B;
	double C;
	double Y_1_pp;
	double Y_2_pp;
	double eta_1;
	double eta_2;

	virtual void SetDefaultParameterValues()
	{
		A = -3.7258004;
		B = 0.79164674;
		C = -0.0018160896;
		Y_1_pp = 109.05955;
		Y_2_pp = 33.077245;
		eta_1 = 0.20641782;
		eta_2 = 0.5428034;
	}

	Model_RRPL2qc_18()
	{
		name = "Model_RRPL2qc_18";
		label = "RR (PL2)^{qc}\\ (18)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(7);
 
		par_unc(0) = 0;	// A
		par_unc(1) = 0;	// B
		par_unc(2) = 0;	// C
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
		A      += de(0);
		B      += de(1);
		C      += de(2);
		Y_1_pp += de(3);
		Y_2_pp += de(4);
		eta_1  += de(5);
		eta_2  += de(6);

		return true;
	}

	double si_p_p(double s) const override
	{
		return 9.*(A + B*log(s) + C*pow(log(s), 2)) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return 9.*(A + B*log(s) + C*pow(log(s), 2)) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = 9.*M_PI * (B/2. + C*log(s)) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) - Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_p(s);
		return 0;
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = 9.*M_PI * (B/2. + C*log(s)) - Y_1_pp * pow(s, -eta_1) / tan((1.-eta_1)/2.*M_PI) + Y_2_pp * pow(s, -eta_2) * tan((1.-eta_2)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
