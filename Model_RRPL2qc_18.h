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
 
		par_unc(0) = 0.89981109;	// A
		par_unc(1) = 0.090048287;	// B
		par_unc(2) = 0.0028144971;	// C
		par_unc(3) = 7.5702999;	// Y_1_pp
		par_unc(4) = 0.96873121;	// Y_2_pp
		par_unc(5) = 0.010884864;	// eta_1
		par_unc(6) = 0.0064143149;	// eta_2

		double corr_data[] = {
			100, -98.2, 89.2, -99.7, 28.1, 94.6, 26.6,
			-98.2, 100, -96, 99.2, -25.5, -87, -24.1,
			89.2, -96, 100, -91.8, 20.6, 70.5, 19.4,
			-99.7, 99.2, -91.8, 100, -24.9, -91.8, -23.7,
			28.1, -25.5, 20.6, -24.9, 100, 33.1, 97.5,
			94.6, -87, 70.5, -91.8, 33.1, 100, 31.2,
			26.6, -24.1, 19.4, -23.7, 97.5, 31.2, 100
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
