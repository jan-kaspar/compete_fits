#include "Model.h"

class Model_RRPEu_19 : public Model
{
	public:

	double Z_pp;
	double X;
	double ep;
	double Y_1_pp;
	double Y_2_pp;
	double eta_1;
	double eta_2;

	virtual void SetDefaultParameterValues()
	{
		Z_pp = 4.1824033;
		X = 15.685953;
		ep = 0.10134777;
		Y_1_pp = 57.450783;
		Y_2_pp = 33.468752;
		eta_1 = 0.34309186;
		eta_2 = 0.54487001;
	}

	Model_RRPEu_19()
	{
		name = "Model_RRPEu_19";
		label = "RR PE_u\\ (19)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(7);
 
		par_unc(0) = 0;	// Z_pp
		par_unc(1) = 0;	// X
		par_unc(2) = 0;	// ep
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
		Z_pp   += de(0);
		X      += de(1);
		ep     += de(2);
		Y_1_pp += de(3);
		Y_2_pp += de(4);
		eta_1  += de(5);
		eta_2  += de(6);

		return true;
	}

	double si_p_p(double s) const override
	{
		return Z_pp + X * pow(s, ep) + Y_1_pp * pow(s, -eta_1) - Y_2_pp * pow(s, -eta_2);
	}

	double si_p_ap(double s) const override
	{
		return Z_pp + X * pow(s, ep) + Y_1_pp * pow(s, -eta_1) + Y_2_pp * pow(s, -eta_2);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = - X*pow(s, ep) / tan((1.+ep)/2*M_PI) - Y_1_pp*pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) - Y_2_pp*pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = - X*pow(s, ep) / tan((1.+ep)/2*M_PI) - Y_1_pp*pow(s, -eta_1) / tan((1.-eta_1)/2*M_PI) + Y_2_pp*pow(s, -eta_2) * tan((1.-eta_2)/2*M_PI);
		return rho_si / si_p_ap(s);
	}
};
