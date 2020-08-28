#include "Model.h"

class Model_RRdPqcL2u_16 : public Model
{
	public:

	double Z;
	double B;
	double s0;
	double Y_1_pp;
	double Y_2_pp;
	double eta;   

	virtual void SetDefaultParameterValues()
	{
		Z = 4.1266114;
		B = 0.31640911;
		s0 = 46.272116;
		Y_1_pp = 46.412235;
		Y_2_pp = 34.157737;
		eta = 0.55062892;
	}

	Model_RRdPqcL2u_16()
	{
		name = "Model_RRdPqcL2u_16";
		label = "(RR)^d P^{qc} L2_u\\ (16)";

		SetDefaultParameterValues();

		par_unc.ResizeTo(6);
 
		par_unc(0) = 0;	// Z
		par_unc(1) = 0;	// B
		par_unc(2) = 0;	// s0
		par_unc(3) = 0;	// Y_1_pp
		par_unc(4) = 0;	// Y_2_pp
		par_unc(5) = 0;	// eta   

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
		Z      += de(0);
		B      += de(1);
		s0     += de(2);
		Y_1_pp += de(3);
		Y_2_pp += de(4);
		eta    += de(5);

		return true;
	}

	double si_p_p(double s) const override
	{
		return 9.*Z + B * pow(log(s/s0), 2) + (Y_1_pp - Y_2_pp) * pow(s, -eta);
	}

	double si_p_ap(double s) const override
	{
		return 9.*Z + B * pow(log(s/s0), 2) + (Y_1_pp + Y_2_pp) * pow(s, -eta);
	}

	double rho_p_p(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) - Y_2_pp * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_p(s);
	}

	double rho_p_ap(double s) const override
	{
		double rho_si = M_PI * B * log(s/s0) - Y_1_pp * pow(s, -eta) / tan((1.-eta)/2.*M_PI) + Y_2_pp * pow(s, -eta) * tan((1.-eta)/2.*M_PI);
		return rho_si / si_p_ap(s);
	}
};
