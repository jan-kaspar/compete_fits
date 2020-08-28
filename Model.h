#ifndef _Model_h_
#define _Model_h_

#include <string>
#include <cmath>

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TRandom2.h"

class Model
{
	public:

	std::string name;
	std::string label;
	double ranking;

	// parameter uncertainties
	TVectorD par_unc;

	// correlation matrix for parameter uncertainties
	TMatrixDSym par_unc_corr;

	// generator matrix for parameter errors
	TMatrixD par_unc_gen;

	virtual void SetDefaultParameterValues() {}

	virtual bool ApplyParameterChange(const TVectorD &de)
	{
		return false;
	}

	virtual ~Model() {}

	virtual double si_p_p(double s) const = 0;
	virtual double si_p_ap(double s) const = 0;

	virtual double rho_p_p(double s) const = 0;
	virtual double rho_p_ap(double s) const = 0;

	virtual void PrepareParameterErrorGeneratorMatrix()
	{
		const unsigned int dim = par_unc.GetNrows();

		// build covariance matrix
		TMatrixDSym cov(dim);
		for (unsigned int i = 0; i < dim; i++)
		{
			for (unsigned int j = 0; j < dim; j++)
			{
				cov(i, j) = par_unc_corr(i, j)/100. * par_unc(i) * par_unc(j);
			}
		}

		// build generator matrix
		TMatrixDSymEigen eig_decomp(cov);
		TVectorD eig_values(eig_decomp.GetEigenValues());
		TMatrixDSym S(dim);
		for (unsigned int i = 0; i < dim; i++)
			S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;


		par_unc_gen.ResizeTo(dim, dim);
		par_unc_gen = eig_decomp.GetEigenVectors() * S;
	}

	virtual bool GenerateRandomParameterErrors()
	{
		const unsigned int dim = par_unc.GetNrows();

		TVectorD de(dim);
		for (unsigned int i = 0; i < dim; i++)
			de(i) = gRandom->Gaus();
		TVectorD de_P = par_unc_gen * de;

		return ApplyParameterChange(de_P);
	}
};

#endif
