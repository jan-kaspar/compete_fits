#include "Model.h"
#include "all_models.h"

#include "stat.h"

#include <vector>

#include "TFile.h"
#include "TGraph.h"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TRandom2.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// model with central values
	Model_RRPL2u_21 *model = new Model_RRPL2u_21;

	// prepare output file
	TFile *f_out = TFile::Open("uncertainty_RRPL2u_21.root", "recreate");

	// define grid of s values
	vector<double> values_W;
	for (double sqrt_s = 10; sqrt_s < 2E5; sqrt_s *= 1.05)
		values_W.push_back(sqrt_s);

	// order of quantities: Z_pp, B, s0, Y_1_pp, Y_2_pp, eta_1, eta_2

	// sigmas
	TVectorD unc(7);
	unc(0) = 0.47062778;
	unc(1) = 0.0098010848;
	unc(2) = 5.3760124;
	unc(3) = 1.3539773;
	unc(4) = 1.0388395;
	unc(5) = 0.016521274;
	unc(6) = 0.0067507812;

	// correlation matrix from PDF
	TMatrixDSym corr(7);
	double corr_data[] = {
		100, 74.1, 94.1, -11.4, 2.54, 86, 4.3,
		74.1, 100, 91.4, -56.6, 5.03, 39.1, 4.99,
		94.1, 91.4, 100, -35.9, 4.84, 66.8, 5.84,
		-11.4, -51.6, -35.9, 100, 25.1, 40.1, 25.8,
		2.54, 5.03, 4.84, 25.1, 100, 11.1, 97.5,
		86, 39.1, 66.8, 40.1, 11.1, 100, 13.8,
		4.3, 4.99, 5.84, 25.8, 97.5, 13.8, 100
	};
	/*
	double corr_data[] = {
		0, 0, 0, 0, 0, 0, 0,
		0, 100, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0,
	};
	*/
	corr.SetMatrixArray(corr_data);

	// build covariance matrix
	TMatrixDSym cov(7);
	for (unsigned int i = 0; i < 7; i++)
	{
		for (unsigned int j = 0; j < 7; j++)
		{
			cov(i, j) = corr(i, j)/100. * unc(i) * unc(j);
		}
	}

	// build generator matrix
	TMatrixDSymEigen eig_decomp(cov);
	TVectorD eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(7);
	for (unsigned int i = 0; i < 7; i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	auto gen_mat = eig_decomp.GetEigenVectors() * S;

	// evaluate uncertainty
	vector<Stat> stat_si_p_p(values_W.size(), Stat(1));
	vector<Stat> stat_si_p_ap(values_W.size(), Stat(1));
	vector<Stat> stat_rho_p_p(values_W.size(), Stat(1));
	vector<Stat> stat_rho_p_ap(values_W.size(), Stat(1));
	
	for (unsigned int ci = 0; ci < 10000; ci++)
	{
		// generate model parameter errors
		TVectorD de(7);
		for (unsigned int i = 0; i < 7; i++)
			de(i) = gRandom->Gaus();
		TVectorD de_P = gen_mat * de;

		// get model with biased parameters;
		Model_RRPL2u_21 *model_bias = new Model_RRPL2u_21;
		model_bias->Z_pp += de_P(0);
		model_bias->B += de_P(1);
		model_bias->s0 += de_P(2);
		model_bias->Y_1_pp += de_P(3);
		model_bias->Y_2_pp += de_P(4);
		model_bias->eta_1 += de_P(5);
		model_bias->eta_2 += de_P(6);

		for (unsigned int si = 0; si < values_W.size(); si++)
		{
			const double s = values_W[si] * values_W[si];
			const double de_si_p_p = model_bias->si_p_p(s) - model->si_p_p(s);
			const double de_si_p_ap = model_bias->si_p_ap(s) - model->si_p_ap(s);
			const double de_rho_p_p = model_bias->rho_p_p(s) - model->rho_p_p(s);
			const double de_rho_p_ap = model_bias->rho_p_ap(s) - model->rho_p_ap(s);

			stat_si_p_p[si].Fill(de_si_p_p);
			stat_si_p_ap[si].Fill(de_si_p_ap);
			stat_rho_p_p[si].Fill(de_rho_p_p);
			stat_rho_p_ap[si].Fill(de_rho_p_ap);
		}

		delete model_bias;
	}

	// make graphs
	TGraph *g_mean_si_p_p = new TGraph();
	TGraph *g_stddev_si_p_p = new TGraph();
	TGraph *g_band_up_si_p_p = new TGraph();
	TGraph *g_band_dw_si_p_p = new TGraph();

	TGraph *g_mean_si_p_ap = new TGraph();
	TGraph *g_stddev_si_p_ap = new TGraph();
	TGraph *g_band_up_si_p_ap = new TGraph();
	TGraph *g_band_dw_si_p_ap = new TGraph();

	TGraph *g_mean_rho_p_p = new TGraph();
	TGraph *g_stddev_rho_p_p = new TGraph();
	TGraph *g_band_up_rho_p_p = new TGraph();
	TGraph *g_band_dw_rho_p_p = new TGraph();

	TGraph *g_mean_rho_p_ap = new TGraph();
	TGraph *g_stddev_rho_p_ap = new TGraph();
	TGraph *g_band_up_rho_p_ap = new TGraph();
	TGraph *g_band_dw_rho_p_ap = new TGraph();

	for (unsigned int si = 0; si < values_W.size(); si++)
	{
		const double W = values_W[si];
		const double s = values_W[si] * values_W[si];

		const double si_p_p = model->si_p_p(s);
		const double si_p_ap = model->si_p_ap(s);
		const double rho_p_p = model->rho_p_p(s);
		const double rho_p_ap = model->rho_p_ap(s);

		int idx = g_mean_si_p_p->GetN();

		g_mean_si_p_p->SetPoint(idx, W, stat_si_p_p[si].GetMean(0));
		g_stddev_si_p_p->SetPoint(idx, W, stat_si_p_p[si].GetStdDev(0));
		g_band_up_si_p_p->SetPoint(idx, W, si_p_p + stat_si_p_p[si].GetMean(0) + stat_si_p_p[si].GetStdDev(0));
		g_band_dw_si_p_p->SetPoint(idx, W, si_p_p + stat_si_p_p[si].GetMean(0) - stat_si_p_p[si].GetStdDev(0));

		g_mean_si_p_ap->SetPoint(idx, W, stat_si_p_ap[si].GetMean(0));
		g_stddev_si_p_ap->SetPoint(idx, W, stat_si_p_ap[si].GetStdDev(0));
		g_band_up_si_p_ap->SetPoint(idx, W, si_p_ap + stat_si_p_ap[si].GetMean(0) + stat_si_p_ap[si].GetStdDev(0));
		g_band_dw_si_p_ap->SetPoint(idx, W, si_p_ap + stat_si_p_ap[si].GetMean(0) - stat_si_p_ap[si].GetStdDev(0));

		g_mean_rho_p_p->SetPoint(idx, W, stat_rho_p_p[si].GetMean(0));
		g_stddev_rho_p_p->SetPoint(idx, W, stat_rho_p_p[si].GetStdDev(0));
		g_band_up_rho_p_p->SetPoint(idx, W, rho_p_p + stat_rho_p_p[si].GetMean(0) + stat_rho_p_p[si].GetStdDev(0));
		g_band_dw_rho_p_p->SetPoint(idx, W, rho_p_p + stat_rho_p_p[si].GetMean(0) - stat_rho_p_p[si].GetStdDev(0));

		g_mean_rho_p_ap->SetPoint(idx, W, stat_rho_p_ap[si].GetMean(0));
		g_stddev_rho_p_ap->SetPoint(idx, W, stat_rho_p_ap[si].GetStdDev(0));
		g_band_up_rho_p_ap->SetPoint(idx, W, rho_p_ap + stat_rho_p_ap[si].GetMean(0) + stat_rho_p_ap[si].GetStdDev(0));
		g_band_dw_rho_p_ap->SetPoint(idx, W, rho_p_ap + stat_rho_p_ap[si].GetMean(0) - stat_rho_p_ap[si].GetStdDev(0));
	}

	g_mean_si_p_p->Write("g_mean_si_p_p");
	g_stddev_si_p_p->Write("g_stddev_si_p_p");
	g_band_up_si_p_p->Write("g_band_up_si_p_p");
	g_band_dw_si_p_p->Write("g_band_dw_si_p_p");

	g_mean_si_p_ap->Write("g_mean_si_p_ap");
	g_stddev_si_p_ap->Write("g_stddev_si_p_ap");
	g_band_up_si_p_ap->Write("g_band_up_si_p_ap");
	g_band_dw_si_p_ap->Write("g_band_dw_si_p_ap");

	g_mean_rho_p_p->Write("g_mean_rho_p_p");
	g_stddev_rho_p_p->Write("g_stddev_rho_p_p");
	g_band_up_rho_p_p->Write("g_band_up_rho_p_p");
	g_band_dw_rho_p_p->Write("g_band_dw_rho_p_p");

	g_mean_rho_p_ap->Write("g_mean_rho_p_ap");
	g_stddev_rho_p_ap->Write("g_stddev_rho_p_ap");
	g_band_up_rho_p_ap->Write("g_band_up_rho_p_ap");
	g_band_dw_rho_p_ap->Write("g_band_dw_rho_p_ap");

	// cross-check
	TGraph *g_stddev_check_si_p_p = new TGraph();
	TGraph *g_stddev_check_rho_p_p = new TGraph();
	for (unsigned int si = 0; si < values_W.size(); si++)
	{
		const double W = values_W[si];
		const double s = values_W[si] * values_W[si];

		double L = log(s/model->s0);
		double den = model->Z_pp + model->B*L*L;
		double unc_si_p_p = L*L * unc(1);
		double unc_rho_p_p = M_PI*model->B*L/den * (1./model->B - L*L/den) * unc(1);

		int idx = g_stddev_check_si_p_p->GetN();
		g_stddev_check_si_p_p->SetPoint(idx, W, unc_si_p_p);
		g_stddev_check_rho_p_p->SetPoint(idx, W, unc_rho_p_p);
	}

	g_stddev_check_si_p_p->SetLineColor(2);
	g_stddev_check_si_p_p->Write("g_stddev_check_si_p_p");

	g_stddev_check_rho_p_p->SetLineColor(2);
	g_stddev_check_rho_p_p->Write("g_stddev_check_rho_p_p");

	// clean up
	delete f_out;
	return 0;
}
