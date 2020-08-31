#include "Model.h"
#include "all_models.h"

#include "stat.h"

#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void ProcessOneModel(Model *model)
{
	// settings
	unsigned int n_error_samples = (unsigned int) 1E4;

	unsigned int si_for_corr = 194;

	// define grid of s values
	vector<double> values_W;
	for (double sqrt_s = 10; sqrt_s < 2E5; sqrt_s *= 1.05)
		values_W.push_back(sqrt_s);

	printf("    si_for_corr = %u, sqrt s = %E\n", si_for_corr, values_W[si_for_corr]);

	// evaluate uncertainty
	vector<Stat> stat_si_p_p(values_W.size(), Stat(1));
	vector<Stat> stat_si_p_ap(values_W.size(), Stat(1));
	vector<Stat> stat_rho_p_p(values_W.size(), Stat(1));
	vector<Stat> stat_rho_p_ap(values_W.size(), Stat(1));

	TGraph *g_de_rho_p_p_vs_de_si_p_p = new TGraph();
	TH2D *h_de_rho_p_p_vs_de_si_p_p = new TH2D("", ";#Delta#sigma_{pp}   (mb);#Delta#rho_{pp}", 100, -20., +20., 100, -0.02, +0.02);
	
	for (unsigned int ci = 0; ci < n_error_samples; ci++)
	{
		const unsigned int dim = values_W.size();
		vector<double> v_si_p_p(dim), v_si_p_ap(dim), v_rho_p_p(dim), v_rho_p_ap(dim);

		// generate model parameter errors
		model->SetDefaultParameterValues();
		bool acceptable_configuration = model->GenerateRandomParameterErrors();

		if (!acceptable_configuration)
			continue;

		// evaluate biased model on the s grid
		for (unsigned int si = 0; si < dim; si++)
		{
			const double s = values_W[si] * values_W[si];
			v_si_p_p[si] = model->si_p_p(s);
			v_si_p_ap[si] = model->si_p_ap(s);
			v_rho_p_p[si] = model->rho_p_p(s);
			v_rho_p_ap[si] = model->rho_p_ap(s);
		}

		// subtract central values
		model->SetDefaultParameterValues();
		for (unsigned int si = 0; si < dim; si++)
		{
			const double s = values_W[si] * values_W[si];
			v_si_p_p[si] -= model->si_p_p(s);
			v_si_p_ap[si] -= model->si_p_ap(s);
			v_rho_p_p[si] -= model->rho_p_p(s);
			v_rho_p_ap[si] -= model->rho_p_ap(s);
		}

		// fill in statistics
		for (unsigned int si = 0; si < dim; si++)
		{
			stat_si_p_p[si].Fill(v_si_p_p[si]);
			stat_si_p_ap[si].Fill(v_si_p_ap[si]);
			stat_rho_p_p[si].Fill(v_rho_p_p[si]);
			stat_rho_p_ap[si].Fill(v_rho_p_ap[si]);
		}

		// fill in correlation plots
		int idx = g_de_rho_p_p_vs_de_si_p_p->GetN();
		g_de_rho_p_p_vs_de_si_p_p->SetPoint(idx, v_si_p_p[si_for_corr], v_rho_p_p[si_for_corr]);
		h_de_rho_p_p_vs_de_si_p_p->Fill(v_si_p_p[si_for_corr], v_rho_p_p[si_for_corr]);
	}

	// make graphs
	struct GraphSet
	{
		TGraph *g_cen_val, *g_unc_mean, *g_unc_stddev, *g_band_up, *g_band_dw;

		GraphSet()
		{
			g_cen_val = new TGraph();
			g_unc_mean = new TGraph();
			g_unc_stddev = new TGraph();
			g_band_up = new TGraph(); g_band_up->SetLineColor(2);
			g_band_dw = new TGraph(); g_band_dw->SetLineColor(2);
		}

		void Fill(double W, double cv, double mean, double stddev)
		{
			int idx = g_cen_val->GetN();

			g_cen_val->SetPoint(idx, W, cv);
			g_unc_mean->SetPoint(idx, W, mean);
			g_unc_stddev->SetPoint(idx, W, stddev);
			g_band_up->SetPoint(idx, W, cv + mean + stddev);
			g_band_dw->SetPoint(idx, W, cv + mean - stddev);
		}

		void Write() const
		{
			g_cen_val->Write("g_cen_val");
			g_unc_mean->Write("g_unc_mean");
			g_unc_stddev->Write("g_unc_stddev");
			g_band_up->Write("g_band_up");
			g_band_dw->Write("g_band_dw");
		}
	};

	GraphSet gs_si_p_p, gs_si_p_ap, gs_rho_p_p, gs_rho_p_ap;

	model->SetDefaultParameterValues();

	for (unsigned int si = 0; si < values_W.size(); si++)
	{
		const double W = values_W[si];
		const double s = values_W[si] * values_W[si];

		gs_si_p_p.Fill(W, model->si_p_p(s), stat_si_p_p[si].GetMean(0), stat_si_p_p[si].GetStdDev(0));
		gs_si_p_ap.Fill(W, model->si_p_ap(s), stat_si_p_ap[si].GetMean(0), stat_si_p_ap[si].GetStdDev(0));
		gs_rho_p_p.Fill(W, model->rho_p_p(s), stat_rho_p_p[si].GetMean(0), stat_rho_p_p[si].GetStdDev(0));
		gs_rho_p_ap.Fill(W, model->rho_p_ap(s), stat_rho_p_ap[si].GetMean(0), stat_rho_p_ap[si].GetStdDev(0));
	}

	// save graphs
	TDirectory *d_top = gDirectory;

	g_de_rho_p_p_vs_de_si_p_p->Write("g_de_rho_p_p_vs_de_si_p_p");
	h_de_rho_p_p_vs_de_si_p_p->Write("h_de_rho_p_p_vs_de_si_p_p");

	gDirectory = d_top->mkdir("si_p_p");
	gs_si_p_p.Write();

	gDirectory = d_top->mkdir("si_p_ap");
	gs_si_p_ap.Write();

	gDirectory = d_top->mkdir("rho_p_p");
	gs_rho_p_p.Write();

	gDirectory = d_top->mkdir("rho_p_ap");
	gs_rho_p_ap.Write();

	gDirectory = d_top;

	// print some benchmarks
	for (double W : { 2760., 7000., 8000., 13000.} )
	{
		printf("    sqrt s = %5.0f GeV: sigma_tot = (%6.2f +- %.2f) mb, rho = (%.4f +- %.4f)\n", W,
			gs_si_p_p.g_cen_val->Eval(W), gs_si_p_p.g_unc_stddev->Eval(W),
			gs_rho_p_p.g_cen_val->Eval(W), gs_rho_p_p.g_unc_stddev->Eval(W)
		);
	}
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// define models
	vector<Model *> models;
	// blue
	models.push_back(new Model_RRdPL2_20);
	models.push_back(new Model_RRdPL2u_17);
	models.push_back(new Model_RRdPL2u_19);
	models.push_back(new Model_RRdPqcL2u_16);
	models.push_back(new Model_RRcdPL2u_15);
	models.push_back(new Model_RRcdPqcL2u_14);
	models.push_back(new Model_RRPL2u_19);
	models.push_back(new Model_RRPnfL2u_21);

	// blue dashed
	models.push_back(new Model_RRPEu_19);

	// magenta
	models.push_back(new Model_RqcRcL2qc_12);
	models.push_back(new Model_RRcL2qc_15);
	models.push_back(new Model_RRL2_18);
	models.push_back(new Model_RRL2qc_17);

	// green
	models.push_back(new Model_RqcRcLqc_12);
	models.push_back(new Model_RqcRLqc_14);
	models.push_back(new Model_RRL_18);
	models.push_back(new Model_RRLnf_19);
	models.push_back(new Model_RRLqc_17);
	models.push_back(new Model_RRPL_21);
	models.push_back(new Model_RRcLqc_15);
	models.push_back(new Model_RRcPL_19);

	// green dashed
	models.push_back(new Model_RRPL2_20);
	models.push_back(new Model_RRPL2qc_18);

	// prepare output file
	TFile *f_out = TFile::Open("uncertainties.root", "recreate");

	// configure random-number engine
	gRandom->SetSeed(1);

	// evaluate model uncertainties
	for (auto &model : models)
	{
		printf("* %s\n", model->name.c_str());

		gDirectory = f_out->mkdir(model->name.c_str());

		ProcessOneModel(model);
	}

	// clean up
	delete f_out;
	return 0;
}
