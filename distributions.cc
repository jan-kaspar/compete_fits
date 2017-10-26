#include "Model.h"

#include "all_models.h"

#include <vector>

#include "TFile.h"
#include "TGraph.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// define models
	vector<Model *> models;
	models.push_back(new Model_RPdPL2_20);
	models.push_back(new Model_RPdPL2u_17);
	models.push_back(new Model_RPdPL2u_19);
	models.push_back(new Model_RPdPqcL2u_16);
	models.push_back(new Model_RqcRcL2qc_12);
	models.push_back(new Model_RqcRcLqc_12);
	models.push_back(new Model_RqcRLqc_14);
	models.push_back(new Model_RRcdPL2u_15);
	models.push_back(new Model_RRcdPqcL2u_14);
	models.push_back(new Model_RRcL2qc_15);
	models.push_back(new Model_RRcLqc_15);
	models.push_back(new Model_RRcPL_19);
	models.push_back(new Model_RRL_18);
	models.push_back(new Model_RRL_19);
	models.push_back(new Model_RRL2_18);
	models.push_back(new Model_RRL2qc_17);
	models.push_back(new Model_RRLqc_17);
	models.push_back(new Model_RRPEu_19);
	models.push_back(new Model_RRPL_21);
	models.push_back(new Model_RRPL2_20);
	models.push_back(new Model_RRPL2qc_18);
	models.push_back(new Model_RRPL2u_19);
	models.push_back(new Model_RRPL2u_21);

	// prepare output file
	TFile *f_out = TFile::Open("distributions.root", "recreate");

	// evaluate models
	for (const auto &model : models)
	{
		printf("* %s\n", model->name.c_str());

		// prepare graphs
		gDirectory = f_out->mkdir(model->name.c_str());

		TGraph *g_si_p_p = new TGraph(); g_si_p_p->SetName("g_si_p_p");
		TGraph *g_si_p_ap = new TGraph(); g_si_p_ap->SetName("g_si_p_ap");
		TGraph *g_rho_p_p = new TGraph(); g_rho_p_p->SetName("g_rho_p_p");
		TGraph *g_rho_p_ap = new TGraph(); g_rho_p_ap->SetName("g_rho_p_ap");

		// loop over sqrt_s (GeV)
		for (double sqrt_s = 10; sqrt_s < 2E5; sqrt_s *= 1.05)
		{
			const double s = sqrt_s * sqrt_s;

			const int idx = g_si_p_p->GetN();
			g_si_p_p->SetPoint(idx, sqrt_s, model->si_p_p(s));
			g_si_p_ap->SetPoint(idx, sqrt_s, model->si_p_ap(s));
			g_rho_p_p->SetPoint(idx, sqrt_s, model->rho_p_p(s));
			g_rho_p_ap->SetPoint(idx, sqrt_s, model->rho_p_ap(s));
		}

		// save graphs
		g_si_p_p->Write();
		g_si_p_ap->Write();
		g_rho_p_p->Write();
		g_rho_p_ap->Write();
	}

	// clean up
	delete f_out;
	return 0;
}
