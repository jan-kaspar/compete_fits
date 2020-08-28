import root;
import pad_layout;

string f = "../uncertainties.root";

string models[];
/*
models.push("Model_RRdPL2_20");
models.push("Model_RRdPL2u_17");
models.push("Model_RRdPL2u_19");
models.push("Model_RRdPqcL2u_16");
models.push("Model_RqcRcLqc_12");
models.push("Model_RqcRLqc_14");
models.push("Model_RRcdPL2u_15");
models.push("Model_RRcdPqcL2u_14");
models.push("Model_RRcLqc_15");
models.push("Model_RRcPL_19");
models.push("Model_RRL_18");
models.push("Model_RRLnf_19");
models.push("Model_RRLqc_17");
models.push("Model_RRPEu_19");
models.push("Model_RRPL_21");
models.push("Model_RRPL2_20");
models.push("Model_RRPL2qc_18");
models.push("Model_RRPL2u_19");
*/

models.push("Model_RRPnfL2u_21");

models.push("Model_RqcRcL2qc_12");
//models.push("Model_RRL2_18");
//models.push("Model_RRL2qc_17");
//models.push("Model_RRcL2qc_15");

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TGraph_x_max = 50e3;

for (int mi : models.keys)
{
	NewRow();

	NewPad(false);
	label(replace(models[mi], "_", "\_"));

	NewPad("$\sqrt s\ung{GeV}$", "uncertainty of $\si_{\rm pp}\ung{mb}$");
	scale(Log, Linear);
	draw(RootGetObject(f, models[mi] + "/si_p_p/g_unc_stddev"), blue);
	limits((1e1, 0.), (5e4, 2.), Crop);
	yaxis(XEquals(13e3, false), dashed);

	NewPad("$\sqrt s\ung{GeV}$", "uncertainty of $\rh_{\rm pp}$");
	scale(Log, Linear);
	draw(RootGetObject(f, models[mi] + "/rho_p_p/g_unc_stddev"), blue);
	limits((1e1, 0.), (5e4, 0.004), Crop);
	yaxis(XEquals(13e3, false), dashed);

	NewPad("$\De\si_{\rm pp}\ung{mb}$", "$\De \rh_{\rm pp}$");
	draw(RootGetObject(f, models[mi] + "/h_de_rho_p_p_vs_de_si_p_p"), "def");
	limits((-10, -0.005), (+10, +0.005), Crop);

	NewPad("$\sqrt s\ung{GeV}$", "uncertainty band for $\si_{\rm pp}\ung{mb}$");
	scale(Log, Linear);
	draw(RootGetObject(f, models[mi] + "/si_p_p/g_band_up"), red);
	draw(RootGetObject(f, models[mi] + "/si_p_p/g_band_dw"), red);
	limits((1e1, 40.), (5e4, 130.), Crop);
	yaxis(XEquals(13e3, false), dashed);

	NewPad("$\sqrt s\ung{GeV}$", "uncertainty band for $\rh_{\rm pp}$");
	scale(Log, Linear);
	draw(RootGetObject(f, models[mi] + "/rho_p_p/g_band_up"), red);
	draw(RootGetObject(f, models[mi] + "/rho_p_p/g_band_dw"), red);
	limits((1e1, -0.10), (5e4, +0.16), Crop);
	yaxis(XEquals(13e3, false), dashed);
}
