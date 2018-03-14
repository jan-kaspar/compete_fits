import root;
import pad_layout;

string f = "../distributions.root";

string models[];
models.push("Model_RRdPL2_20");
models.push("Model_RRdPL2u_17");
models.push("Model_RRdPL2u_19");
models.push("Model_RRdPqcL2u_16");
models.push("Model_RqcRcL2qc_12");
models.push("Model_RqcRcLqc_12");
models.push("Model_RqcRLqc_14");
models.push("Model_RRcdPL2u_15");
models.push("Model_RRcdPqcL2u_14");
models.push("Model_RRcL2qc_15");
models.push("Model_RRcLqc_15");
models.push("Model_RRcPL_19");
models.push("Model_RRL_18");
models.push("Model_RRLnf_19");
models.push("Model_RRL2_18");
models.push("Model_RRL2qc_17");
models.push("Model_RRLqc_17");
models.push("Model_RRPEu_19");
models.push("Model_RRPL_21");
models.push("Model_RRPL2_20");
models.push("Model_RRPL2qc_18");
models.push("Model_RRPL2u_19");
models.push("Model_RRPnfL2u_21");

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TGraph_x_max = 30e3;

for (int mi : models.keys)
{
	NewPage();

	NewPad(false);
	RootObject obj_label = RootGetObject(f, models[mi] + "/g_label");
	string l = obj_label.sExec("GetTitle");
	label("{\SetFontSizesXX$\rm " + l +  "$}");

	NewPad("$\sqrt s\ung{GeV}$", "$\si_{\rm tot}\ung{mb}$");
	scale(Log, Linear);
	draw(RootGetObject(f, models[mi] + "/g_si_p_p"), red, "$\rm pp$");
	draw(RootGetObject(f, models[mi] + "/g_si_p_ap"), red+dashed, "$\rm \bar pp$");
	ylimits(20, 130, Crop);
	AttachLegend(NW, NW);

	NewPad("$\sqrt s\ung{GeV}$", "$\rh_{\rm pp}$");
	currentpad.yTicks = RightTicks(0.05, 0.01);
	scale(Log, Linear);
	draw(RootGetObject(f, models[mi] + "/g_rho_p_p"), blue, "$\rm pp$");
	draw(RootGetObject(f, models[mi] + "/g_rho_p_ap"), blue+dashed, "$\rm \bar pp$");
	ylimits(-0., 0.2, Crop);
	AttachLegend(NW, NW);
}
