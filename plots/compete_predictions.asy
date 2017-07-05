import root;
import pad_layout;

string f = "../distributions.root";

string models[] = {
	"Model_RPdPL2_20",
	"Model_RPdPL2u_17",
	"Model_RPdPL2u_19",
	"Model_RPdPqcL2u_16",

	"Model_RqcRcLqc_12",
	"Model_RqcRcL2qc_12",
};

//----------------------------------------------------------------------------------------------------

void DrawAll(string obj)
{
	for (int mi : models.keys)
		draw(RootGetObject(f, models[mi] + "/" + obj), StdPen(mi), replace(models[mi], "_", "\_"));
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TGraph_x_max = 30e3;

NewPad("$\sqrt s\ung{GeV}$", "$\si_{\rm pp}\ung{mb}$");
scale(Log, Linear);
DrawAll("g_si_p_p");

NewPad("$\sqrt s\ung{GeV}$", "$\si_{\rm \bar pp}\ung{mb}$");
scale(Log, Linear);
DrawAll("g_si_p_ap");

NewRow();

NewPad("$\sqrt s\ung{GeV}$", "$\rh_{\rm pp}$");
scale(Log, Linear);
DrawAll("g_rho_p_p");

NewPad("$\sqrt s\ung{GeV}$", "$\rh_{\rm \bar pp}$");
scale(Log, Linear);
DrawAll("g_rho_p_ap");

//----------------------------------------------------------------------------------------------------

frame f_legend = BuildLegend();

NewPad(false);
attach(f_legend);
