import root;
import pad_layout;

string topDir = "../";

string f = topDir + "uncertainty_RRPL2u_21.root";
string f_B = topDir + "uncertainty_RRPL2u_21_B_only.root";

//----------------------------------------------------------------------------------------------------

NewPad("$\sqrt s\ung{GeV}$", "$\si[\si_{\rm pp}]\ung{mb}$");
scale(Log, Linear);

draw(RootGetObject(f_B, "si_p_p/g_unc_stddev"), "l", heavygreen, "only B unc., numerical propatation");
draw(RootGetObject(f_B, "check/g_stddev_check_si_p_p"), "l", red+dashed, "only B unc., analytic propagation");

draw(RootGetObject(f, "si_p_p/g_unc_stddev"), "l", blue, mCi+1pt+blue, "all unc., numerical propagation");

frame f_legend = BuildLegend();


NewPad("$\sqrt s\ung{GeV}$", "$\si[\si_{\rm p\bar p}]\ung{mb}$");
scale(Log, Linear);

draw(RootGetObject(f, "si_p_ap/g_unc_stddev"), "l", blue, "full");

//----------------------------------------------------------------------------------------------------

NewRow();

NewPad("$\sqrt s\ung{GeV}$", "$\si[\rh_{\rm pp}]\ung{mb}$");
scale(Log, Linear);

draw(RootGetObject(f_B, "rho_p_p/g_unc_stddev"), "l", heavygreen, "only B unc., numerical propatation");
draw(RootGetObject(f_B, "check/g_stddev_check_rho_p_p"), "l", red+dashed, "only B unc., analytic propagation");

draw(RootGetObject(f, "rho_p_p/g_unc_stddev"), "l", blue, mCi+1pt+blue, "all unc., numerical propagation");


NewPad("$\sqrt s\ung{GeV}$", "$\si[\rh_{\rm p\bar p}]\ung{mb}$");
scale(Log, Linear);

draw(RootGetObject(f, "rho_p_ap/g_unc_stddev"), "l", blue, "full");

NewPad(false);
attach(f_legend);
