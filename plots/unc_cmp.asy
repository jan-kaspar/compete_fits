import root;
import pad_layout;

string topDir = "../";

string f = topDir + "uncertainty_RRPL2u_21.root";
string f_B = topDir + "uncertainty_RRPL2u_21_B_only.root";

//----------------------------------------------------------------------------------------------------

NewPad("$\sqrt s\ung{GeV}$", "$\si[\si_{\rm pp}]\ung{mb}$");
scale(Log, Linear);

draw(RootGetObject(f_B, "g_stddev_check_si_p_p"), "l", heavygreen, "only B unc., numerical propatation");
draw(RootGetObject(f_B, "g_stddev_si_p_p"), "l", red+dashed, "only B unc., analytic propagation");

draw(RootGetObject(f, "g_stddev_si_p_p"), "l", blue, mCi+1pt+blue, "all unc., numerical propagation");

frame f_legend = BuildLegend();


NewPad("$\sqrt s\ung{GeV}$", "$\si[\si_{\rm p\bar p}]\ung{mb}$");
scale(Log, Linear);

draw(RootGetObject(f, "g_stddev_si_p_ap"), "l", blue, "full");

//----------------------------------------------------------------------------------------------------

NewRow();

NewPad("$\sqrt s\ung{GeV}$", "$\si[\rh_{\rm pp}]\ung{mb}$");
scale(Log, Linear);

draw(RootGetObject(f_B, "g_stddev_check_rho_p_p"), "l", heavygreen, "only B unc., numerical propatation");
draw(RootGetObject(f_B, "g_stddev_rho_p_p"), "l", red+dashed, "only B unc., analytic propagation");

draw(RootGetObject(f, "g_stddev_rho_p_p"), "l", blue, mCi+1pt+blue, "all unc., numerical propagation");


NewPad("$\sqrt s\ung{GeV}$", "$\si[\rh_{\rm p\bar p}]\ung{mb}$");
scale(Log, Linear);

draw(RootGetObject(f, "g_stddev_rho_p_ap"), "l", blue, "full");

NewPad(false);
attach(f_legend);
