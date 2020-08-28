import root;
import pad_layout;

string topDir = "../";

string objects[] = {
	"g_cen_val",
	"g_unc_mean",
	"g_unc_stddev",
	"g_band_up",
	"g_band_dw",
};

string dirs[] = {
	"si_p_p/",
	"si_p_ap/",
	"rho_p_p/",
	"rho_p_ap/",
};


//----------------------------------------------------------------------------------------------------

for (int di : dirs.keys)
{
	string dir = dirs[di];

	NewRow();

	for (int oi : objects.keys)
	{
		NewPad();
		scale(Log, Linear);
	
		draw(RootGetObject(topDir + "uncertainty_RRPnfL2u_21.root", dir + objects[oi]), blue);
	
		draw(RootGetObject(topDir + "uncertainties.root", "Model_RRPnfL2u_21/" + dir + objects[oi]), red+dashed);
	}
}
