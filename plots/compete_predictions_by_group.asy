import root;
import pad_layout;

string f = "../distributions.root";

string models[];
pen m_pens[];

models.push("Model_RPdPL2_20"); m_pens.push(blue);
models.push("Model_RPdPL2u_17"); m_pens.push(blue);
models.push("Model_RPdPL2u_19"); m_pens.push(blue);
models.push("Model_RPdPqcL2u_16"); m_pens.push(blue);
models.push("Model_RqcRcL2qc_12"); m_pens.push(magenta);
models.push("Model_RqcRcLqc_12"); m_pens.push(heavygreen);
models.push("Model_RqcRLqc_14"); m_pens.push(heavygreen);
models.push("Model_RRcdPL2u_15"); m_pens.push(blue);
models.push("Model_RRcdPqcL2u_14"); m_pens.push(blue);
models.push("Model_RRcL2qc_15"); m_pens.push(magenta);
models.push("Model_RRcLqc_15"); m_pens.push(heavygreen);
models.push("Model_RRcPL_19"); m_pens.push(heavygreen);
models.push("Model_RRL_18"); m_pens.push(heavygreen);
models.push("Model_RRL_19"); m_pens.push(heavygreen);
models.push("Model_RRL2_18"); m_pens.push(magenta);
models.push("Model_RRL2qc_17"); m_pens.push(magenta);
models.push("Model_RRLqc_17"); m_pens.push(heavygreen);
models.push("Model_RRPEu_19"); m_pens.push(blue + dashed);
models.push("Model_RRPL_21"); m_pens.push(heavygreen);
models.push("Model_RRPL2_20"); m_pens.push(heavygreen);
models.push("Model_RRPL2qc_18"); m_pens.push(heavygreen);
models.push("Model_RRPL2u_19"); m_pens.push(blue);
models.push("Model_RRPL2u_21"); m_pens.push(blue);



//----------------------------------------------------------------------------------------------------

void DrawAll(string obj)
{
	for (int mi : models.keys)
	{
		draw(RootGetObject(f, models[mi] + "/" + obj), m_pens[mi], replace(models[mi], "_", "\_"));
	}
}

//----------------------------------------------------------------------------------------------------

void DrawPoint(real W, real si, real em, real ep, pen col=red, marker m)
{
	draw((Scale((W, si-em))--Scale((W, si+ep))), col);
	draw(Scale((W, si)), m);
}

//----------------------------------------------------------------------------------------------------

void DrawPointRel(real W, real si, real re, pen col=red, marker m)
{
	draw((Scale((W, si*(1-re/100)))--Scale((W, si*(1+re/100)))), col);
	draw(Scale((W, si)), m);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TGraph_x_max = 30e3;

NewPad("$\sqrt s\ung{GeV}$", "$\si_{\rm pp}\ung{mb}$");
scale(Log, Linear);
DrawAll("g_si_p_p");

DrawPoint(2.76e3, 84.7, 3.3, 3.3, red+0.8pt, mCi+true+1.8pt+red);
DrawPointRel(7e3, 98.1, 2.4, red+0.8pt, mCi+true+1.8pt+red);
DrawPointRel(8e3, 102, 2.8, red+0.8pt, mCi+true+1.8pt+red);	

ylimits(40, 130, Crop);

yTicksDef = RightTicks(0.05, 0.01);

NewPad("$\sqrt s\ung{GeV}$", "$\rh_{\rm pp}$");
scale(Log, Linear);
DrawAll("g_rho_p_p");
ylimits(-0., 0.2, Crop);

DrawPoint(13e3, 0.10, 0.01, 0.01, red+0.8pt, mCi+true+1.8pt+red);

//----------------------------------------------------------------------------------------------------

frame f_legend = BuildLegend();

NewPad(false);
attach(f_legend);
