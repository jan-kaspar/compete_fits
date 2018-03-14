import root;
import pad_layout;

string f = "../distributions.root";

string models[];
pen m_pens[];

models.push("Model_RRdPL2_20"); m_pens.push(blue);
models.push("Model_RRdPL2u_17"); m_pens.push(blue);
models.push("Model_RRdPL2u_19"); m_pens.push(blue);
models.push("Model_RRdPqcL2u_16"); m_pens.push(blue);
models.push("Model_RqcRcL2qc_12"); m_pens.push(magenta);
models.push("Model_RqcRcLqc_12"); m_pens.push(heavygreen);
models.push("Model_RqcRLqc_14"); m_pens.push(heavygreen);
models.push("Model_RRcdPL2u_15"); m_pens.push(blue);
models.push("Model_RRcdPqcL2u_14"); m_pens.push(blue);
models.push("Model_RRcL2qc_15"); m_pens.push(magenta);
models.push("Model_RRcLqc_15"); m_pens.push(heavygreen);
models.push("Model_RRcPL_19"); m_pens.push(heavygreen);
models.push("Model_RRL_18"); m_pens.push(heavygreen);
models.push("Model_RRLnf_19"); m_pens.push(heavygreen);
models.push("Model_RRL2_18"); m_pens.push(magenta);
models.push("Model_RRL2qc_17"); m_pens.push(magenta);
models.push("Model_RRLqc_17"); m_pens.push(heavygreen);
models.push("Model_RRPEu_19"); m_pens.push(blue + dashed);
models.push("Model_RRPL_21"); m_pens.push(heavygreen);
models.push("Model_RRPL2_20"); m_pens.push(heavygreen);
models.push("Model_RRPL2qc_18"); m_pens.push(heavygreen);
models.push("Model_RRPL2u_19"); m_pens.push(blue);
models.push("Model_RRPnfL2u_21"); m_pens.push(blue);

xSizeDef = 10cm;
ySizeDef = xSizeDef * 2/3;

legendLabelPen = fontcommand("\SetFontSizesVI");

drawGridDef = false;


//----------------------------------------------------------------------------------------------------

void DrawAll(string obj)
{
	for (int mi : models.keys)
	{
		RootObject obj_label = RootGetObject(f, models[mi] + "/g_label");
		string l = obj_label.sExec("GetTitle");
		draw(RootGetObject(f, models[mi] + "/" + obj), m_pens[mi], "$\rm " + l + "$");
	}
}

//----------------------------------------------------------------------------------------------------

// fine shift
real fsh=0;

void DrawPoint(real W, real si, real em, real ep=em, pen col=red, marker m)
{
	col += squarecap;

	draw(shift(fsh, 0)*(Scale((W, si-em))--Scale((W, si+ep))), col);
	draw(shift(fsh, 0)*Scale((W, si)), m);

	// reset fine shift
	fsh = 0;
}

//----------------------------------------------------------------------------------------------------

void DrawAxes(real y_label)
{
	yaxis(XEquals(0.546e3, false), dotted + roundcap);
	yaxis(XEquals(0.9e3, false), dotted + roundcap);
	yaxis(XEquals(1.8e3, false), dotted + roundcap);
	yaxis(XEquals(2.76e3, false), dotted + roundcap);
	yaxis(XEquals(7e3, false), dotted + roundcap);
	yaxis(XEquals(8e3, false), dotted + roundcap);
	yaxis(XEquals(13e3, false), dotted + roundcap);

	label(rotate(90)*Label("\SmallerFonts$0.546\un{TeV}$"), Scale((0.546e3, y_label)), Fill(white));
	label(rotate(90)*Label("\SmallerFonts$0.9\un{TeV}$"), Scale((0.9e3, y_label)), Fill(white));
	label(rotate(90)*Label("\SmallerFonts$1.8\un{TeV}$"), Scale((1.8e3, y_label)), Fill(white));
	label(rotate(90)*Label("\SmallerFonts$2.76\un{TeV}$"), Scale((2.76e3, y_label)), Fill(white));
	label(rotate(90)*Label("\SmallerFonts$7\un{TeV}$"), Scale((7e3-300, y_label)), Fill(white));
	label(rotate(90)*Label("\SmallerFonts$8\un{TeV}$"), Scale((8e3+300, y_label)), Fill(white));
	label(rotate(90)*Label("\SmallerFonts$13\un{TeV}$"), Scale((13e3, y_label)), Fill(white));
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TGraph_x_max = 1.1e5;

NewPad("$\sqrt s\ung{GeV}$", "$\si_{\rm pp}\ung{mb}$");
currentpad.drawGridY = true;
currentpad.yTicks = RightTicks(10., 5.);
scale(Log, Linear);

DrawAll("g_si_p_p");

real fshu = 0.02;

fsh = +0fshu; DrawPoint(2.76e3, 84.7, 3.3, red, mCi+2pt+red);

fsh = -1.0fshu; DrawPoint(7e3, 98.0, 2.5, red, mCi+2pt+red);
fsh = +1.0fshu; DrawPoint(7e3, 98.6, 2.2, red, mCi+2pt+red);

fsh = -1.0fshu; DrawPoint(8e3, 101.5, 2.1, red, mCi+2pt+red);
fsh = +1.0fshu; DrawPoint(8e3, 102.9, 2.3, red, mCi+2pt+red);

fsh = +0fshu; DrawPoint(13e3, 110.627, 3.0, red, mCi+2pt+red);

limits((1e1, 40), (1e5, 150), Crop);

DrawAxes(130);


//----------------------------------------------------------------------------------------------------

NewRow();

NewPad("$\sqrt s\ung{GeV}$", "$\rh_{\rm pp}$");
currentpad.drawGridY = true;
currentpad.yTicks = RightTicks(0.05, 0.01);
scale(Log, Linear);

AddToLegend("TOTEM measurements", red, mCi+2pt+red);

AddToLegend("<{\it COMPETE models:}");

DrawAll("g_rho_p_p");
limits((1e1, 0.), (1e5, 0.2), Crop);

DrawAxes(0.05);

DrawPoint(8e3, 0.12, 0.03, red+0.8pt, mCi+true+2pt+red);

DrawPoint(13e3, 0.10, 0.01, 0.01, red+0.8pt, mCi+true+2pt+red);

//----------------------------------------------------------------------------------------------------

frame f_legend = BuildLegend();

NewPad(false);
//attach(shift(0, 120)*f_legend);
attach(f_legend);
FixPad(270, 0);

GShipout(margin=3mm, vSkip=1mm);
