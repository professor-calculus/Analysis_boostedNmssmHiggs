#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TAxis.h>

TStyle * TDRStyle()
{
	TStyle * tdrStyle = new TStyle("tdrStyle","");
	gROOT->SetStyle("tdrStyle"); 
	//tdrStyle->SetPalette(palette);

	// For the canvas:
	tdrStyle->SetCanvasBorderMode(0);
	tdrStyle->SetCanvasColor(kWhite);
	tdrStyle->SetCanvasDefH(600); //Height of canvas
	tdrStyle->SetCanvasDefW(800); //Width of canvas
	tdrStyle->SetCanvasDefX(0);   //POsition on screen
	tdrStyle->SetCanvasDefY(0);

	// For the Pad:
	tdrStyle->SetPadBorderMode(0);
	// tdrStyle->SetPadBorderSize(Width_t size = 1);
	tdrStyle->SetPadColor(kWhite);
	tdrStyle->SetPadGridX(false);
	tdrStyle->SetPadGridY(false);
	tdrStyle->SetGridColor(0);
	tdrStyle->SetGridStyle(3);
	tdrStyle->SetGridWidth(1);
	tdrStyle->SetPadGridX(false);
	tdrStyle->SetPadGridY(false);

	// For the frame:
	tdrStyle->SetFrameBorderMode(0);
	tdrStyle->SetFrameBorderSize(1);
	tdrStyle->SetFrameFillColor(0);
	tdrStyle->SetFrameFillStyle(0);
	tdrStyle->SetFrameLineColor(1);
	tdrStyle->SetFrameLineStyle(1);
	tdrStyle->SetFrameLineWidth(1);

	// For the histo:
	// tdrStyle->SetHistFillColor(1);
	// tdrStyle->SetHistFillStyle(0);
	tdrStyle->SetHistLineColor(1);
	tdrStyle->SetHistLineStyle(0);
	tdrStyle->SetHistLineWidth(1);
	// tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
	// tdrStyle->SetNumberContours(Int_t number = 20);

	tdrStyle->SetEndErrorSize(2);
	//  tdrStyle->SetErrorMarker(20);
	tdrStyle->SetErrorX(0.);

	tdrStyle->SetMarkerStyle(1); // EDITFROM 20

	//For the legend
	tdrStyle->SetLegendBorderSize(0);
	tdrStyle->SetLegendFillColor(0);
	tdrStyle->SetLegendFont(42); // EDITFROM 42
	tdrStyle->SetLegendTextSize(0.04);

	//For the fit/function:
	tdrStyle->SetOptFit(0);
	tdrStyle->SetFitFormat("5.4g");
	tdrStyle->SetFuncColor(2);
	tdrStyle->SetFuncStyle(1);
	tdrStyle->SetFuncWidth(1);

	//For the date:
	tdrStyle->SetOptDate(0);
	// tdrStyle->SetDateX(Float_t x = 0.01);
	// tdrStyle->SetDateY(Float_t y = 0.01);

	// For the statistics box:
	tdrStyle->SetOptFile(0);
	tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
	tdrStyle->SetStatColor(kWhite);
	tdrStyle->SetStatFont(42);
	tdrStyle->SetStatFontSize(0.025);
	tdrStyle->SetStatTextColor(1);
	tdrStyle->SetStatFormat("6.4g");
	tdrStyle->SetStatBorderSize(1);
	tdrStyle->SetStatH(0.1);
	tdrStyle->SetStatW(0.15);
	// tdrStyle->SetStatStyle(Style_t style = 1001);
	// tdrStyle->SetStatX(Float_t x = 0);
	// tdrStyle->SetStatY(Float_t y = 0);

	// Margins:
	//tdrStyle->SetPadTopMargin(0.05);
	tdrStyle->SetPadTopMargin(0.10);
	tdrStyle->SetPadBottomMargin(0.13);
	// tdrStyle->SetPadLeftMargin(0.16);
	tdrStyle->SetPadLeftMargin(0.14);
	// tdrStyle->SetPadRightMargin(0.02);
	tdrStyle->SetPadRightMargin(0.07);
	// tdrStyle->SetPadRightMargin(0.10); // only really want to be changing this for the scatters

	// For the Global title:
	// tdrStyle->SetOptTitle(1); //EDITFROM 0
	tdrStyle->SetTitleFont(42);
	tdrStyle->SetTitleColor(1);
	tdrStyle->SetTitleTextColor(1);
	tdrStyle->SetTitleFillColor(10);
	tdrStyle->SetTitleFontSize(0.05);
	// tdrStyle->SetTitleH(0); // Set the height of the title box
	// tdrStyle->SetTitleW(0); // Set the width of the title box
	// tdrStyle->SetTitleX(0); // Set the position of the title box
	// tdrStyle->SetTitleY(0.985); // Set the position of the title box
	// tdrStyle->SetTitleStyle(Style_t style = 1001);
	// tdrStyle->SetTitleBorderSize(2);

	// For the axis titles:
	tdrStyle->SetTitleColor(1, "XYZ");
	tdrStyle->SetTitleFont(42, "XYZ");
	tdrStyle->SetTitleSize(0.045, "XYZ");
	// tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
	// tdrStyle->SetTitleYSize(Float_t size = 0.02);
	tdrStyle->SetTitleXOffset(0.95);//EDITFROM 0.9
	tdrStyle->SetTitleYOffset(1.35);//EDITFROM 1.25
	// tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

	// For the axis labels:
	tdrStyle->SetLabelColor(1, "XYZ");
	tdrStyle->SetLabelFont(42, "XYZ");
	tdrStyle->SetLabelOffset(0.007, "XYZ");
	// tdrStyle->SetLabelSize(0.05, "XYZ");
	tdrStyle->SetLabelSize(0.045, "XYZ");

	// For the axis:
	tdrStyle->SetAxisColor(1, "XYZ");
	tdrStyle->SetStripDecimals(kTRUE);
	tdrStyle->SetTickLength(0.03, "XYZ");
	tdrStyle->SetNdivisions(510, "XYZ");
	tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	tdrStyle->SetPadTickY(1);

	// Change for log plots:
	tdrStyle->SetOptLogx(0);
	tdrStyle->SetOptLogy(0);
	tdrStyle->SetOptLogz(0);

	// Postscript options:
	tdrStyle->SetPaperSize(20.,20.);
	// tdrStyle->SetLineScalePS(Float_t scale = 3);
	// tdrStyle->SetLineStyleString(Int_t i, const char* text);
	// tdrStyle->SetHeaderPS(const char* header);
	// tdrStyle->SetTitlePS(const char* pstitle);

	// tdrStyle->SetBarOffset(Float_t baroff = 0.5);
	// tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
	// tdrStyle->SetPaintTextFormat(const char* format = "g");
	// tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
	// tdrStyle->SetTimeOffset(Double_t toffset);
	// tdrStyle->SetHistMinimumZero(kTRUE);

	//tdrStyle->cd();
	return tdrStyle;
}

void drawAndSave(TMultiGraph * mg, TLegend * leg, std::string saveName, std::string yAxisTitle, double mass_squark, double mass_higgs, double legXmin, double legXmax, double legYmin, double legYmax)
{
	TStyle * tdrStyle = TDRStyle();
	TCanvas* cPlot=new TCanvas("cPlot","cPlot");
	mg->Draw("a3L");
	mg->GetXaxis()->SetTitle("mass_Higgs / mass_NLSP");
	mg->GetYaxis()->SetTitle(yAxisTitle.c_str());
	// mg->SetTitle(Form("Mass_{squark} = %.1f, Mass_{higgs} = %.1f", mass_squark, mass_higgs));
	mg->Draw("a3L");
	leg->SetX1NDC(legXmin);
	leg->SetX2NDC(legXmax);
	leg->SetY1NDC(legYmin);
	leg->SetY2NDC(legYmax);
	leg->Draw("same");
	TLatex * latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);
	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,Form("Mass_Squark = %.0fGeV; Mass_Higgs =  %.0fGeV", mass_squark, mass_higgs));
	cPlot->SaveAs(saveName.c_str());
	cPlot->Close();
}

int SetColor(int posititon, int maxColors)
{	
	gStyle->SetPalette(55); // sets what sort of colours we will use (this number could be changed)

	// modifier is an offset in the colour spectrum
	double modifier = 0.00;
    double colorIndex;
    int colour(1);
    double fraction = (double)(posititon)/(double)(maxColors);

    if( posititon > maxColors || posititon < 0 || maxColors < 0 ) colour = 1;
    else
    {
        colorIndex = (fraction + modifier) * gStyle->GetNumberOfColors();
        colour = gStyle->GetColorPalette(colorIndex);
    }
    return colour;
}

void SetGraphOptions(TGraphAsymmErrors * gr, size_t index, size_t numberOfBins){
	gr->SetLineWidth(2);
	gr->SetLineColor(SetColor(index, numberOfBins));
	gr->SetFillColorAlpha(SetColor(index, numberOfBins), 0.35);
	gr->SetFillStyle(3001);
	gr->SetTitle("testingTitleGR");
}






void massParamsToKinematics(){

	/*
	We calculate kintematic properties of the following interaction:
	pp -> squark + squark

	squark -> quark + NLSP;
	NLSP -> LSP + higgs;
	higgs -> bb;

	interesting properties:

	quark: energy/momentum, transverse momentum
	NLSP: energy, momentum, transverse momentum

	LSP:
	in rest frame of NLSP - energy, momentum
	in lab frame - energy, momentum, transverse momentum

	Higgs:
	in rest frame of NLSP - energy, momentum
	in lab frame - energy, momentum, transverse momentum

	Splitting angle of the bbar pair
	(if LSP pt small, it will carry the NLSP momentum)
	*/


	// --- USER INPUT --- //
	// setup the four mass parameters
	double mass_squark = 1000.0;
	double mass_higgs = 83.0;
	double mass_ratio_beginPoint = 0.7; // mass_ratio = mass_higgs / mass_NLSP
	double mass_ratio_stepSize = 0.01;
	std::vector<double> mass_delta_vec = {15.0, 10.0, 5.0, 1.0}; // mass_delta = mass_NLSP - mass_higgs - mass_LSP

	TMultiGraph * mg_mass_NLSP = new TMultiGraph();
	TMultiGraph * mg_mass_LSP = new TMultiGraph();
	TMultiGraph * mg_energy_quark = new TMultiGraph();
	TMultiGraph * mg_momentum_quark = new TMultiGraph();
	TMultiGraph * mg_energy_NLSP = new TMultiGraph();
	TMultiGraph * mg_momentum_NLSP = new TMultiGraph();
	TMultiGraph * mg_energy_LSP = new TMultiGraph();
	TMultiGraph * mg_momentum_LSP = new TMultiGraph();
	TMultiGraph * mg_energy_higgs = new TMultiGraph();
	TMultiGraph * mg_momentum_higgs = new TMultiGraph();
	TMultiGraph * mg_seperation_bb = new TMultiGraph();

	TLegend * legend = new TLegend();
	legend->SetHeader("Mass Gap (GeV)");
	TCanvas* c=new TCanvas("c","c");



	for (size_t i_delta = 0; i_delta < mass_delta_vec.size(); ++i_delta){
		double mass_delta = mass_delta_vec[i_delta];

		// WHAT TO DO ABOUT TRANVERSE MOMENTUM???
		std::vector<double> mass_ratio_vec;
		std::vector<double> null_vec;
		std::vector<double> mass_NLSP_vec;
		std::vector<double> mass_LSP_vec;
		std::vector<double> energy_quark_vec;
		std::vector<double> momentum_quark_vec;
		std::vector<double> energy_NLSP_vec;
		std::vector<double> momentum_NLSP_vec;
		std::vector<double> energy_LSP_mean_vec;
		std::vector<double> energy_LSP_dev_vec;
		std::vector<double> momentum_LSP_mean_vec;
		std::vector<double> momentum_LSP_devUp_vec;
		std::vector<double> momentum_LSP_devDown_vec;
		std::vector<double> energy_higgs_mean_vec;
		std::vector<double> energy_higgs_dev_vec;
		std::vector<double> momentum_higgs_mean_vec;
		std::vector<double> momentum_higgs_devUp_vec;
		std::vector<double> momentum_higgs_devDown_vec;
		std::vector<double> seperation_bb_vec;

		for (double mass_ratio = mass_ratio_beginPoint; mass_ratio < mass_higgs / (mass_delta + mass_higgs); mass_ratio += mass_ratio_stepSize){  
			mass_ratio_vec.push_back(mass_ratio);
			null_vec.push_back(0);  
			
			//////////////////
			// ------------ // 
			// CALCULATIONS //
			// ------------ //
			//////////////////

			// calculate the LSP and NLSP masses
			double mass_NLSP = mass_higgs / mass_ratio;
			double mass_LSP = mass_NLSP - mass_higgs - mass_delta;

			// squark -> quark + NLSP
			double energy_quark = (mass_squark * mass_squark - mass_NLSP * mass_NLSP) / (2 * mass_squark);
			double momentum_quark = energy_quark;
			double energy_NLSP = mass_squark - energy_quark;
			double momentum_NLSP = momentum_quark;

			// NLSP -> LSP + higgs
			double boost_gamma = energy_NLSP / mass_NLSP;
			double boost_betaGamma = momentum_NLSP / mass_NLSP;

			// LSP: rest frame
			double energy_LSP_restFrame = (mass_NLSP * mass_NLSP + mass_LSP * mass_LSP - mass_higgs * mass_higgs) / (2 * mass_NLSP);
			double momentum_LSP_restFrame = sqrt(energy_LSP_restFrame * energy_LSP_restFrame - mass_LSP * mass_LSP);
			// higgs: rest frame
			double energy_higgs_restFrame = mass_NLSP - energy_LSP_restFrame;
			double momentum_higgs_restFrame = momentum_LSP_restFrame;
			// LSP: lab frame                                                     //ARE THE CALCULATIONS RIGHT??
			double energy_LSP_mean = boost_gamma * energy_LSP_restFrame;
			double energy_LSP_dev = boost_betaGamma * momentum_LSP_restFrame; // energy extremums can be +/- this value depending on boost angle
			double momentum_LSP_mean = sqrt(energy_LSP_mean * energy_LSP_mean - mass_LSP * mass_LSP);
			double momentum_LSP_devUp = sqrt( (energy_LSP_mean+energy_LSP_dev) * (energy_LSP_mean+energy_LSP_dev) - mass_LSP * mass_LSP ) - momentum_LSP_mean;
			double momentum_LSP_devDown = -1 * sqrt( (energy_LSP_mean-energy_LSP_dev) * (energy_LSP_mean-energy_LSP_dev) - mass_LSP * mass_LSP ) + momentum_LSP_mean;
			// higgs: lab frame
			double energy_higgs_mean = boost_gamma * energy_higgs_restFrame;
			double energy_higgs_dev = boost_betaGamma * momentum_higgs_restFrame; // energy extremums can be +/- this value depending on boost angle
			double momentum_higgs_mean = sqrt(energy_higgs_mean * energy_higgs_mean - mass_higgs * mass_higgs);
			double momentum_higgs_devUp = sqrt( (energy_higgs_mean+energy_higgs_dev) * (energy_higgs_mean+energy_higgs_dev) - mass_higgs * mass_higgs ) - momentum_higgs_mean;
			double momentum_higgs_devDown = -1 * sqrt( (energy_higgs_mean-energy_higgs_dev) * (energy_higgs_mean-energy_higgs_dev) - mass_higgs * mass_higgs ) + momentum_higgs_mean;

			// higgs -> bbar
			double seperation_bb = 2 * mass_higgs / momentum_higgs_mean; //transverse momentum???

			//////////////////
			// -----END---- // 
			// CALCULATIONS //
			// ------------ //
			//////////////////

			mass_NLSP_vec.push_back(mass_NLSP);
			mass_LSP_vec.push_back(mass_LSP);
			energy_quark_vec.push_back(energy_quark);
			momentum_quark_vec.push_back(momentum_quark);
			energy_NLSP_vec.push_back(energy_NLSP);
			momentum_NLSP_vec.push_back(momentum_NLSP);
			energy_LSP_mean_vec.push_back(energy_LSP_mean);
			energy_LSP_dev_vec.push_back(energy_LSP_dev);
			momentum_LSP_mean_vec.push_back(momentum_LSP_mean);
			momentum_LSP_devUp_vec.push_back(momentum_LSP_devUp);
			momentum_LSP_devDown_vec.push_back(momentum_LSP_devDown);
			energy_higgs_mean_vec.push_back(energy_higgs_mean);
			energy_higgs_dev_vec.push_back(energy_higgs_dev);
			momentum_higgs_mean_vec.push_back(momentum_higgs_mean);
			momentum_higgs_devUp_vec.push_back(momentum_higgs_devUp);
			momentum_higgs_devDown_vec.push_back(momentum_higgs_devDown);
			seperation_bb_vec.push_back(seperation_bb);

		} // closes loop through the mass_ratio values

		TGraphAsymmErrors * gr_mass_NLSP = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(mass_NLSP_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]));
		TGraphAsymmErrors * gr_mass_LSP = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(mass_LSP_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]));
		TGraphAsymmErrors * gr_energy_quark = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(energy_quark_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]));
		TGraphAsymmErrors * gr_momentum_quark = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(momentum_quark_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]));
		TGraphAsymmErrors * gr_energy_NLSP = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(energy_NLSP_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]));
		TGraphAsymmErrors * gr_momentum_NLSP = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(momentum_NLSP_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]));
		TGraphAsymmErrors * gr_energy_LSP = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(energy_LSP_mean_vec[0]), &(null_vec[0]), &(null_vec[0]), &(energy_LSP_dev_vec[0]), &(energy_LSP_dev_vec[0]));
		TGraphAsymmErrors * gr_momentum_LSP = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(momentum_LSP_mean_vec[0]), &(null_vec[0]), &(null_vec[0]), &(momentum_LSP_devDown_vec[0]), &(momentum_LSP_devUp_vec[0]));
		TGraphAsymmErrors * gr_energy_higgs = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(energy_higgs_mean_vec[0]), &(null_vec[0]), &(null_vec[0]), &(energy_higgs_dev_vec[0]), &(energy_higgs_dev_vec[0]));
		TGraphAsymmErrors * gr_momentum_higgs = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(momentum_higgs_mean_vec[0]), &(null_vec[0]), &(null_vec[0]), &(momentum_higgs_devDown_vec[0]), &(momentum_higgs_devUp_vec[0]));
		TGraphAsymmErrors * gr_seperation_bb = new TGraphAsymmErrors(mass_ratio_vec.size(), &(mass_ratio_vec[0]), &(seperation_bb_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]), &(null_vec[0]));

		SetGraphOptions(gr_mass_NLSP, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_mass_LSP, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_energy_quark, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_momentum_quark, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_energy_NLSP, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_momentum_NLSP, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_energy_LSP, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_momentum_LSP, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_energy_higgs, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_momentum_higgs, i_delta, mass_delta_vec.size());
		SetGraphOptions(gr_seperation_bb, i_delta, mass_delta_vec.size());

		TH1F * h_dummy = new TH1F("h_dummy", "", 100, 0, 100);
		h_dummy->SetLineColor(SetColor(i_delta, mass_delta_vec.size()));
		legend->AddEntry(h_dummy, Form("%.1f", mass_delta_vec[i_delta]), "L");
	    legend->SetLineColor(0);

		// add the TGraphs for this mass_delta bin to the TMultiGraphs
		mg_mass_NLSP->Add(gr_mass_NLSP);
		mg_mass_LSP->Add(gr_mass_LSP);
		mg_energy_quark->Add(gr_energy_quark);
		mg_momentum_quark->Add(gr_momentum_quark);
		mg_energy_NLSP->Add(gr_energy_NLSP);
		mg_momentum_NLSP->Add(gr_momentum_NLSP);
		mg_energy_LSP->Add(gr_energy_LSP);
		mg_momentum_LSP->Add(gr_momentum_LSP);
		mg_energy_higgs->Add(gr_energy_higgs);
		mg_momentum_higgs->Add(gr_momentum_higgs);
		mg_seperation_bb->Add(gr_seperation_bb);

	} // closes loop through the mass_delta entries

	drawAndSave(mg_mass_NLSP, legend, "NLSP_mass.pdf", "mass_NLSP (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.65, 0.86);
	drawAndSave(mg_mass_LSP, legend, "LSP_mass.pdf", "mass_LSP (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.65, 0.86);
	// drawAndSave(mg_energy_quark, legend, "quark_energy.pdf", "energy_quark (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.65, 0.86);
	drawAndSave(mg_momentum_quark, legend, "quark_momentum.pdf", "momentum_quark (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.25, 0.46);
	drawAndSave(mg_energy_NLSP, legend, "NLSP_energy.pdf", "energy_NLSP (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.65, 0.86);
	drawAndSave(mg_momentum_NLSP, legend, "NLSP_momentum.pdf", "momentum_NLSP (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.25, 0.46);
	// drawAndSave(mg_energy_LSP, legend, "LSP_energy.pdf", "energy_LSP (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.65, 0.86);
	drawAndSave(mg_momentum_LSP, legend, "LSP_momentum.pdf", "momentum_LSP (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.65, 0.86);
	// drawAndSave(mg_energy_higgs, legend, "higgs_energy.pdf", "energy_higgs (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.65, 0.86);
	drawAndSave(mg_momentum_higgs, legend, "higgs_momentum.pdf", "momentum_higgs (GeV)", mass_squark, mass_higgs, 0.70, 0.88, 0.25, 0.46);
	drawAndSave(mg_seperation_bb, legend, "speration_bb.pdf", "seperation_bb", mass_squark, mass_higgs, 0.70, 0.88, 0.65, 0.86);

}































