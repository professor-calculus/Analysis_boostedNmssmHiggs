#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TAxis.h>



void drawAndSave(TMultiGraph * mg, TLegend * leg, std::string saveName)
{

	TCanvas* cPlot=new TCanvas("cPlot","cPlot");


// a function to draw and save!  use tdr style also
 	
// add latex stamps

	mg->Draw("ALP");
	mg->GetXaxis()->SetTitle("massHigss / massNLSP");
	mg->GetYaxis()->SetTitle("massLSP (GeV)");
	mg->Draw("ALP");

leg->Draw();
cPlot->SaveAs(saveName.c_str());
}










// maybe get these out of the plotting class...
int SetColor(int posititon, int maxColors)
{
// the plot has 'maxColors' number of colours involved
// this function gives you a colour for the nth histogram
	
	gStyle->SetPalette(55); // sets what sort of colours we will use (this number could be changed)
	
	// nice for four inputs
	// for three inputs use 'position+1' and 'size+1' for best outputs

	// modifier is an offset in the colour spectrum
	double modifier = 0.00;
	// double modifier = 0.05;
	// double modifier = 0.10;
    // double modifier = 0.15;
    // double modifier = 0.20;
	// double modifier = 0.25;	
    // double modifier = 0.30;
    double colorIndex;
    int colour(1);
    // double fraction = (double)(posititon)/(double)(maxColors-1);
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

	gr->SetLineColor(SetColor(index, numberOfBins));
	gr->SetLineWidth(2);
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
	std::vector<double> mass_delta_vec = {1.0, 5.0, 10.0, 15.0}; // mass_delta = mass_NLSP - mass_higgs - mass_LSP

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

	// maybe just use one legend as they are all the same anyway
	TLegend * legend = new TLegend(0.1, 0.7, 0.48, 0.9);
	legend->SetHeader("Mass Gap (GeV)");

	// TLegend * leg_mass_NLSP = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_mass_LSP = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_energy_quark = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_momentum_quark = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_energy_NLSP = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_momentum_NLSP = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_energy_LSP = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_momentum_LSP = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_energy_higgs = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_momentum_higgs = new TLegend(0.1, 0.7, 0.48, 0.9);
	// TLegend * leg_seperation_bb = new TLegend(0.1, 0.7, 0.48, 0.9);

	// leg_mass_NLSP->SetHeader("Mass Gap (GeV)");
	// leg_mass_LSP->SetHeader("Mass Gap (GeV)");
	// leg_energy_quark->SetHeader("Mass Gap (GeV)");
	// leg_momentum_quark->SetHeader("Mass Gap (GeV)");
	// leg_energy_NLSP->SetHeader("Mass Gap (GeV)");
	// leg_momentum_NLSP->SetHeader("Mass Gap (GeV)"); 
	// leg_energy_LSP->SetHeader("Mass Gap (GeV)");
	// leg_momentum_LSP->SetHeader("Mass Gap (GeV)");
	// leg_energy_higgs->SetHeader("Mass Gap (GeV)");
	// leg_momentum_higgs->SetHeader("Mass Gap (GeV)");
	// leg_seperation_bb->SetHeader("Mass Gap (GeV)");

	TCanvas* c=new TCanvas("c","c");
	// TCanvas * c_mass_NLSP = new TCanvas("c", "c");
	// TCanvas * c_mass_LSP = new TCanvas("c", "c");
	// TCanvas * c_energy_quark = new TCanvas("c", "c");
	// TCanvas * c_momentum_quark = new TCanvas("c", "c");
	// TCanvas * c_energy_NLSP = new TCanvas("c", "c");
	// TCanvas * c_momentum_NLSP = new TCanvas("c", "c");
	// TCanvas * c_energy_LSP = new TCanvas("c", "c");
	// TCanvas * c_momentum_LSP = new TCanvas("c", "c");
	// TCanvas * c_energy_higgs = new TCanvas("c", "c");
	// TCanvas * c_momentum_higgs = new TCanvas("c", "c");
	// TCanvas * c_seperation_bb = new TCanvas("c", "c");



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
			// LSP: lab frame
			double energy_LSP_mean = boost_gamma * energy_LSP_restFrame;
			double energy_LSP_dev = boost_betaGamma * momentum_LSP_restFrame; // energy extremums can be +/- this value depending on boost angle
			double momentum_LSP_mean = sqrt(energy_LSP_mean * energy_LSP_mean - mass_LSP * mass_LSP);
			double momentum_LSP_devUp = sqrt( (energy_LSP_mean+energy_LSP_dev) * (energy_LSP_mean+energy_LSP_dev) - mass_LSP * mass_LSP );
			double momentum_LSP_devDown = sqrt( (energy_LSP_mean-energy_LSP_dev) * (energy_LSP_mean-energy_LSP_dev) - mass_LSP * mass_LSP );
			// higgs: lab frame
			double energy_higgs_mean = boost_gamma * energy_higgs_restFrame;
			double energy_higgs_dev = boost_betaGamma * momentum_higgs_restFrame; // energy extremums can be +/- this value depending on boost angle
			double momentum_higgs_mean = sqrt(energy_higgs_mean * energy_higgs_mean - mass_higgs * mass_higgs);
			double momentum_higgs_devUp = sqrt( (energy_higgs_mean+energy_higgs_dev) * (energy_higgs_mean+energy_higgs_dev) - mass_higgs * mass_higgs );
			double momentum_higgs_devDown = sqrt( (energy_higgs_mean-energy_higgs_dev) * (energy_higgs_mean-energy_higgs_dev) - mass_higgs * mass_higgs );

			// higgs -> bbar
			double seperation_bb = 2 * mass_higgs / momentum_higgs_mean; //transverse momentum???

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

		// could get the histograms right with some dummy histograms!!!
		legend->AddEntry("gr_mass_NLSP", Form("%.1f", mass_delta_vec[i_delta]), "L");

		// function to set up the graph, eg line colour and line width, the legend thing also...
		// plot styles need to have been updated first joe
		// leg_mass_NLSP->AddEntry("gr_mass_NLSP",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_mass_LSP->AddEntry("gr_mass_LSP",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_energy_quark->AddEntry("gr_energy_quark",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_momentum_quark->AddEntry("gr_momentum_quark",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_energy_NLSP->AddEntry("gr_energy_NLSP",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_momentum_NLSP->AddEntry("gr_momentum_NLSP",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_energy_LSP->AddEntry("gr_energy_LSP",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_momentum_LSP->AddEntry("gr_momentum_LSP",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_energy_higgs->AddEntry("gr_energy_higgs",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_momentum_higgs->AddEntry("gr_momentum_higgs",Form("%.2f", mass_delta_vec[i_delta]),"L");
		// leg_seperation_bb->AddEntry("gr_seperation_bb",Form("%.2f", mass_delta_vec[i_delta]),"L");

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




drawAndSave(mg_mass_LSP, legend, "mg_mass_LSP.pdf");
drawAndSave(mg_momentum_higgs, legend, "mg_momHiggs.pdf");









}