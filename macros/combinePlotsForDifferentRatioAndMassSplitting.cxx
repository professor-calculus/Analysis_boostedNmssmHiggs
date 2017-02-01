// CPP headers
#include <string>
#include <vector>
#include <iostream>

// ROOT headers
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Plotting.h"

// script to submit McSignalStudies many times
// run with
// $ root -q -b -l $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/macros/combinePlotsForDifferentRatioAndMassSplitting.cxx

void Plotter_standard(std::string motherDir, std::string outputDir, std::string plotTitle, std::vector<std::string> inputDirs, std::vector<std::string> legendNames, std::string histogramName, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool normalPlotOrdering);
void Plotter_log(std::string motherDir, std::string outputDir, std::string plotTitle, std::vector<std::string> inputDirs, std::vector<std::string> legendNames, std::string histogramName, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool normalPlotOrdering);
void Plotter_add(std::string motherDir, std::string outputDir, std::string plotTitle, std::vector<std::string> inputDirs, std::vector<std::string> legendNames, std::string histogramName, std::string histogramName2, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool normalPlotOrdering);

void combinePlotsForDifferentRatioAndMassSplitting()
{
	std::string motherDir = "/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/paramCardType_03/version_varyRatioAndSplitting/"; // needs the forward slash on the end
	std::string outputDir = motherDir + "combinedWithExtraRatioKEEP3/";  // needs the forward slash on the end

	std::vector<std::string> vecInputDirSplitting = {"mH70p0_mSusy1200p0_ratio0p8_splitting0p1_25000events", "mH70p0_mSusy1200p0_ratio0p8_splitting5p0_25000events", "mH70p0_mSusy1200p0_ratio0p8_splitting15p0_25000events"};
	std::string commonTitleSplitting = "mH70p0_mSusy1200p0_ratio0p8_25000events";
	std::vector<std::string> vecLegendSplitting = {"Mass Splitting" ,"0.1 (GeV)", "5.0 (GeV)", "15.0 (GeV)"}; // legend title then should follow vecInputDirSplitting entries

	std::vector<std::string> vecInputDirRatio = {"mH70p0_mSusy1200p0_ratio0p8_splitting0p1_25000events", "mH70p0_mSusy1200p0_ratio0p9_splitting0p1_25000events", "mH70p0_mSusy1200p0_ratio0p95_splitting0p1_25000events", "mH70p0_mSusy1200p0_ratio0p99_splitting0p1_25000events"};
	std::string commonTitleRatio= "mH70p0_mSusy1200p0_splitting0p1_25000events";
	std::vector<std::string> vecLegendRatio = {"Mass Ratio", "0.80", "0.90", "0.95", "0.99"}; // legend title then should follow vecInputDirRatio entries

	bool makeDir = !(std::system(Form("mkdir %s",outputDir.c_str())));
	if (makeDir == false){
		std::cout << "The chosen output directory already exists, or the parent directory does not exist:" << std::endl;
		std::cout << "Do not wish to overwrite: Exiting y'all..." << std::endl;
		return;
	}

	Plotter_standard(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "lspMET", "constantRatioVarySplitting_lspMET.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_log(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "lspMET", "constantRatioVarySplitting_lspMETlog.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "detectorMET", "constantRatioVarySplitting_detectorMET.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "detectorMHT", "constantRatioVarySplitting_detectorMHT.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "detectorHT", "constantRatioVarySplitting_detectorHT.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "detectorLeadingAk4JetPt", "constantRatioVarySplitting_detectorLeadingAk4JetPt.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "detectorSecondaryAk4JetPt", "constantRatioVarySplitting_detectorSecondaryAk4JetPt.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "detectorLeadingAk4JetEta", "constantRatioVarySplitting_detectorLeadingAk4JetEta.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "detectorSecondaryAk4JetEta", "constantRatioVarySplitting_detectorSecondaryAk4JetEta.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_add(motherDir, outputDir, commonTitleSplitting, vecInputDirSplitting, vecLegendSplitting, "leadingBBbarSeperation", "secondaryBBbarSeperation", "constantRatioVarySplitting_bBbarSeperationBothArms.pdf", 0.69, 0.88, 0.67, 0.87, true);

	Plotter_standard(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "lspMET", "constantSplittingVaryRatio_lspMET.pdf", 0.69, 0.88, 0.67, 0.87, false);
	Plotter_log(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "lspMET", "constantSplittingVaryRatio_lspMETlog.pdf", 0.69, 0.88, 0.67, 0.87, false);
	Plotter_standard(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "detectorMET", "constantSplittingVaryRatio_detectorMET.pdf", 0.69, 0.88, 0.67, 0.87, false);
	Plotter_standard(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "detectorMHT", "constantSplittingVaryRatio_detectorMHT.pdf", 0.69, 0.88, 0.67, 0.87, false);
	Plotter_standard(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "detectorHT", "constantSplittingVaryRatio_detectorHT.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "detectorLeadingAk4JetPt", "constantSplittingVaryRatio_detectorLeadingAk4JetPt.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "detectorSecondaryAk4JetPt", "constantSplittingVaryRatio_detectorSecondaryAk4JetPt.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "detectorLeadingAk4JetEta", "constantSplittingVaryRatio_detectorLeadingAk4JetEta.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_standard(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "detectorSecondaryAk4JetEta", "constantSplittingVaryRatio_detectorSecondaryAk4JetEta.pdf", 0.69, 0.88, 0.67, 0.87, true);
	Plotter_add(motherDir, outputDir, commonTitleRatio, vecInputDirRatio, vecLegendRatio, "leadingBBbarSeperation", "secondaryBBbarSeperation", "constantSplittingVaryRatio_bBbarSeperationBothArms.pdf", 0.69, 0.88, 0.67, 0.87, false);

} // closes function "combine_McSignalStudies"






void Plotter_standard(std::string motherDir, std::string outputDir, std::string plotTitle, std::vector<std::string> inputDirs, std::vector<std::string> legendNames, std::string histogramName, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool normalPlotOrdering)
{
	// set the TStyle for the plots
	TStyle * tdrStyle = TDRStyle();
	gROOT->SetStyle("tdrStyle");
	// set the latex defaults
	TLatex * latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);

	TCanvas* c=new TCanvas("c","c"); 

	TLegend * legend = new TLegend();
	legend->SetHeader(legendNames[0].c_str());
	legend->SetX1NDC(legXmin);
	legend->SetX2NDC(legXmax);
	legend->SetY1NDC(legYmin);
	legend->SetY2NDC(legYmax);

	if (normalPlotOrdering){
		for (size_t iD = 0; iD<inputDirs.size(); ++iD){
		// for (size_t iD = 0; iD<1; ++iD){

			std::string inputHistoFile = motherDir + inputDirs[iD] + "/output.root";
			TFile * f = TFile::Open(inputHistoFile.c_str());
			TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
			h->SetLineWidth(2);
			// h->SetLineColor(SetColor(iD, inputDirs.size() + 1));
			h->SetLineColor(SetColor(iD, inputDirs.size()));
			// h->GetXaxis()->SetTitle("");
			h->GetXaxis()->SetTitleSize(0.06);	
			h->GetXaxis()->SetLabelSize(0.05);
			// h->GetYaxis()->SetTitle("");
			h->GetYaxis()->SetTitleSize(0.06);
			h->GetYaxis()->SetLabelSize(0.05);
			h->Draw("same");
			legend->AddEntry(h, legendNames[iD+1].c_str(), "L");		
		}
	}
	else{
		for (size_t iD=inputDirs.size()-1; iD != -1; --iD){
		// for (size_t iD=0; iD != -1; --iD){

			std::string inputHistoFile = motherDir + inputDirs[iD] + "/output.root";
			TFile * f = TFile::Open(inputHistoFile.c_str());
			TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
			h->SetLineWidth(2);
			// h->SetLineColor(SetColor(iD, inputDirs.size() + 1));
			h->SetLineColor(SetColor(iD, inputDirs.size()));
			// h->GetXaxis()->SetTitle("");
			h->GetXaxis()->SetTitleSize(0.06);	
			h->GetXaxis()->SetLabelSize(0.05);
			// h->GetYaxis()->SetTitle("");
			h->GetYaxis()->SetTitleSize(0.06);
			h->GetYaxis()->SetLabelSize(0.05);
			h->Draw("same");
			legend->AddEntry(h, legendNames[iD+1].c_str(), "L");
		}
	}

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,plotTitle.c_str());
	legend->Draw("same");
	std::string saveName = outputDir + baseSaveName;
	c->SaveAs(saveName.c_str());
	c->Close();
} // closes function "Plotter_standard"







void Plotter_log(std::string motherDir, std::string outputDir, std::string plotTitle, std::vector<std::string> inputDirs, std::vector<std::string> legendNames, std::string histogramName, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool normalPlotOrdering)
{
	// set the TStyle for the plots
	TStyle * tdrStyle = TDRStyle();
	gROOT->SetStyle("tdrStyle");
	// set the latex defaults
	TLatex * latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);

	TCanvas* c=new TCanvas("c","c"); 
	gPad->SetLogy();
	TLegend * legend = new TLegend();
	legend->SetHeader(legendNames[0].c_str());
	legend->SetX1NDC(legXmin);
	legend->SetX2NDC(legXmax);
	legend->SetY1NDC(legYmin);
	legend->SetY2NDC(legYmax);

	if (normalPlotOrdering){
		for (size_t iD=inputDirs.size()-1; iD != -1; --iD){
		// for (size_t iD=0; iD != -1; --iD){

			std::string inputHistoFile = motherDir + inputDirs[iD] + "/output.root";
			TFile * f = TFile::Open(inputHistoFile.c_str());
			TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
			h->SetLineWidth(2);
			// h->SetLineColor(SetColor(iD, inputDirs.size() + 1));
			h->SetLineColor(SetColor(iD, inputDirs.size()));
			// h->GetXaxis()->SetTitle("");
			h->GetXaxis()->SetTitleSize(0.06);	
			h->GetXaxis()->SetLabelSize(0.05);
			// h->GetYaxis()->SetTitle("");
			h->GetYaxis()->SetTitleSize(0.06);
			h->GetYaxis()->SetLabelSize(0.05);
			h->Draw("same");
			legend->AddEntry(h, legendNames[iD+1].c_str(), "L");		
		}
	}
	else{
		for (size_t iD=inputDirs.size()-1; iD != -1; --iD){
		// for (size_t iD=0; iD != -1; --iD){

			std::string inputHistoFile = motherDir + inputDirs[iD] + "/output.root";
			TFile * f = TFile::Open(inputHistoFile.c_str());
			TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
			h->SetLineWidth(2);
			// h->SetLineColor(SetColor(iD, inputDirs.size() + 1));
			h->SetLineColor(SetColor(iD, inputDirs.size()));
			// h->GetXaxis()->SetTitle("");
			h->GetXaxis()->SetTitleSize(0.06);	
			h->GetXaxis()->SetLabelSize(0.05);
			// h->GetYaxis()->SetTitle("");
			h->GetYaxis()->SetTitleSize(0.06);
			h->GetYaxis()->SetLabelSize(0.05);
			h->Draw("same");
			legend->AddEntry(h, legendNames[iD+1].c_str(), "L");
		}
	}

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,plotTitle.c_str());
	legend->Draw("same");
	std::string saveName = outputDir + baseSaveName;
	c->SaveAs(saveName.c_str());
	c->Close();
} // closes function "Plotter_standard"






void Plotter_add(std::string motherDir, std::string outputDir, std::string plotTitle, std::vector<std::string> inputDirs, std::vector<std::string> legendNames, std::string histogramName, std::string histogramName2, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool normalPlotOrdering)
{
	// set the TStyle for the plots
	TStyle * tdrStyle = TDRStyle();
	gROOT->SetStyle("tdrStyle");
	// set the latex defaults
	TLatex * latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);

	TCanvas* c=new TCanvas("c","c"); 

	TLegend * legend = new TLegend();
	legend->SetHeader(legendNames[0].c_str());
	legend->SetX1NDC(legXmin);
	legend->SetX2NDC(legXmax);
	legend->SetY1NDC(legYmin);
	legend->SetY2NDC(legYmax);

	if (normalPlotOrdering){
		for (size_t iD = 0; iD<inputDirs.size(); ++iD){
		// for (size_t iD = 0; iD<1; ++iD){

			std::string inputHistoFile = motherDir + inputDirs[iD] + "/output.root";
			TFile * f = TFile::Open(inputHistoFile.c_str());
			TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
			TH1F * hDummy = (TH1F*)f->Get(Form("%s", histogramName2.c_str()));
			h->Add(hDummy);
			h->SetLineWidth(2);
			// h->SetLineColor(SetColor(iD, inputDirs.size() + 1));
			h->SetLineColor(SetColor(iD, inputDirs.size()));
			// h->GetXaxis()->SetTitle("");
			h->GetXaxis()->SetTitleSize(0.06);	
			h->GetXaxis()->SetLabelSize(0.05);
			// h->GetYaxis()->SetTitle("");
			h->GetYaxis()->SetTitleSize(0.06);
			h->GetYaxis()->SetLabelSize(0.05);
			h->Draw("same");
			legend->AddEntry(h, legendNames[iD+1].c_str(), "L");		
		}
	}
	else{
		for (size_t iD=inputDirs.size()-1; iD != -1; --iD){

			std::string inputHistoFile = motherDir + inputDirs[iD] + "/output.root";
			TFile * f = TFile::Open(inputHistoFile.c_str());
			TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
			TH1F * hDummy = (TH1F*)f->Get(Form("%s", histogramName2.c_str()));
			h->Add(hDummy);
			h->SetLineWidth(2);
			// h->SetLineColor(SetColor(iD, inputDirs.size() + 1));
			h->SetLineColor(SetColor(iD, inputDirs.size()));
			// h->GetXaxis()->SetTitle("");
			h->GetXaxis()->SetTitleSize(0.06);	
			h->GetXaxis()->SetLabelSize(0.05);
			// h->GetYaxis()->SetTitle("");
			h->GetYaxis()->SetTitleSize(0.06);
			h->GetYaxis()->SetLabelSize(0.05);
			h->Draw("same");
			legend->AddEntry(h, legendNames[iD+1].c_str(), "L");
		}
	}

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,plotTitle.c_str());
	legend->Draw("same");
	std::string saveName = outputDir + baseSaveName;
	c->SaveAs(saveName.c_str());
	c->Close();
} // closes function "Plotter_standard"
