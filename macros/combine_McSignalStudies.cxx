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
// $ root -q -b -l combine_McSignalStudies.cxx

void Plotter_standard(std::string motherDir, std::string outputDir, std::vector<std::string> vecHiggsMassToUse, std::vector<std::string> vecSusyMassToUse, std::string tailArgument, std::string histogramName, std::string baseSaveName, double, double, double, double, bool, bool);
void Plotter_add(std::string motherDir, std::string outputDir, std::vector<std::string> vecHiggsMassToUse, std::vector<std::string> vecSusyMassToUse, std::string tailArgument, std::string histogramName1, std::string histogramName2, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool, bool);

void combine_McSignalStudies()
{
	// * U S E R * I N P U T * //
	std::string motherDir = "/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/"; // needs the forward slash on the end
	std::vector<std::string> vecHiggsMassToUse = {"70p0", "90p0", "110p0", "130p0"}; // has to match my already defined naming convention
	std::vector<std::string> vecSusyMassToUse = {"1000p0", "1200p0", "1600p0", "2000p0"}; // has to match my already defined naming convention
	std::string tailArgument = "_ratio0p95_splitting1p0_10000events"; // needs the dash at the beginning
	std::string outputDir = motherDir + "combinedPlots_v1/";  // needs the forward slash on the end
	// * * * * * * * * * * * * //

	bool makeDir = !(std::system(Form("mkdir %s",outputDir.c_str())));
	if (makeDir == false){
		std::cout << "The chosen output directory already exists, or the parent directory does not exist:" << std::endl;
		std::cout << "Do not wish to overwrite: Exiting y'all..." << std::endl;
		return;
	}


	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "lspMET", "lspMET", 0.69, 0.88, 0.67, 0.87, true, false);
	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "detectorHT", "detectorHT", 0.69, 0.88, 0.67, 0.87, true, true);
	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "detectorLeadingJet", "detectorLeadingJet", 0.69, 0.88, 0.67, 0.87, true, true);
	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "detectorSecondaryJet", "detectorSecondaryJet", 0.69, 0.88, 0.67, 0.87, true, true);
	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "detectorMET", "detectorMET", 0.69, 0.88, 0.67, 0.87, true, false);
	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "numberOfGluinos", "numberOfGluinos", 0.69, 0.88, 0.67, 0.87, false, true);
	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "numberOfDetectorJets", "numberOfDetectorJets", 0.69, 0.88, 0.67, 0.87, false, true);
	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "leadingHiggsQjetDphi", "leadingHiggsQjetDphi", 0.69, 0.88, 0.67, 0.87, false, true);
	Plotter_standard(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "secondaryHiggsQjetDphi", "secondaryHiggsQjetDphi", 0.69, 0.88, 0.67, 0.87, false, true);
	Plotter_add(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "leadingQjetPt", "secondaryQjetPt", "qjetPtBothArms", 0.69, 0.88, 0.67, 0.87, true, true);
	Plotter_add(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "leadingQjetEta", "secondaryQjetEta", "qjetEtaBothArms", 0.69, 0.88, 0.67, 0.87, false, true);
	Plotter_add(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "leadingSquarkPt", "secondarySquarkPt", "squarkPtBothArms", 0.69, 0.88, 0.67, 0.87, true, true);
	Plotter_add(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "leadingSquarkEta", "secondarySquarkEta", "squarkEtaBothArms", 0.69, 0.88, 0.67, 0.87, false, true);
	Plotter_add(motherDir, outputDir, vecHiggsMassToUse, vecSusyMassToUse, tailArgument, "leadingBBbarSeperation", "secondaryBBbarSeperation", "BBbarSeperationBothArms", 0.69, 0.88, 0.67, 0.87, false, true);

} // closes function "combine_McSignalStudies"






void Plotter_standard(std::string motherDir, std::string outputDir, std::vector<std::string> vecHiggsMassToUse, std::vector<std::string> vecSusyMassToUse, std::string tailArgument, std::string histogramName, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool normalLoopOrderSusy, bool normalLoopOrderHiggs)
{

	// set the TStyle for the plots
	TStyle * tdrStyle = TDRStyle();
	gROOT->SetStyle("tdrStyle");
	// set the latex defaults
	TLatex * latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);

	// Fix Higgs mass, change Susy mass
	for (size_t iH=0; iH<vecHiggsMassToUse.size(); ++iH){

		TCanvas* c=new TCanvas("c","c"); 	
		TLegend * legend = new TLegend();
		legend->SetHeader("Mass SUSY");
		legend->SetX1NDC(legXmin);
		legend->SetX2NDC(legXmax);
		legend->SetY1NDC(legYmin);
		legend->SetY2NDC(legYmax);
		std::string saveName = outputDir + baseSaveName + "_fixMassHiggs" + vecHiggsMassToUse[iH] + "_varyMassSusy.pdf";
		std::string plotTitle = "mH" + vecHiggsMassToUse[iH] + tailArgument;

		if (normalLoopOrderSusy){
			for (size_t iSu=0; iSu<vecSusyMassToUse.size(); ++iSu){

				std::string inputHistoFile = motherDir + "mH" + vecHiggsMassToUse[iH] + "_mSusy" + vecSusyMassToUse[iSu] + tailArgument + "/output.root";
				TFile * f = TFile::Open(inputHistoFile.c_str());
				TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
				h->SetLineWidth(2);
				h->SetLineColor(SetColor(iSu, vecSusyMassToUse.size()));
				// h->GetXaxis()->SetTitle("");
				h->GetXaxis()->SetTitleSize(0.06);	
				h->GetXaxis()->SetLabelSize(0.05);
				// h->GetYaxis()->SetTitle("");
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetLabelSize(0.05);
				h->Draw("same");
				std::string legendName = vecSusyMassToUse[iSu].substr(0,vecSusyMassToUse[iSu].size()-2);
				legendName = legendName + " (GeV)";
				legend->AddEntry(h, legendName.c_str(), "L");
			}
		}
		else{
			for (size_t iSu=vecSusyMassToUse.size()-1; iSu != -1; --iSu){

				std::string inputHistoFile = motherDir + "mH" + vecHiggsMassToUse[iH] + "_mSusy" + vecSusyMassToUse[iSu] + tailArgument + "/output.root";
				TFile * f = TFile::Open(inputHistoFile.c_str());
				TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
				h->SetLineWidth(2);
				h->SetLineColor(SetColor(iSu, vecSusyMassToUse.size()));
				// h->GetXaxis()->SetTitle("");
				h->GetXaxis()->SetTitleSize(0.06);	
				h->GetXaxis()->SetLabelSize(0.05);
				// h->GetYaxis()->SetTitle("");
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetLabelSize(0.05);
				h->Draw("same");
				std::string legendName = vecSusyMassToUse[iSu].substr(0,vecSusyMassToUse[iSu].size()-2);
				legendName = legendName + " (GeV)";
				legend->AddEntry(h, legendName.c_str(), "L");
			}
		}

		latex->SetTextAlign(11); // align from left
		latex->DrawLatex(0.15,0.92,plotTitle.c_str());
		legend->Draw("same");
		c->SaveAs(saveName.c_str());
		c->Close();
	}

	// Fix Susy mass, change Higgs mass
	for (size_t iSu=0; iSu<vecSusyMassToUse.size(); ++iSu){

		TCanvas* c=new TCanvas("c","c"); 	
		TLegend * legend = new TLegend();
		legend->SetHeader("Mass Higgs");
		legend->SetX1NDC(legXmin);
		legend->SetX2NDC(legXmax);
		legend->SetY1NDC(legYmin);
		legend->SetY2NDC(legYmax);
		std::string saveName = outputDir + baseSaveName + "_fixSusyMass" + vecSusyMassToUse[iSu] + "_varyMassHiggs.pdf";
		std::string plotTitle = "mSusy" + vecSusyMassToUse[iSu] + tailArgument;

		if (normalLoopOrderHiggs){
			for (size_t iH=0; iH<vecHiggsMassToUse.size(); ++iH){

				std::string inputHistoFile = motherDir + "mH" + vecHiggsMassToUse[iH] + "_mSusy" + vecSusyMassToUse[iSu] + tailArgument + "/output.root";
				TFile * f = TFile::Open(inputHistoFile.c_str());
				TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
				h->SetLineWidth(2);
				h->SetLineColor(SetColor(iH, vecSusyMassToUse.size()));
				// h->GetXaxis()->SetTitle("");
				h->GetXaxis()->SetTitleSize(0.06);	
				h->GetXaxis()->SetLabelSize(0.05);
				// h->GetYaxis()->SetTitle("");
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetLabelSize(0.05);
				h->Draw("same");
				std::string legendName = vecHiggsMassToUse[iH].substr(0,vecHiggsMassToUse[iH].size()-2);
				legendName = legendName + " (GeV)";
				legend->AddEntry(h, legendName.c_str(), "L");
			}
		}
		else{
			for (size_t iH=vecHiggsMassToUse.size()-1; iH != -1; --iH){

				std::string inputHistoFile = motherDir + "mH" + vecHiggsMassToUse[iH] + "_mSusy" + vecSusyMassToUse[iSu] + tailArgument + "/output.root";
				TFile * f = TFile::Open(inputHistoFile.c_str());
				TH1F * h = (TH1F*)f->Get(Form("%s", histogramName.c_str()));
				h->SetLineWidth(2);
				h->SetLineColor(SetColor(iH, vecSusyMassToUse.size()));
				// h->GetXaxis()->SetTitle("");
				h->GetXaxis()->SetTitleSize(0.06);	
				h->GetXaxis()->SetLabelSize(0.05);
				// h->GetYaxis()->SetTitle("");
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetLabelSize(0.05);
				h->Draw("same");
				std::string legendName = vecHiggsMassToUse[iH].substr(0,vecHiggsMassToUse[iH].size()-2);
				legendName = legendName + " (GeV)";
				legend->AddEntry(h, legendName.c_str(), "L");
			}		
		}

		latex->SetTextAlign(11); // align from left
		latex->DrawLatex(0.15,0.92,plotTitle.c_str());
		legend->Draw("same");
		c->SaveAs(saveName.c_str());
		c->Close();
	}

} // closes function "Plotter_standard"








void Plotter_add(std::string motherDir, std::string outputDir, std::vector<std::string> vecHiggsMassToUse, std::vector<std::string> vecSusyMassToUse, std::string tailArgument, std::string histogramName1, std::string histogramName2, std::string baseSaveName, double legXmin, double legXmax, double legYmin, double legYmax, bool normalLoopOrderSusy, bool normalLoopOrderHiggs)
{

	// set the TStyle for the plots
	TStyle * tdrStyle = TDRStyle();
	gROOT->SetStyle("tdrStyle");
	// set the latex defaults
	TLatex * latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);

	// Fix Higgs mass, change Susy mass
	for (size_t iH=0; iH<vecHiggsMassToUse.size(); ++iH){

		TCanvas* c=new TCanvas("c","c"); 	
		TLegend * legend = new TLegend();
		legend->SetHeader("Mass SUSY");
		legend->SetX1NDC(legXmin);
		legend->SetX2NDC(legXmax);
		legend->SetY1NDC(legYmin);
		legend->SetY2NDC(legYmax);
		std::string saveName = outputDir + baseSaveName + "_fixMassHiggs" + vecHiggsMassToUse[iH] + "_varyMassSusy.pdf";
		std::string plotTitle = "mH" + vecHiggsMassToUse[iH] + tailArgument;

		if (normalLoopOrderSusy){
			for (size_t iSu=0; iSu<vecSusyMassToUse.size(); ++iSu){

				std::string inputHistoFile = motherDir + "mH" + vecHiggsMassToUse[iH] + "_mSusy" + vecSusyMassToUse[iSu] + tailArgument + "/output.root";
				TFile * f = TFile::Open(inputHistoFile.c_str());
				TH1F * h = (TH1F*)f->Get(Form("%s", histogramName1.c_str()));
				TH1F * hDummy = (TH1F*)f->Get(Form("%s", histogramName2.c_str()));
				h->Add(hDummy);
				h->SetLineWidth(2);
				h->SetLineColor(SetColor(iSu, vecSusyMassToUse.size()));
				// h->GetXaxis()->SetTitle("");
				h->GetXaxis()->SetTitleSize(0.06);	
				h->GetXaxis()->SetLabelSize(0.05);
				// h->GetYaxis()->SetTitle("");
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetLabelSize(0.05);
				h->Draw("same");
				std::string legendName = vecSusyMassToUse[iSu].substr(0,vecSusyMassToUse[iSu].size()-2);
				legendName = legendName + " (GeV)";
				legend->AddEntry(h, legendName.c_str(), "L");
			}
		}
		else{
			for (size_t iSu = vecSusyMassToUse.size()-1; iSu != -1; --iSu){

				std::string inputHistoFile = motherDir + "mH" + vecHiggsMassToUse[iH] + "_mSusy" + vecSusyMassToUse[iSu] + tailArgument + "/output.root";
				TFile * f = TFile::Open(inputHistoFile.c_str());
				TH1F * h = (TH1F*)f->Get(Form("%s", histogramName1.c_str()));
				TH1F * hDummy = (TH1F*)f->Get(Form("%s", histogramName2.c_str()));
				h->Add(hDummy);
				h->SetLineWidth(2);
				h->SetLineColor(SetColor(iSu, vecSusyMassToUse.size()));
				// h->GetXaxis()->SetTitle("");
				h->GetXaxis()->SetTitleSize(0.06);	
				h->GetXaxis()->SetLabelSize(0.05);
				// h->GetYaxis()->SetTitle("");
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetLabelSize(0.05);
				h->Draw("same");
				std::string legendName = vecSusyMassToUse[iSu].substr(0,vecSusyMassToUse[iSu].size()-2);
				legendName = legendName + " (GeV)";
				legend->AddEntry(h, legendName.c_str(), "L");
			}			
		}	
		latex->SetTextAlign(11); // align from left
		latex->DrawLatex(0.15,0.92,plotTitle.c_str());
		legend->Draw("same");
		c->SaveAs(saveName.c_str());
		c->Close();
	}

	// Fix Susy mass, change Higgs mass
	for (size_t iSu=0; iSu<vecSusyMassToUse.size(); ++iSu){

		TCanvas* c=new TCanvas("c","c"); 	
		TLegend * legend = new TLegend();
		legend->SetHeader("Mass Higgs");
		legend->SetX1NDC(legXmin);
		legend->SetX2NDC(legXmax);
		legend->SetY1NDC(legYmin);
		legend->SetY2NDC(legYmax);
		std::string saveName = outputDir + baseSaveName + "_fixSusyMass" + vecSusyMassToUse[iSu] + "_varyMassHiggs.pdf";
		std::string plotTitle = "mSusy" + vecSusyMassToUse[iSu] + tailArgument;

		if (normalLoopOrderHiggs){
			for (size_t iH=0; iH<vecHiggsMassToUse.size(); ++iH){

				std::string inputHistoFile = motherDir + "mH" + vecHiggsMassToUse[iH] + "_mSusy" + vecSusyMassToUse[iSu] + tailArgument + "/output.root";
				TFile * f = TFile::Open(inputHistoFile.c_str());
				TH1F * h = (TH1F*)f->Get(Form("%s", histogramName1.c_str()));
				TH1F * hDummy = (TH1F*)f->Get(Form("%s", histogramName2.c_str()));
				h->Add(hDummy);
				h->SetLineWidth(2);
				h->SetLineColor(SetColor(iH, vecSusyMassToUse.size()));
				// h->GetXaxis()->SetTitle("");
				h->GetXaxis()->SetTitleSize(0.06);	
				h->GetXaxis()->SetLabelSize(0.05);
				// h->GetYaxis()->SetTitle("");
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetLabelSize(0.05);
				h->Draw("same");
				std::string legendName = vecHiggsMassToUse[iH].substr(0,vecHiggsMassToUse[iH].size()-2);
				legendName = legendName + " (GeV)";
				legend->AddEntry(h, legendName.c_str(), "L");
			}
		}
		else{
			for (size_t iH=vecHiggsMassToUse.size()-1; iH != -1; --iH){

				std::string inputHistoFile = motherDir + "mH" + vecHiggsMassToUse[iH] + "_mSusy" + vecSusyMassToUse[iSu] + tailArgument + "/output.root";
				TFile * f = TFile::Open(inputHistoFile.c_str());
				TH1F * h = (TH1F*)f->Get(Form("%s", histogramName1.c_str()));
				TH1F * hDummy = (TH1F*)f->Get(Form("%s", histogramName2.c_str()));
				h->Add(hDummy);
				h->SetLineWidth(2);
				h->SetLineColor(SetColor(iH, vecHiggsMassToUse.size()));
				// h->GetXaxis()->SetTitle("");
				h->GetXaxis()->SetTitleSize(0.06);	
				h->GetXaxis()->SetLabelSize(0.05);
				// h->GetYaxis()->SetTitle("");
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetLabelSize(0.05);
				h->Draw("same");
				std::string legendName = vecHiggsMassToUse[iH].substr(0,vecHiggsMassToUse[iH].size()-2);
				legendName = legendName + " (GeV)";
				legend->AddEntry(h, legendName.c_str(), "L");
			}
		}

		latex->SetTextAlign(11); // align from left
		latex->DrawLatex(0.15,0.92,plotTitle.c_str());
		legend->Draw("same");
		c->SaveAs(saveName.c_str());
		c->Close();
	}

} // closes function "Plotter_standard"
