#ifndef Analysis_Analysis_boostedNmssmHiggs_PlottingDoubleBTaggerEfficiencyStudies_h
#define Analysis_Analysis_boostedNmssmHiggs_PlottingDoubleBTaggerEfficiencyStudies_h

// CPP headers
#include <string>
#include <vector>

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

// CMSSW headers

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Plotting.h"


class PlottingMcSignalStudies
{

public:

	// constructor
	PlottingMcSignalStudies(std::string,std::string);

private:

	// objects common for all plots 
	TFile * f;
	std::string outputDirectory;
	std::string massTitle;

	// individual plots
	void oneDimension_standard(std::string, std::string);
	// void oneDimension_twoPlotsAdd();
	void oneDimension_twoPlotsSeparate(std::string, std::string, std::string, std::string, std::string, double, double, double, double);
	void oneDimension_threePlotsSeparate(std::string, std::string, std::string, std::string, std::string, std::string, std::string, double, double, double, double);
	void twoDimension_standard(std::string, std::string);

	// TStyle * TDRStyle(); // now get this function from seperate file
	TStyle * tdrStyle;
	TLatex * latex = new TLatex();
};

////////////////////////////////
// class function definitions //
////////////////////////////////


//--------constructor---------//

PlottingMcSignalStudies::PlottingMcSignalStudies(std::string inputHistoFile, std::string titleName)
{
	// open input .root file containing the histograms
 	f = TFile::Open(inputHistoFile.c_str());

 	massTitle = titleName;
 	// get the name of directory holding the .root file, so we can save .pdfs here also
	for (size_t c = inputHistoFile.size()-1; c > 0; --c){
		std::string forwardSlash = "/";
		if (inputHistoFile[c] == forwardSlash[0]){
			outputDirectory = inputHistoFile.substr(0, c+1);
			break;
		}
	}
	// set the TStyle for the plots
	tdrStyle = TDRStyle();
	gROOT->SetStyle("tdrStyle");

	// set the latex defaults
	latex->SetNDC();
	latex->SetTextFont(42);
	// latex->SetTextAlign(31);

	// make the .pdfs
	oneDimension_standard("numberOfGluinos", "numberOfGluinos.pdf");

 	oneDimension_standard("detectorLeadingJet", "detectorLeadingJet.pdf");
 	oneDimension_standard("detectorSecondaryJet", "detectorSecondaryJet.pdf");
	oneDimension_standard("detectorMET", "detectorMET.pdf");
	oneDimension_standard("detectorHT", "detectorHT.pdf");

	oneDimension_twoPlotsSeparate("leadingSquarkPt", "secondarySquarkPt", "Leading Arm", "Secondary Arm", "squarkPt.pdf", 0.60, 0.88, 0.73, 0.88);
	oneDimension_twoPlotsSeparate("leadingSquarkEta", "secondarySquarkEta", "Leading Arm", "Secondary Arm", "squarkEta.pdf", 0.60, 0.88, 0.73, 0.88);
	oneDimension_threePlotsSeparate("leadingSquarkPt_zeroGluinos", "leadingSquarkPt_oneGluinos", "leadingSquarkPt_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkPt_gluLeading.pdf", 0.65, 0.88, 0.73, 0.88);
	oneDimension_threePlotsSeparate("leadingSquarkEta_zeroGluinos", "leadingSquarkEta_oneGluinos", "leadingSquarkEta_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkEta_gluLeading.pdf", 0.65, 0.88, 0.73, 0.88);
	oneDimension_threePlotsSeparate("secondarySquarkPt_zeroGluinos", "secondarySquarkPt_oneGluinos", "secondarySquarkPt_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkPt_gluSecondary.pdf", 0.65, 0.88, 0.73, 0.88);
	oneDimension_threePlotsSeparate("secondarySquarkEta_zeroGluinos", "secondarySquarkEta_oneGluinos", "secondarySquarkEta_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkEta_gluSecondary.pdf", 0.65, 0.88, 0.73, 0.88);
	twoDimension_standard("leadingSquarkPt_SecondarySquarkPt", "leadingSquarkPt_SecondarySquarkPt.pdf");
	twoDimension_standard("leadingSquarkEta_SecondarySquarkEta", "leadingSquarkEta_SecondarySquarkEta.pdf");
	twoDimension_standard("leadingSquarkPhi_SecondarySquarkPhi", "leadingSquarkPhi_SecondarySquarkPhi.pdf");

	oneDimension_twoPlotsSeparate("leadingQjetPt", "secondaryQjetPt", "Leading Arm", "Secondary Arm", "QjetPt.pdf", 0.60, 0.88, 0.73, 0.88);
	oneDimension_twoPlotsSeparate("leadingQjetEta", "secondaryQjetEta", "Leading Arm", "Secondary Arm", "QjetEta.pdf", 0.60, 0.88, 0.73, 0.88);

	twoDimension_standard("leadingQjetPt_secondaryQjetPt", "leadingQjetPt_SecondaryQjetPt.pdf");
	twoDimension_standard("leadingQjetEta_secondaryQjetEta", "leadingQjetEta_SecondaryQjetEta.pdf");
	twoDimension_standard("leadingQjetPhi_secondaryQjetPhi", "leadingQjetPhi_SecondaryQjetPhi.pdf");

	oneDimension_twoPlotsSeparate("leadingNlspPt", "secondaryNlspPt", "Leading Arm", "Secondary Arm", "nlspPt.pdf", 0.60, 0.88, 0.73, 0.88);
	oneDimension_twoPlotsSeparate("leadingNlspEta", "secondaryNlspEta", "Leading Arm", "Secondary Arm", "nlspEta.pdf", 0.60, 0.88, 0.73, 0.88);
	twoDimension_standard("leadingNlspPt_secondaryNlspPt", "leadingNlspPt_secondaryNlspPt.pdf");
	twoDimension_standard("leadingNlspEta_secondaryNlspEta", "leadingNlspEta_secondaryNlspEta.pdf");
	twoDimension_standard("leadingNlspPhi_secondaryNlspPhi", "leadingNlspPhi_secondaryNlspPhi.pdf");

	oneDimension_twoPlotsSeparate("leadingHiggsPt", "secondaryHiggsPt", "Leading Arm", "Secondary Arm", "higgsPt.pdf", 0.60, 0.88, 0.73, 0.88);
	oneDimension_twoPlotsSeparate("leadingHiggsEta", "secondaryHiggsEta", "Leading Arm", "Secondary Arm", "higgsEta.pdf", 0.60, 0.88, 0.73, 0.88);

	twoDimension_standard("leadingQjetPt_leadingHiggsPt", "leadingQjetPt_leadingHiggsPt.pdf");
	twoDimension_standard("leadingQjetEta_leadingHiggsEta", "leadingQjetEta_leadingHiggsEta.pdf");
	twoDimension_standard("leadingQjetPhi_leadingHiggsPhi", "leadingQjetPhi_leadingHiggsPhi.pdf");
	twoDimension_standard("secondaryQjetPt_secondaryHiggsPt", "secondaryQjetPt_secondaryHiggsPt.pdf");
	twoDimension_standard("secondaryQjetEta_secondaryHiggsEta", "secondaryQjetEta_secondaryHiggsEta.pdf");
	twoDimension_standard("secondaryQjetPhi_secondaryHiggsPhi", "secondaryQjetPhi_secondaryHiggsPhi.pdf");

	oneDimension_twoPlotsSeparate("leadingLspPt", "secondaryLspPt", "Leading Arm", "Secondary Arm", "lspPt.pdf", 0.60, 0.88, 0.73, 0.88);
	oneDimension_twoPlotsSeparate("leadingLspEta", "secondaryLspEta", "Leading Arm", "Secondary Arm", "lspEta.pdf", 0.60, 0.88, 0.73, 0.88);
	oneDimension_standard("lspMET", "lspMET.pdf");

	oneDimension_twoPlotsSeparate("leadingBBbarSeperation", "secondaryBBbarSeperation", "Leading Arm", "Secondary Arm", "BBbarSeperation.pdf", 0.60, 0.88, 0.73, 0.88);
	oneDimension_twoPlotsSeparate("leadingBBbarInvmass", "secondaryBBbarInvmass", "Leading Arm", "Secondary Arm", "BBbarInvmass.pdf", 0.60, 0.88, 0.73, 0.88);
	twoDimension_standard("leadingBBbarSeperation_massHiggsOverPt", "leadingBBbarSeperation_massHiggsOverPt.pdf");
	twoDimension_standard("secondaryBBbarSeperation_massHiggsOverPt", "secondaryBBbarSeperation_massHiggsOverPt.pdf");

// ADD TO WORKFLOW
// DPHI
// DR DR SCATTER

}

//-----------public-----------//


//-----------private----------//

void PlottingMcSignalStudies::oneDimension_standard(std::string histoname, std::string saveName)
{
    TCanvas* c=new TCanvas("c","c"); 	
	TH1F * h = (TH1F*)f->Get(Form("%s", histoname.c_str()));
	h->SetLineWidth(2);
	// h->SetLineColor(2);
	// h->GetXaxis()->SetTitle("");
	h->GetXaxis()->SetTitleSize(0.06);	
	h->GetXaxis()->SetLabelSize(0.05);
	// h->GetYaxis()->SetTitle("");
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.05);
	h->Draw();
	
	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
}





void PlottingMcSignalStudies::oneDimension_twoPlotsSeparate(std::string histoname1, std::string histoname2, std::string legendName1, std::string legendName2, std::string saveName, double legXmin, double legXmax, double legYmin, double legYmax)
{
    TCanvas* c=new TCanvas("c","c"); 	
	TH1F * h1 = (TH1F*)f->Get(Form("%s", histoname1.c_str()));
	TH1F * h2 = (TH1F*)f->Get(Form("%s", histoname2.c_str()));	
	
	h1->SetLineWidth(2);
	// h1->SetLineColor(2);
	// h1->GetXaxis()->SetTitle("");
	h1->GetXaxis()->SetTitleSize(0.06);	
	h1->GetXaxis()->SetLabelSize(0.05);
	// h1->GetYaxis()->SetTitle("");
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->Draw();
	
	h2->SetLineWidth(2);
	h2->SetLineColor(2);
	// h2->GetXaxis()->SetTitle("");
	h2->GetXaxis()->SetTitleSize(0.06);	
	h2->GetXaxis()->SetLabelSize(0.05);
	// h2->GetYaxis()->SetTitle("");
	h2->GetYaxis()->SetTitleSize(0.06);
	h2->GetYaxis()->SetLabelSize(0.05);
	h2->Draw("same");

	TLegend * legend = new TLegend();
	legend->SetX1NDC(legXmin);
	legend->SetX2NDC(legXmax);
	legend->SetY1NDC(legYmin);
	legend->SetY2NDC(legYmax);
	legend->AddEntry(h1, legendName1.c_str(), "L");
	legend->AddEntry(h2, legendName2.c_str(), "L");
	legend->Draw("same");

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
}





void PlottingMcSignalStudies::oneDimension_threePlotsSeparate(std::string histoname1, std::string histoname2, std::string histoname3, std::string legendName1, std::string legendName2, std::string legendName3, std::string saveName, double legXmin, double legXmax, double legYmin, double legYmax)
{
    TCanvas* c=new TCanvas("c","c"); 	
	TH1F * h1 = (TH1F*)f->Get(Form("%s", histoname1.c_str()));
	TH1F * h2 = (TH1F*)f->Get(Form("%s", histoname2.c_str()));	
	TH1F * h3 = (TH1F*)f->Get(Form("%s", histoname3.c_str()));	
	Double_t norm = h1->GetEntries();
	h1->Scale(1/norm);
	norm = h2->GetEntries();
	h2->Scale(1/norm);
	norm = h3->GetEntries();
	h3->Scale(1/norm);

	h3->SetLineWidth(2);
	h3->SetLineColor(kGreen+1);
	// h3->GetXaxis()->SetTitle("");
	h3->GetXaxis()->SetTitleSize(0.06);	
	h3->GetXaxis()->SetLabelSize(0.05);
	// h3->GetYaxis()->SetTitle("");
	h3->GetYaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetLabelSize(0.05);
	h3->Draw();
	
	h2->SetLineWidth(2);
	h2->SetLineColor(2);
	// h2->GetXaxis()->SetTitle("");
	h2->GetXaxis()->SetTitleSize(0.06);	
	h2->GetXaxis()->SetLabelSize(0.05);
	// h2->GetYaxis()->SetTitle("");
	h2->GetYaxis()->SetTitleSize(0.06);
	h2->GetYaxis()->SetLabelSize(0.05);
	h2->Draw("same");

	h1->SetLineWidth(2);
	// h1->SetLineColor(2);
	// h1->GetXaxis()->SetTitle("");
	h1->GetXaxis()->SetTitleSize(0.06);	
	h1->GetXaxis()->SetLabelSize(0.05);
	// h1->GetYaxis()->SetTitle("");
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->Draw("same");

	TLegend * legend = new TLegend();
	legend->SetX1NDC(legXmin);
	legend->SetX2NDC(legXmax);
	legend->SetY1NDC(legYmin);
	legend->SetY2NDC(legYmax);
	legend->AddEntry(h1, legendName1.c_str(), "L");
	legend->AddEntry(h2, legendName2.c_str(), "L");
	legend->AddEntry(h3, legendName3.c_str(), "L");
	legend->Draw("same");

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
}





void PlottingMcSignalStudies::twoDimension_standard(std::string histoname, std::string saveName)
{
    TCanvas* c=new TCanvas("c","c"); 	
	TH2F * h = (TH2F*)f->Get(Form("%s", histoname.c_str()));
	h->SetLineWidth(2);
	// h->SetLineColor(2);
	// h->GetXaxis()->SetTitle("");
	h->GetXaxis()->SetTitleSize(0.06);	
	h->GetXaxis()->SetLabelSize(0.05);
	// h->GetYaxis()->SetTitle("");
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.05);
	h->Draw("colz");
	
	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
}

#endif