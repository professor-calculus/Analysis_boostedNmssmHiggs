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
	void oneDimension_twoPlotsSeparate(std::string, std::string, std::string, std::string, std::string, double, double, double, double, int);
	void oneDimension_threePlotsSeparate(std::string, std::string, std::string, std::string, std::string, std::string, std::string, double, double, double, double);
	void twoDimension_standard(std::string, std::string);
	void twoDimension_addPlots(std::string, std::string, std::string);
	void twoDimension_addPlots_withGradientLines(std::string, std::string, std::string);


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
	// legXmin, legXmax, legYmin, legYmax, whichHistoFirst?
 	oneDimension_twoPlotsSeparate("detectorLeadingJet", "detectorSecondaryJet", "Leading", "Secondary", "detectorLeadingJet.pdf",  0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_standard("detectorMET", "detectorMET.pdf");
	oneDimension_standard("detectorHT", "detectorHT.pdf");
	twoDimension_standard("detectorSecondaryJet_detectorLeadingJet", "detectorSecondaryJet_detectorLeadingJet.pdf");
	twoDimension_standard("detectorSecondaryJet_detectorHT", "detectorSecondaryJet_detectorHT.pdf");
	twoDimension_standard("detectorLeadingJet_detectorHT", "detectorLeadingJet_detectorHT.pdf");
	oneDimension_standard("numberOfDetectorJets", "numberOfDetectorJets.pdf");

	oneDimension_twoPlotsSeparate("leadingSquarkPt", "secondarySquarkPt", "Leading Arm", "Secondary Arm", "squarkPt.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingSquarkEta", "secondarySquarkEta", "Leading Arm", "Secondary Arm", "squarkEta.pdf", 0.63, 0.88, 0.76, 0.88, 1);
	oneDimension_threePlotsSeparate("leadingSquarkPt_zeroGluinos", "leadingSquarkPt_oneGluinos", "leadingSquarkPt_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkPt_gluLeading.pdf", 0.63, 0.88, 0.73, 0.88);
	oneDimension_threePlotsSeparate("leadingSquarkEta_zeroGluinos", "leadingSquarkEta_oneGluinos", "leadingSquarkEta_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkEta_gluLeading.pdf", 0.63, 0.88, 0.73, 0.88);
	oneDimension_threePlotsSeparate("secondarySquarkPt_zeroGluinos", "secondarySquarkPt_oneGluinos", "secondarySquarkPt_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkPt_gluSecondary.pdf", 0.63, 0.88, 0.73, 0.88);
	oneDimension_threePlotsSeparate("secondarySquarkEta_zeroGluinos", "secondarySquarkEta_oneGluinos", "secondarySquarkEta_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkEta_gluSecondary.pdf", 0.63, 0.88, 0.73, 0.88);
	twoDimension_standard("leadingSquarkPt_SecondarySquarkPt", "leadingSquarkPt_SecondarySquarkPt.pdf");
	twoDimension_standard("leadingSquarkEta_SecondarySquarkEta", "leadingSquarkEta_SecondarySquarkEta.pdf");
	twoDimension_standard("leadingSquarkPhi_SecondarySquarkPhi", "leadingSquarkPhi_SecondarySquarkPhi.pdf");

	oneDimension_twoPlotsSeparate("leadingQjetPt", "secondaryQjetPt", "Leading Arm", "Secondary Arm", "QjetPt.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingQjetEta", "secondaryQjetEta", "Leading Arm", "Secondary Arm", "QjetEta.pdf", 0.63, 0.88, 0.76, 0.88, 2);

	twoDimension_standard("leadingQjetPt_secondaryQjetPt", "leadingQjetPt_SecondaryQjetPt.pdf");
	twoDimension_standard("leadingQjetEta_secondaryQjetEta", "leadingQjetEta_SecondaryQjetEta.pdf");
	twoDimension_standard("leadingQjetPhi_secondaryQjetPhi", "leadingQjetPhi_SecondaryQjetPhi.pdf");

	oneDimension_twoPlotsSeparate("leadingNlspPt", "secondaryNlspPt", "Leading Arm", "Secondary Arm", "nlspPt.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingNlspEta", "secondaryNlspEta", "Leading Arm", "Secondary Arm", "nlspEta.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	twoDimension_standard("leadingNlspPt_secondaryNlspPt", "leadingNlspPt_secondaryNlspPt.pdf");
	twoDimension_standard("leadingNlspEta_secondaryNlspEta", "leadingNlspEta_secondaryNlspEta.pdf");
	twoDimension_standard("leadingNlspPhi_secondaryNlspPhi", "leadingNlspPhi_secondaryNlspPhi.pdf");

	oneDimension_twoPlotsSeparate("leadingHiggsPt", "secondaryHiggsPt", "Leading Arm", "Secondary Arm", "higgsPt.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingHiggsEta", "secondaryHiggsEta", "Leading Arm", "Secondary Arm", "higgsEta.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingHiggsQjetDphi", "secondaryHiggsQjetDphi", "Leading Arm", "Secondary Arm", "higgsQjetDphi.pdf", 0.63, 0.88, 0.76, 0.88, 2);

	twoDimension_addPlots("leadingQjetPt_leadingHiggsPt", "secondaryQjetPt_secondaryHiggsPt", "QjetPt_HiggsPt_bothArms.pdf");
	twoDimension_addPlots("leadingQjetEta_leadingHiggsEta", "secondaryQjetEta_secondaryHiggsEta", "QjetEta_HiggsEta_bothArms.pdf");
	twoDimension_addPlots("leadingQjetPhi_leadingHiggsPhi", "secondaryQjetPhi_secondaryHiggsPhi", "QjetPhi_HiggsPhi_bothArms.pdf");

	oneDimension_twoPlotsSeparate("leadingLspPt", "secondaryLspPt", "Leading Arm", "Secondary Arm", "lspPt.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingLspEta", "secondaryLspEta", "Leading Arm", "Secondary Arm", "lspEta.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_standard("lspMET", "lspMET.pdf");

	oneDimension_twoPlotsSeparate("leadingBBbarSeperation", "secondaryBBbarSeperation", "Leading Arm", "Secondary Arm", "BBbarSeperation.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingBBbarInvmass", "secondaryBBbarInvmass", "Leading Arm", "Secondary Arm", "BBbarInvmass.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	twoDimension_addPlots_withGradientLines("leadingBBbarSeperation_massHiggsOverPt", "secondaryBBbarSeperation_massHiggsOverPt", "BBbarSeperation_massHiggsOverPt_bothArms.pdf");
	twoDimension_standard("secondaryBBbarSeperation_leadingBBbarSeperation", "secondaryBBbarSeperation_leadingBBbarSeperation.pdf");

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





void PlottingMcSignalStudies::oneDimension_twoPlotsSeparate(std::string histoname1, std::string histoname2, std::string legendName1, std::string legendName2, std::string saveName, double legXmin, double legXmax, double legYmin, double legYmax, int plotWhichHistoFirst)
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
	
	h2->SetLineWidth(2);
	h2->SetLineColor(2);
	// h2->GetXaxis()->SetTitle("");
	h2->GetXaxis()->SetTitleSize(0.06);	
	h2->GetXaxis()->SetLabelSize(0.05);
	// h2->GetYaxis()->SetTitle("");
	h2->GetYaxis()->SetTitleSize(0.06);
	h2->GetYaxis()->SetLabelSize(0.05);

	if (plotWhichHistoFirst==2){
		h2->Draw();
		h1->Draw("same");
	}
	else{
		h1->Draw();
		h2->Draw("same");
	}

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
	double defaultParam = tdrStyle->GetPadRightMargin();
	tdrStyle->SetPadRightMargin(0.10);
    TCanvas* c=new TCanvas("c","c"); 	
	TH2F * h = (TH2F*)f->Get(Form("%s", histoname.c_str()));
	// h->SetLineWidth(2);
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
	tdrStyle->SetPadRightMargin(defaultParam);
}





void PlottingMcSignalStudies::twoDimension_addPlots(std::string histoname1, std::string histoname2, std::string saveName)
{
	double defaultParam = tdrStyle->GetPadRightMargin();
	tdrStyle->SetPadRightMargin(0.10);
    TCanvas* c=new TCanvas("c","c"); 	
	TH2F * h1 = (TH2F*)f->Get(Form("%s", histoname1.c_str()));
	TH2F * h2 = (TH2F*)f->Get(Form("%s", histoname2.c_str()));
	h1->Add(h2);
	// h1->SetLineWidth(2);
	// h1->SetLineColor(2);
	// h1->GetXaxis()->SetTitle("");
	h1->GetXaxis()->SetTitleSize(0.06);	
	h1->GetXaxis()->SetLabelSize(0.05);
	// h1->GetYaxis()->SetTitle("");
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->Draw("colz");
	
	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
	tdrStyle->SetPadRightMargin(defaultParam);
}






void PlottingMcSignalStudies::twoDimension_addPlots_withGradientLines(std::string histoname1, std::string histoname2, std::string saveName)
{
	double defaultParam = tdrStyle->GetPadRightMargin();
	tdrStyle->SetPadRightMargin(0.10);
    TCanvas* c=new TCanvas("c","c"); 	
	TH2F * h1 = (TH2F*)f->Get(Form("%s", histoname1.c_str()));
	TH2F * h2 = (TH2F*)f->Get(Form("%s", histoname2.c_str()));
	h1->Add(h2);
	// h1->SetLineWidth(2);
	// h1->SetLineColor(2);
	// h1->GetXaxis()->SetTitle("");
	h1->GetXaxis()->SetTitleSize(0.06);	
	h1->GetXaxis()->SetLabelSize(0.05);
	// h1->GetYaxis()->SetTitle("");
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->Draw("colz");
	
	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	TF1 *line1 = new TF1("line1","2*x",0,10);
	line1->SetLineWidth(3);
	line1->SetLineStyle(2);
	line1->Draw("same");	

	TF1 *line2 = new TF1("line2","1.8*x",0,10);
	line2->SetLineWidth(3);
	line2->SetLineStyle(2);
	line2->SetLineColor(kMagenta);
	line2->Draw("same");


	leg = new TLegend(0.6,0.3,0.88,0.5);
	leg->AddEntry("line1","y=2x","l");
	leg->AddEntry("line2","y=1.8x","l");
	leg->Draw("same");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
	tdrStyle->SetPadRightMargin(defaultParam);
}




#endif