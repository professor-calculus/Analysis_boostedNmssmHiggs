#ifndef Analysis_Analysis_boostedNmssmHiggs_PlottingMcSignalStudiesCMSSW_h
#define Analysis_Analysis_boostedNmssmHiggs_PlottingMcSignalStudiesCMSSW_h

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


class PlottingMcSignalStudiesCMSSW
{

public:

	// constructor
	PlottingMcSignalStudiesCMSSW(std::string,std::string);

private:

	// objects common for all plots 
	TFile * f;
	std::string outputDirectory;
	std::string massTitle;

	// individual plots
	void oneDimension_standard(std::string, std::string);
	void twoDimension_normalisedWithText(std::string, std::string);
	void oneDimension_fourNormalisedPlotsSeparate(std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, double, double, double, double);
	// void oneDimension_twoPlotsAdd();
	void oneDimension_twoPlotsSeparate(std::string, std::string, std::string, std::string, std::string, double, double, double, double, int);
	void oneDimension_threeNormalisedPlotsSeparate(std::string, std::string, std::string, std::string, std::string, std::string, std::string, double, double, double, double);
	void oneDimension_threeNormalisedPlotsSeparateRev(std::string, std::string, std::string, std::string, std::string, std::string, std::string, double, double, double, double);
	void twoDimension_standard(std::string, std::string);
	void twoDimension_forLogic(std::string, std::string);
	void twoDimension_forLogic3x3(std::string, std::string);
	void twoDimension_addPlots(std::string, std::string, std::string);


	// TStyle * TDRStyle(); // now get this function from seperate file
	TStyle * tdrStyle;
	TLatex * latex = new TLatex();
};

////////////////////////////////
// class function definitions //
////////////////////////////////


//--------constructor---------//

PlottingMcSignalStudiesCMSSW::PlottingMcSignalStudiesCMSSW(std::string inputHistoFile, std::string titleName)
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
	oneDimension_standard("qFromGluinoPt_oneGluino", "qFromGluinoPt_oneGluino.pdf");
	
	oneDimension_twoPlotsSeparate("leadingSquarkPt", "secondarySquarkPt", "Leading Arm", "Secondary Arm", "squarkPt.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingSquarkEta", "secondarySquarkEta", "Leading Arm", "Secondary Arm", "squarkEta.pdf", 0.63, 0.88, 0.76, 0.88, 1);

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

	oneDimension_twoPlotsSeparate("drLeadingHiggsQjet", "drSecondaryHiggsQjet", "Leading Arm", "Secondary Arm", "drHiggsQjet.pdf", 0.15, 0.35, 0.76, 0.88, 2);

	oneDimension_twoPlotsSeparate("leadingHiggsPt", "secondaryHiggsPt", "Leading Arm", "Secondary Arm", "higgsPt.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingHiggsEta", "secondaryHiggsEta", "Leading Arm", "Secondary Arm", "higgsEta.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	twoDimension_standard("leadingHiggsPt_secondaryHiggsPt", "leadingHiggsPt_secondaryHiggsPt.pdf");
	twoDimension_standard("leadingHiggsEta_secondaryHiggsEta", "leadingHiggsEta_secondaryHiggsEta.pdf");
	twoDimension_standard("leadingHiggsPhi_secondaryHiggsPhi", "leadingHiggsPhi_secondaryHiggsPhi.pdf");
	oneDimension_twoPlotsSeparate("leadingHiggsQjetDphi", "secondaryHiggsQjetDphi", "Leading Arm", "Secondary Arm", "higgsQjetDphi.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_standard("leadingHiggsSecondaryHiggsDphi", "leadingHiggsSecondaryHiggsDphi.pdf");
	oneDimension_standard("leadingQjetSecondaryQjetDphi", "leadingQjetSecondaryQjetDphi.pdf");

	twoDimension_addPlots("leadingQjetPt_leadingHiggsPt", "secondaryQjetPt_secondaryHiggsPt", "QjetPt_HiggsPt_bothArms.pdf");
	twoDimension_addPlots("leadingQjetEta_leadingHiggsEta", "secondaryQjetEta_secondaryHiggsEta", "QjetEta_HiggsEta_bothArms.pdf");
	twoDimension_addPlots("leadingQjetPhi_leadingHiggsPhi", "secondaryQjetPhi_secondaryHiggsPhi", "QjetPhi_HiggsPhi_bothArms.pdf");

	oneDimension_twoPlotsSeparate("leadingLspPt", "secondaryLspPt", "Leading Arm", "Secondary Arm", "lspPt.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingLspEta", "secondaryLspEta", "Leading Arm", "Secondary Arm", "lspEta.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_standard("lspMET", "lspMET.pdf");

	oneDimension_twoPlotsSeparate("leadingBBbarSeperation", "secondaryBBbarSeperation", "Leading Arm", "Secondary Arm", "BBbarSeperation.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	twoDimension_standard("leadingBBbarSeperation_secondaryBBbarSeperation", "leadingBBbarSeperation_secondaryBBbarSeperation.pdf");
	twoDimension_forLogic("logicBbDr", "logicBbDr.pdf");
	twoDimension_forLogic("logicHiggsPt", "logicHiggsPt.pdf");
	twoDimension_forLogic("logicBbDrHiggsPt", "logicBbDrHiggsPtBothArmsIndividually.pdf");
	twoDimension_forLogic3x3("logicNumberHiggsPtNumberBbDr", "logicNumberHiggsPtNumberBbDr.pdf");

	oneDimension_twoPlotsSeparate("leadingMatchedFatJetPtTwoLooseTagsMatch", "secondaryMatchedFatJetPtTwoLooseTagsMatch", "Leading", "Secondary", "matchedFatJetPtTwoLooseTagsMatch.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingMatchedFatJetPtTwoMediumTagsMatch", "secondaryMatchedFatJetPtTwoMediumTagsMatch", "Leading", "Secondary", "matchedFatJetPtTwoMediumTagsMatch.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	oneDimension_twoPlotsSeparate("leadingMatchedFatJetPtTwoTightTagsMatch", "secondaryMatchedFatJetPtTwoTightTagsMatch", "Leading", "Secondary", "matchedFatJetPtTwoTightTagsMatch.pdf", 0.63, 0.88, 0.76, 0.88, 2);
	twoDimension_standard("leading_secondaryMatchedFatJetPtTwoLooseTagsMatch", "leading_secondaryMatchedFatJetPtTwoLooseTagsMatch.pdf");
	twoDimension_standard("leading_secondaryMatchedFatJetPtTwoMediumTagsMatch", "leading_secondaryMatchedFatJetPtTwoMediumTagsMatch.pdf");
	twoDimension_standard("leading_secondaryMatchedFatJetPtTwoTightTagsMatch", "leading_secondaryMatchedFatJetPtTwoTightTagsMatch.pdf");

 	oneDimension_twoPlotsSeparate("detectorLeadingAk4JetPt", "detectorSecondaryAk4JetPt", "Leading", "Secondary", "detectorAk4JetPt.pdf",  0.63, 0.88, 0.76, 0.88, 2);
 	oneDimension_twoPlotsSeparate("detectorLeadingAk4JetEta", "detectorSecondaryAk4JetEta", "Leading", "Secondary", "detectorAk4JetEta.pdf",  0.63, 0.88, 0.76, 0.88, 1);
	oneDimension_standard("detectorMET", "detectorMET.pdf");
	oneDimension_standard("detectorHT", "detectorHT.pdf");
	oneDimension_standard("detectorMHT", "detectorMHT.pdf");
	oneDimension_standard("detectorMinBiasedDeltaPhi", "detectorMinBiasedDeltaPhi.pdf");
	oneDimension_standard("detectorAlphaT", "detectorAlphaT.pdf");
	oneDimension_standard("detectorMHToverMET", "detectorMHToverMET.pdf");
	twoDimension_standard("detectorMHT_detectorMET", "detectorMHT_detectorMET.pdf"); 
	twoDimension_standard("detectorLeadingAk4JetPt_detectorSecondaryAk4JetPt", "detectorLeadingAk4JetPt_detectorSecondaryAk4JetPt.pdf");
	twoDimension_standard("detectorHT_detectorSecondaryAk4JetPt", "detectorHT_detectorSecondaryAk4JetPt.pdf");
	twoDimension_standard("detectorHT_detectorLeadingAk4JetPt", "detectorHT_detectorLeadingAk4JetPt.pdf");
	twoDimension_forLogic("detectorLogicAlphaTMinBiasedDeltaPhi", "detectorLogicAlphaTMinBiasedDeltaPhi.pdf");
	twoDimension_forLogic("detectorLogicMHToverMETmht", "detectorLogicMHToverMETmht.pdf");

	twoDimension_normalisedWithText("detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52", "detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52.pdf");

	// at end to not mess up the colour scheme
	oneDimension_fourNormalisedPlotsSeparate("detectorNumAk4JetsOver100GeV", "detectorNumAk4JetsOver40GeV", "detectorNumAk4JetsOver20GeV", "detectorNumAk4JetsOver10GeV", "p_{T} > 100 GeV", "p_{T} > 40 GeV", "p_{T} > 20 GeV", "p_{T} > 10 GeV", "detectorNumAk4Jets.pdf", 0.63, 0.88, 0.68, 0.88);

	oneDimension_threeNormalisedPlotsSeparate("fatJetNumberLooseDoubleBTagsNoMatching", "fatJetNumberMediumDoubleBTagsNoMatching", "fatJetNumberTightDoubleBTagsNoMatching", "loose", "medium", "tight", "fatJetNumberDoubleBTagsNoMatching.pdf", 0.66, 0.88, 0.73, 0.88);
	oneDimension_threeNormalisedPlotsSeparate("fatJetNumberLooseDoubleBTagsWithMatching", "fatJetNumberMediumDoubleBTagsWithMatching", "fatJetNumberTightDoubleBTagsWithMatching", "loose", "medium", "tight", "fatJetNumberDoubleBTagsWithMatching.pdf", 0.66, 0.88, 0.73, 0.88);
	oneDimension_threeNormalisedPlotsSeparateRev("leadingSquarkPt_zeroGluinos", "leadingSquarkPt_oneGluinos", "leadingSquarkPt_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkPt_gluLeading.pdf", 0.65, 0.88, 0.73, 0.88);
	oneDimension_threeNormalisedPlotsSeparateRev("leadingSquarkEta_zeroGluinos", "leadingSquarkEta_oneGluinos", "leadingSquarkEta_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkEta_gluLeading.pdf", 0.65, 0.88, 0.73, 0.88);
	oneDimension_threeNormalisedPlotsSeparateRev("secondarySquarkPt_zeroGluinos", "secondarySquarkPt_oneGluinos", "secondarySquarkPt_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkPt_gluSecondary.pdf", 0.65, 0.88, 0.73, 0.88);
	oneDimension_threeNormalisedPlotsSeparateRev("secondarySquarkEta_zeroGluinos", "secondarySquarkEta_oneGluinos", "secondarySquarkEta_twoGluinos", "Zero Gluinos", "One Gluino", "Two Gluinos", "squarkEta_gluSecondary.pdf", 0.65, 0.88, 0.73, 0.88);

}

//-----------public-----------//


//-----------private----------//

void PlottingMcSignalStudiesCMSSW::oneDimension_standard(std::string histoname, std::string saveName)
{
    TCanvas* c=new TCanvas("c","c"); 	

	TH1F * h = (TH1F*)f->Get(Form("%s", histoname.c_str()));
	h->SetLineWidth(2);
	h->SetLineColor(kBlue+1);
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




void PlottingMcSignalStudiesCMSSW::twoDimension_normalisedWithText(std::string histoname, std::string saveName)
{
	double defaultParam = tdrStyle->GetPadRightMargin();
	tdrStyle->SetPadRightMargin(0.10);

    TCanvas* c=new TCanvas("c","c"); 	

	TH2F * h = (TH2F*)f->Get(Form("%s", histoname.c_str()));
	Double_t norm = h->GetEntries();
	h->Scale(1/norm);
	h->SetLineWidth(2);
	h->SetLineColor(kBlue+1);
	h->SetMarkerSize(4);
	// h->GetXaxis()->SetTitle("");
	h->GetXaxis()->SetTitleSize(0.06);	
	h->GetXaxis()->SetLabelSize(0.05);
	// h->GetYaxis()->SetTitle("");
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.05);
	h->Draw("colz, TEXT");

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
	tdrStyle->SetPadRightMargin(defaultParam);
}





void PlottingMcSignalStudiesCMSSW::oneDimension_fourNormalisedPlotsSeparate(std::string histoname1, std::string histoname2, std::string histoname3, std::string histoname4, std::string legendName1, std::string legendName2, std::string legendName3, std::string legendName4, std::string saveName, double legXmin, double legXmax, double legYmin, double legYmax)
{
    TCanvas* c=new TCanvas("c","c"); 	
	TH1F * h1 = (TH1F*)f->Get(Form("%s", histoname1.c_str()));
	TH1F * h2 = (TH1F*)f->Get(Form("%s", histoname2.c_str()));
	TH1F * h3 = (TH1F*)f->Get(Form("%s", histoname3.c_str()));
	TH1F * h4 = (TH1F*)f->Get(Form("%s", histoname4.c_str()));	

	Double_t norm = h1->GetEntries();
	h1->Scale(1/norm);
	norm = h2->GetEntries();
	h2->Scale(1/norm);
	norm = h3->GetEntries();
	h3->Scale(1/norm);
	norm = h4->GetEntries();
	h4->Scale(1/norm);

	int colour1 = SetColor(0, 4);
	int colour2 = SetColor(1, 4);
	int colour3 = SetColor(2, 4);
	int colour4 = SetColor(3, 4);

	h1->SetLineWidth(2);
	h1->SetLineColor(colour1);
	// h1->GetXaxis()->SetTitle("");
	h1->GetXaxis()->SetTitleSize(0.06);	
	h1->GetXaxis()->SetLabelSize(0.05);
	// h1->GetYaxis()->SetTitle("");
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	
	h2->SetLineWidth(2);
	h2->SetLineColor(colour2);
	// h2->GetXaxis()->SetTitle("");
	h2->GetXaxis()->SetTitleSize(0.06);	
	h2->GetXaxis()->SetLabelSize(0.05);
	// h2->GetYaxis()->SetTitle("");
	h2->GetYaxis()->SetTitleSize(0.06);
	h2->GetYaxis()->SetLabelSize(0.05);

	h3->SetLineWidth(2);
	h3->SetLineColor(colour3);
	// h3->GetXaxis()->SetTitle("");
	h3->GetXaxis()->SetTitleSize(0.06);	
	h3->GetXaxis()->SetLabelSize(0.05);
	// h3->GetYaxis()->SetTitle("");
	h3->GetYaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetLabelSize(0.05);

	h4->SetLineWidth(2);
	h4->SetLineColor(colour4);
	// h4->GetXaxis()->SetTitle("");
	h4->GetXaxis()->SetTitleSize(0.06);	
	h4->GetXaxis()->SetLabelSize(0.05);
	// h4->GetYaxis()->SetTitle("");
	h4->GetYaxis()->SetTitleSize(0.06);
	h4->GetYaxis()->SetLabelSize(0.05);

	// if (plotWhichHistoFirst==4){
	// 	h4->Draw();
	// 	h3->Draw("same");
	// 	h2->Draw("same");
	// 	h1->Draw("same");
	// }
	// else{
		h1->Draw();
		h2->Draw("same");
		h3->Draw("same");
		h4->Draw("same");
	// }

	TLegend * legend = new TLegend();
	legend->SetX1NDC(legXmin);
	legend->SetX2NDC(legXmax);
	legend->SetY1NDC(legYmin);
	legend->SetY2NDC(legYmax);
	legend->AddEntry(h1, legendName1.c_str(), "L");
	legend->AddEntry(h2, legendName2.c_str(), "L");
	legend->AddEntry(h3, legendName3.c_str(), "L");
	legend->AddEntry(h4, legendName4.c_str(), "L");
	legend->Draw("same");

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();


}



void PlottingMcSignalStudiesCMSSW::oneDimension_twoPlotsSeparate(std::string histoname1, std::string histoname2, std::string legendName1, std::string legendName2, std::string saveName, double legXmin, double legXmax, double legYmin, double legYmax, int plotWhichHistoFirst)
{
    TCanvas* c=new TCanvas("c","c"); 	
	TH1F * h1 = (TH1F*)f->Get(Form("%s", histoname1.c_str()));
	TH1F * h2 = (TH1F*)f->Get(Form("%s", histoname2.c_str()));	
	
	h1->SetLineWidth(2);
	h1->SetLineColor(kBlue+1);
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






void PlottingMcSignalStudiesCMSSW::oneDimension_threeNormalisedPlotsSeparate(std::string histoname1, std::string histoname2, std::string histoname3, std::string legendName1, std::string legendName2, std::string legendName3, std::string saveName, double legXmin, double legXmax, double legYmin, double legYmax)
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

	int colour1 = SetColor(0, 3);
	int colour2 = SetColor(1, 3);
	int colour3 = SetColor(2, 3);

	h1->SetLineWidth(3);
	h1->SetLineStyle(1);
	h1->SetLineColor(colour1);
	// h1->GetXaxis()->SetTitle("");
	h1->GetXaxis()->SetTitleSize(0.06);	
	h1->GetXaxis()->SetLabelSize(0.05);
	// h1->GetYaxis()->SetTitle("");
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->Draw();

	h2->SetLineWidth(4);
	h2->SetLineColor(colour2);
	h2->SetLineStyle(7);
	// h2->GetXaxis()->SetTitle("");
	h2->GetXaxis()->SetTitleSize(0.06);	
	h2->GetXaxis()->SetLabelSize(0.05);
	// h2->GetYaxis()->SetTitle("");
	h2->GetYaxis()->SetTitleSize(0.06);
	h2->GetYaxis()->SetLabelSize(0.05);
	h2->Draw("same");

	h3->SetLineWidth(5);
	h3->SetLineColor(colour3);
	h3->SetLineStyle(2);
	// h3->GetXaxis()->SetTitle("");
	h3->GetXaxis()->SetTitleSize(0.06);	
	h3->GetXaxis()->SetLabelSize(0.05);
	// h3->GetYaxis()->SetTitle("");
	h3->GetYaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetLabelSize(0.05);
	h3->Draw("same");
	
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









void PlottingMcSignalStudiesCMSSW::oneDimension_threeNormalisedPlotsSeparateRev(std::string histoname1, std::string histoname2, std::string histoname3, std::string legendName1, std::string legendName2, std::string legendName3, std::string saveName, double legXmin, double legXmax, double legYmin, double legYmax)
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

	int colour1 = SetColor(0, 3);
	int colour2 = SetColor(1, 3);
	int colour3 = SetColor(2, 3);

	h3->SetLineWidth(2);
	h3->SetLineColor(colour3);
	// h3->GetXaxis()->SetTitle("");
	h3->GetXaxis()->SetTitleSize(0.06);	
	h3->GetXaxis()->SetLabelSize(0.05);
	// h3->GetYaxis()->SetTitle("");
	h3->GetYaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetLabelSize(0.05);
	h3->Draw();
	
	h2->SetLineWidth(2);
	h2->SetLineColor(colour2);
	// h2->GetXaxis()->SetTitle("");
	h2->GetXaxis()->SetTitleSize(0.06);	
	h2->GetXaxis()->SetLabelSize(0.05);
	// h2->GetYaxis()->SetTitle("");
	h2->GetYaxis()->SetTitleSize(0.06);
	h2->GetYaxis()->SetLabelSize(0.05);
	h2->Draw("same");

	h1->SetLineWidth(2);
	h1->SetLineColor(colour1);
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





void PlottingMcSignalStudiesCMSSW::twoDimension_standard(std::string histoname, std::string saveName)
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





void PlottingMcSignalStudiesCMSSW::twoDimension_forLogic(std::string histoname, std::string saveName)
{
	double defaultParam = tdrStyle->GetPadRightMargin();
	tdrStyle->SetPadRightMargin(0.10);
    TCanvas* c=new TCanvas("c","c"); 	
	TH2F * h = (TH2F*)f->Get(Form("%s", histoname.c_str()));
	Double_t norm = h->GetEntries();
	h->Scale(1/norm);
	// h->SetLineWidth(2);
	// h->SetLineColor(2);
	h->SetMarkerSize(4);
	// h->GetXaxis()->SetTitle("");
	h->GetXaxis()->SetTitleSize(0.06);	
	h->GetXaxis()->SetLabelSize(0.05);
	// h->GetYaxis()->SetTitle("");
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.05);
	h->Draw("colz, TEXT");
	
	c->Update();
	TLine *line = new TLine(0,1,2,1);
	line->SetLineStyle(2);
	line->SetLineWidth(3);
	line->Draw();
	TLine *line2 = new TLine(1,0,1,2);
	line2->SetLineStyle(2);
	line2->SetLineWidth(3);
	line2->Draw();

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
	tdrStyle->SetPadRightMargin(defaultParam);
}





void PlottingMcSignalStudiesCMSSW::twoDimension_forLogic3x3(std::string histoname, std::string saveName)
{
	double defaultParam = tdrStyle->GetPadRightMargin();
	tdrStyle->SetPadRightMargin(0.10);
    TCanvas* c=new TCanvas("c","c"); 	
	TH2F * h = (TH2F*)f->Get(Form("%s", histoname.c_str()));
	Double_t norm = h->GetEntries();
	h->Scale(1/norm);
	// h->SetLineWidth(2);
	// h->SetLineColor(2);
	h->SetMarkerSize(3);
	// h->GetXaxis()->SetTitle("");
	h->GetXaxis()->SetTitleSize(0.06);	
	h->GetXaxis()->SetLabelSize(0.05);
	// h->GetYaxis()->SetTitle("");
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.05);
	h->Draw("colz, TEXT");
	
	c->Update();
	TLine *line = new TLine(0,1,3,1);
	line->SetLineStyle(2);
	line->SetLineWidth(3);
	line->Draw();
	TLine *line2 = new TLine(0,2,3,2);
	line2->SetLineStyle(2);
	line2->SetLineWidth(3);
	line2->Draw();
	TLine *line3 = new TLine(1,0,1,3);
	line3->SetLineStyle(2);
	line3->SetLineWidth(3);
	line3->Draw();
	TLine *line4 = new TLine(2,0,2,3);
	line4->SetLineStyle(2);
	line4->SetLineWidth(3);
	line4->Draw();

	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,massTitle.c_str());
	// latex->SetTextAlign(31); // align from right
	// latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
	c->Close();
	tdrStyle->SetPadRightMargin(defaultParam);
}







void PlottingMcSignalStudiesCMSSW::twoDimension_addPlots(std::string histoname1, std::string histoname2, std::string saveName)
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



#endif