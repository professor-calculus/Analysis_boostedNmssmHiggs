#ifndef Analysis_Analysis_boostedNmssmHiggs_PlottingVersionOne_h
#define Analysis_Analysis_boostedNmssmHiggs_PlottingVersionOne_h

// CPP headers
#include <string>
#include <vector>

// ROOT headers
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
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



class PlottingVersionOne
{

public:

	// constructor
	PlottingVersionOne(std::string, std::vector<std::string>, std::vector<double>);

private:

	// objects common for all plots 
	TFile * f;
	std::vector<double> etaBinning;
	std::vector<std::string> doubleBtagWPname;
	std::string outputDirectory;
	std::vector<int> colorPalette;
	std::vector<int> colorPaletteSubtle;

	void hBbMassDist();
	void deltaRMatchDist();
 	void fatJetEtaDist();
 	void fatJetVsHBbPtDist();
 	void effComparingWPs();
	void effComparingEta();

	int SetColor(int position, int maxColors);
	TStyle * TDRStyle();
	TStyle * tdrStyle;
	TLatex * latex = new TLatex();
};

////////////////////////////////
// class function definitions //
////////////////////////////////


//--------constructor---------//

PlottingVersionOne::PlottingVersionOne(std::string inputHistoFile, std::vector<std::string> doubleBtagWPnameD, std::vector<double> etaBinningD)
{
	// open input .root file containing the histograms
 	f = TFile::Open(inputHistoFile.c_str());
 	// load the binning conventions
 	doubleBtagWPname = doubleBtagWPnameD;
 	etaBinning = etaBinningD;
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
 	hBbMassDist();
 	deltaRMatchDist();
 	fatJetEtaDist();
 	fatJetVsHBbPtDist();
 	effComparingWPs();
 	effComparingEta();
}

//-----------public-----------//


//-----------private----------//

void PlottingVersionOne::hBbMassDist()
{
	std::vector<TH1F*> vecHistos;
    TLegend * legend = new TLegend(0.20, 0.60, 0.45, 0.85); //(xmin, ymin, xmax, ymax)

	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
		
	    TCanvas* c=new TCanvas("c","c"); 	
		TH1F * h = (TH1F*)f->Get(Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str()));

		// SETUP HOW YOU WOULD LIKE THE PLOT (tdrStyle does most of this)
		h->SetLineWidth(2);
		// h->SetLineColor(2);
		// h->GetXaxis()->SetTitle("");
		h->GetXaxis()->SetTitleSize(0.06);	
		h->GetXaxis()->SetLabelSize(0.05);
		// h->GetYaxis()->SetTitle("");
		h->GetYaxis()->SetTitleSize(0.06);
		h->GetYaxis()->SetLabelSize(0.05);

		h->Draw();
		// Add stamps
		latex->SetTextAlign(11); // align from left
		latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
		latex->DrawLatex(0.22, 0.72, Form("Tag > %s WP", doubleBtagWPname[iWP].c_str()));
		latex->SetTextAlign(31); // align from right
		latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

		std::string saveName = Form("fatJetMass_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());
		c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
		c->Close();

		// for the combined plot
		vecHistos.push_back(h);
		vecHistos[iWP]->SetLineColor(iWP+1);
		Double_t norm = vecHistos[iWP]->GetEntries();
		vecHistos[iWP]->Scale(1/norm);
		legend->AddEntry(h, Form("%s", doubleBtagWPname[iWP].c_str()), "L");

	} // closes loop through Btag WP labels
	
	// put all of them on the same entry
    TCanvas* c=new TCanvas("c","c");
    for (int iMass=vecHistos.size()-1; iMass>=0; --iMass){
	    vecHistos[iMass]->Draw("same");
	    // if (iMass==0) break;
	}
	legend->Draw();
	// Add stamps
	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
	latex->SetTextAlign(31); // align from right
	latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%sfatJetMass_allDoubleBTagWPNormalised.pdf", outputDirectory.c_str()));
	c->Close();
} // closes the function 'hBbMassDist'







void PlottingVersionOne::deltaRMatchDist()
{
	std::vector<TH1F*> vecHistos;
    TLegend * legend = new TLegend(0.60, 0.60, 0.85, 0.85); //(xmin, ymin, xmax, ymax)

	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
		
	    TCanvas* c=new TCanvas("c","c"); 	
		TH1F * h = (TH1F*)f->Get(Form("deltaR_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str()));

		// SETUP HOW YOU WOULD LIKE THE PLOT (tdrStyle does most of this)
		h->SetLineWidth(2);
		// h->SetLineColor(2);
		// h->GetXaxis()->SetTitle("");
		h->GetXaxis()->SetTitleSize(0.06);	
		h->GetXaxis()->SetLabelSize(0.05);
		// h->GetYaxis()->SetTitle("");
		h->GetYaxis()->SetTitleSize(0.06);
		h->GetYaxis()->SetLabelSize(0.05);

		h->Draw();
		// Add stamps
		latex->SetTextAlign(11); // align from left
		latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
		latex->SetTextAlign(31); // align from right
		latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");
		latex->DrawLatex(0.85, 0.72, Form("Tag > %s WP", doubleBtagWPname[iWP].c_str()));

		std::string saveName = Form("deltaR_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());
		c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
		c->Close();

		// for the combined plot
		vecHistos.push_back(h);
		vecHistos[iWP]->SetLineColor(iWP+1);
		Double_t norm = vecHistos[iWP]->GetEntries();
		vecHistos[iWP]->Scale(1/norm);
		legend->AddEntry(h, Form("%s", doubleBtagWPname[iWP].c_str()), "L");

	} // closes loop through Btag WP labels
	
	// put all of them on the same entry
    TCanvas* c=new TCanvas("c","c");
    for (int iDR=vecHistos.size()-1; iDR>=0; --iDR){
	    vecHistos[iDR]->Draw("same");
	}
	legend->Draw();
	// Add stamps
	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
	latex->SetTextAlign(31); // align from right
	latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%sdeltaR_allDoubleBTagWPNormalised.pdf", outputDirectory.c_str()));
	c->Close();
} // closes the function 'deltaRMatchDist'







void PlottingVersionOne::fatJetEtaDist()
{
	std::vector<TH1F*> vecHistos;
    TLegend * legend = new TLegend(0.20, 0.60, 0.42, 0.85); //(xmin, ymin, xmax, ymax)

	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
		
	    TCanvas* c=new TCanvas("c","c"); 	
		TH1F * h = (TH1F*)f->Get(Form("fatJetEta_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str()));

		// SETUP HOW YOU WOULD LIKE THE PLOT (tdrStyle does most of this)
		h->SetLineWidth(2);
		// h->SetLineColor(2);
		// h->GetXaxis()->SetTitle("");
		h->GetXaxis()->SetTitleSize(0.06);	
		h->GetXaxis()->SetLabelSize(0.05);
		// h->GetYaxis()->SetTitle("");
		h->GetYaxis()->SetTitleSize(0.06);
		h->GetYaxis()->SetLabelSize(0.05);

		h->Draw();
		// Add stamps
		latex->SetTextAlign(11); // align from left
		latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
		latex->DrawLatex(0.17, 0.80, Form("Tag > %s WP", doubleBtagWPname[iWP].c_str()));
		latex->SetTextAlign(31); // align from right
		latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

		std::string saveName = Form("fatJetEta_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());
		c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
		c->Close();

		// for the combined plot
		vecHistos.push_back(h);
		vecHistos[iWP]->SetLineColor(iWP+1);
		Double_t norm = vecHistos[iWP]->GetEntries();
		vecHistos[iWP]->Scale(1/norm);
		legend->AddEntry(h, Form("%s", doubleBtagWPname[iWP].c_str()), "L");

	} // closes loop through Btag WP labels
	
	// put all of them on the same entry
    TCanvas* c=new TCanvas("c","c");
    for (int iDR=vecHistos.size()-1; iDR>=0; --iDR){
	    vecHistos[iDR]->Draw("same");
	}
	legend->Draw();
	// Add stamps
	latex->SetTextAlign(11); // align from left
	latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
	latex->SetTextAlign(31); // align from right
	latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");

	c->SaveAs(Form("%sfatJetEta_allDoubleBTagWPNormalised.pdf", outputDirectory.c_str()));
	c->Close();
} // closes the function 'fatJetEtaDist'







void PlottingVersionOne::fatJetVsHBbPtDist()
{
	double defaultParam = tdrStyle->GetPadRightMargin();
	tdrStyle->SetPadRightMargin(0.10);
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
		
		TCanvas* c=new TCanvas("c","c");
		TH2F * h2 = (TH2F*)f->Get(Form("ptScatter_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str()));

		// SETUP HOW YOU WOULD LIKE THE PLOT (tdrStyle does most of this)
		// h2->GetXaxis()->SetTitle("");
		h2->GetXaxis()->SetTitleSize(0.06);
		h2->GetXaxis()->SetLabelSize(0.05);	
		// h2->GetYaxis()->SetTitle("");
		h2->GetYaxis()->SetTitleSize(0.06);
		h2->GetYaxis()->SetLabelSize(0.05);	
		
		h2->Draw("colz");

		// Add diagnol line		
		c->Update();
		TLine *line = new TLine(0,0,gPad->GetUxmax(),gPad->GetUymax());
		line->SetLineStyle(2);
		line->SetLineWidth(2);
		line->Draw();
		// Add stamps
		latex->SetTextAlign(11);
		latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
		latex->SetTextAlign(31);
		latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");
		latex->DrawLatex(0.85, 0.20, Form("Tag > %s WP", doubleBtagWPname[iWP].c_str()));

		std::string saveName = Form("fatJetVsHBbPtScatter_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());
		c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
		c->Close();

	} // closes loop through Btag WP labels
	tdrStyle->SetPadRightMargin(defaultParam);
} // closes the function 'fatJetVsHBbDist'







void PlottingVersionOne::effComparingWPs()
{
	for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinning.size()-1; ++iEtaBin){
	
	    TCanvas* c=new TCanvas("c","c");
	    TLegend * legend = new TLegend(0.6, 0.15, 0.91, 0.30); //(xmin, ymin, xmax, ymax)
	    legend->SetLineColor(0);
	    legend-> SetNColumns(2);
		TH1F * hDen = (TH1F*)f->Get( Form("effDenominator_eta%.2f-%.2f", etaBinning[iEtaBin], etaBinning[iEtaBin+1]) );

		for (std::vector<std:: string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
	
			TH1F * hNum = (TH1F*)f->Get( Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1]) );
			TEfficiency * hEff = new TEfficiency(*hNum,*hDen);		
			TH1F *hEffHist = (TH1F*) hNum->Clone();
			hEffHist->Divide(hDen);

			// SETUP HOW YOU WOULD LIKE THE PLOT
			int colour = SetColor(iWP, doubleBtagWPname.size());
			hEff->SetLineColor(colour);
			hEff->SetMarkerColor(colour);
			hEffHist->SetMarkerColor(colour);
			hEff->SetLineWidth(2);
			hEffHist->GetXaxis()->SetLabelSize(0.05);
			hEffHist->GetYaxis()->SetLabelSize(0.05);
			// hEffHist->GetXaxis()->SetTitle("");
			hEffHist->GetXaxis()->SetTitleSize(0.06);	
			hEffHist->GetYaxis()->SetTitle("Efficiency");
			hEffHist->GetYaxis()->SetTitleSize(0.06);
			hEffHist->GetYaxis()->SetRangeUser(0,1.11);

			hEffHist->Draw("same, P");
			hEff->Draw("same");
			legend->AddEntry(hEff, Form("%s", doubleBtagWPname[iWP].c_str()), "L");

		} // closes loop through WPs

		legend->Draw();
		// Add Stamps
		latex->SetTextAlign(11);
		latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
		if (iEtaBin==0) latex->DrawLatex(0.18,0.84, Form("|#eta| < %.2f", etaBinning[iEtaBin+1]));
		else latex->DrawLatex(0.18,0.84, Form("%.2f #leq |#eta| < %.2f", etaBinning[iEtaBin], etaBinning[iEtaBin+1]));
		latex->SetTextAlign(31);
		latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");
		// Add line across max efficiency		
		c->Update();
		TLine *line = new TLine(0,1,gPad->GetUxmax(),1);
		line->SetLineStyle(2);
		line->SetLineWidth(2);
		line->Draw();

		std::string saveName = Form("efficiency_eta%.2f-%.2f.pdf", etaBinning[iEtaBin], etaBinning[iEtaBin+1]);
		c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
		c->Close();

	} // closes loop through eta bins 
} // closes the function effComparingWPs








void PlottingVersionOne::effComparingEta()
{
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){

    TCanvas* c=new TCanvas("c","c");
    TLegend * legend = new TLegend(0.30, 0.20, 0.61, 0.36); //(xmin, ymin, xmax, ymax)
    legend->SetLineColor(0);
    // legend-> SetNColumns(2);

		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinning.size()-1; ++iEtaBin){

			TH1F * hDen = (TH1F*)f->Get( Form("effDenominator_eta%.2f-%.2f", etaBinning[iEtaBin], etaBinning[iEtaBin+1]) );
			TH1F * hNum = (TH1F*)f->Get( Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1]) );		
			TEfficiency * hEff = new TEfficiency(*hNum,*hDen);
			TH1F *hEffHist = (TH1F*) hNum->Clone();
			hEffHist->Divide(hDen);

			// SETUP HOW YOU WOULD LIKE THE PLOT			
			int colour = SetColor(iEtaBin, etaBinning.size());
			hEff->SetLineColor(colour);
			hEff->SetMarkerColor(colour);
			hEffHist->SetMarkerColor(colour);
			hEff->SetLineWidth(2);
			hEffHist->GetXaxis()->SetLabelSize(0.05);
			hEffHist->GetYaxis()->SetLabelSize(0.05);
			hEffHist->GetXaxis()->SetTitleSize(0.06);	
			hEffHist->GetYaxis()->SetTitle("Efficiency");
			hEffHist->GetYaxis()->SetTitleSize(0.06);

			hEffHist->GetYaxis()->SetRangeUser(0,1.11);

			hEffHist->Draw("same, P");
			hEff->Draw("same");
			if (iEtaBin==0) legend->AddEntry(hEff, Form("|#eta| < %.2f", etaBinning[iEtaBin+1]), "L");
			else legend->AddEntry(hEff, Form("%.2f #leq |#eta| < %.2f", etaBinning[iEtaBin], etaBinning[iEtaBin+1]), "L");

		} // closes loop through eta bins 

		if (iWP>1){
		legend->SetX1(0.18);
		legend->SetX2(0.49);
		legend->SetY1(0.6);
		legend->SetY2(0.75);
		}
		legend->Draw();

		// Add Stamps
		latex->SetTextAlign(11);
		latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation} Work In Progress");
		latex->DrawLatex(0.18,0.84, Form("Tag > %s WP", doubleBtagWPname[iWP].c_str()));
		latex->SetTextAlign(31);
		latex->DrawLatex(0.92,0.92,"#sqrt{s} = 13 TeV");
		// Add line across max efficiency		
		c->Update();
		TLine *line = new TLine(0,1,gPad->GetUxmax(),1);
		line->SetLineStyle(2);
		line->SetLineWidth(2);
		line->Draw();

		std::string saveName = Form("efficiency_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());
		c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveName.c_str()));
		c->Close();	

	} // closes loop through WPs
} // closes the function effComparingEta





int PlottingVersionOne::SetColor(int posititon, int maxColors)
{
// the plot has 'maxColors' number of colours involved
// this function gives you a colour for the nth histogram

	// the smaller the modifer value the greater the contrast
	// double modifier = 0.05;
	// double modifier = 0.10;
    // double modifier = 0.15;
    // double modifier = 0.20;
	// double modifier = 0.25;	
    double modifier = 0.30;
    double colorIndex;
    int colour(1);
    double fraction = (double)(posititon)/(double)(maxColors-1);

    if( posititon > maxColors-1 || posititon < 0 || maxColors < 0 ) colour = 1;
    else
    {
    	// colorIndex = fraction * gStyle->GetNumberOfColors();
        colorIndex = (fraction * (1.0-2.0*modifier) + modifier) * gStyle->GetNumberOfColors();
        colour = gStyle->GetColorPalette(colorIndex);
    }
    return colour;
}





TStyle * PlottingVersionOne::TDRStyle()
{
	TStyle * tdrStyle = new TStyle("tdrStyle","");

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

	tdrStyle->SetMarkerStyle(20);

	//For the legend
	tdrStyle->SetLegendBorderSize(0);
	tdrStyle->SetLegendFillColor(0);
	tdrStyle->SetLegendFont(42);
	tdrStyle->SetLegendTextSize(0.05);

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
	tdrStyle->SetOptTitle(0);
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
	tdrStyle->SetTitleSize(0.06, "XYZ");
	// tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
	// tdrStyle->SetTitleYSize(Float_t size = 0.02);
	tdrStyle->SetTitleXOffset(0.95);//EDITFROM 0.9
	tdrStyle->SetTitleYOffset(1.25);
	// tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

	// For the axis labels:
	tdrStyle->SetLabelColor(1, "XYZ");
	tdrStyle->SetLabelFont(42, "XYZ");
	tdrStyle->SetLabelOffset(0.007, "XYZ");
	// tdrStyle->SetLabelSize(0.05, "XYZ");
	tdrStyle->SetLabelSize(0.06, "XYZ");

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

#endif