#ifndef Analysis_Analysis_boostedNmssmHiggs_PlottingVersionOne_h
#define Analysis_Analysis_boostedNmssmHiggs_PlottingVersionOne_h

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

// CMSSW headers

// Headers from this package




class PlottingVersionOne
{

public:

	// constructor
	PlottingVersionOne(std::string, std::vector<std::string>, std::vector<double>);

private:
	
	TFile * f;
	std::vector<double> etaBinning;
	std::vector<std::string> doubleBtagWPname;
	std::string outputDirectory;

	void hBbMassDist();
	void hBbMassDistDraw(TH1F *, std::string);

 	void fatJetVsHBbDist();
 	void fatJetVsHBbDistDraw(TH2F *, std::string);


 	// void effComparingWPs();
 	// void effComparingEta();
};



////////////////////////////////
// class function definitions //
////////////////////////////////

//--------constructor---------//

PlottingVersionOne::PlottingVersionOne(std::string inputHistoFile, std::vector<std::string> doubleBtagWPnameD, std::vector<double> etaBinningD)
{
 	f = TFile::Open(inputHistoFile.c_str());
 	doubleBtagWPname = doubleBtagWPnameD;
 	etaBinning = etaBinningD;

 	// std::cout << "debug1" << std::endl;

	for (size_t c = inputHistoFile.size()-1; c > 0; --c){
		std::string forwardSlash = "/";
		if (inputHistoFile[c] == forwardSlash[0]){
			outputDirectory = inputHistoFile.substr(0, c+1);
			break;
		}
	}
 	// std::cout << "debug2" << std::endl;


 	hBbMassDist();
 	fatJetVsHBbDist();
 	// effComparingWPs();
 	// effComparingEta();

}

//-----------public-----------//



//-----------private----------//

void PlottingVersionOne::hBbMassDist()
{
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
		
		TH1F * h = (TH1F*)f->Get(Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str()));
		std::string saveName = Form("fatJetMass_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());

		// etc. do lot's of things here
		h->SetLineWidth(4);
		
		// make plot
		hBbMassDistDraw(h, saveName);
		// delete h;
	} // closes loop through Btag WP labels
} // closes the function 'hBbMassDist'


void PlottingVersionOne::hBbMassDistDraw(TH1F * hD, std::string saveNameD)
{
    TCanvas* c=new TCanvas("c","c",650,600);
	hD->Draw();
	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveNameD.c_str()));
	c->Close();
}




void PlottingVersionOne::fatJetVsHBbDist()
{
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
		
		TH2F * h = (TH2F*)f->Get(Form("ptScatter_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str()));
		std::string saveName = Form("fatJetVsHBbPtScatter_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());

		// etc. do lot's of things here
		// h->SetLineWidth(4);
		
		// make plot
		fatJetVsHBbDistDraw(h, saveName);
		// delete h;
	} // closes loop through Btag WP labels
} // closes the function 'hBbMassDist'


void PlottingVersionOne::fatJetVsHBbDistDraw(TH2F * h2D, std::string saveNameD)
{
    TCanvas* c=new TCanvas("c","c",650,600);
	h2D->Draw("colz");
	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveNameD.c_str()));
	c->Close();
}








#endif