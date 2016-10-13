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

	void hBbMassDist();
	void hBbMassDistDraw(TH1F *, std::string);

 	void fatJetVsHBbDist();
 	void fatJetVsHBbDistDraw(TH2F *, std::string);

 	void effComparingWPs();
 	void effComparingWPsDraw(std::vector<TEfficiency>, std::vector<TH1F*>, std::string);

	void effComparingEta();
	void effComparingEtaDraw(std::vector<TEfficiency> hEffD, std::vector<TH1F*> hDivD, std::string saveNameD);

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
	// make the .pdfs
 	hBbMassDist();
 	fatJetVsHBbDist();
 	effComparingWPs();
 	effComparingEta();

}

//-----------public-----------//


//-----------private----------//

void PlottingVersionOne::hBbMassDist()
{
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
		
		TH1F * h = (TH1F*)f->Get(Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str()));
		std::string saveName = Form("fatJetMass_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());

		// SETUP HOW YOU WOULD LIKE THE PLOT
		h->SetLineWidth(4);
	

		// make plot with complimentary function
		hBbMassDistDraw(h, saveName);

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
		
		TH2F * h2 = (TH2F*)f->Get(Form("ptScatter_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str()));
		std::string saveName = Form("fatJetVsHBbPtScatter_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());

		// SETUP HOW YOU WOULD LIKE THE PLOT
		// h2->SetLineWidth(4);
		

		// make plot with complimentary function
		fatJetVsHBbDistDraw(h2, saveName);

	} // closes loop through Btag WP labels
} // closes the function 'fatJetVsHBbDist'

void PlottingVersionOne::fatJetVsHBbDistDraw(TH2F * h2D, std::string saveNameD)
{
    TCanvas* c=new TCanvas("c","c",650,600);
	h2D->Draw("colz");
	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveNameD.c_str()));
	c->Close();
}






void PlottingVersionOne::effComparingWPs()
{
	for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinning.size()-1; ++iEtaBin){
	
		std::vector<TEfficiency> hEff;
		std::vector<TH1F*> hDiv;

		// get the denominator for this eta bin
		TH1F * hDen = (TH1F*)f->Get( Form("effDenominator_eta%.2f-%.2f", etaBinning[iEtaBin], etaBinning[iEtaBin+1]) );

		for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){
	
			TH1F * hNum = (TH1F*)f->Get( Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1]) );		
			hEff.push_back(TEfficiency(*hNum,*hDen));
			hNum->Divide(hDen);
			hDiv.push_back(hNum);

			// SETUP HOW YOU WOULD LIKE THE PLOT
			hEff[iWP].SetLineColor(iWP+1);
			hDiv[iWP]->SetLineColor(iWP+1);

			// WILL NEED A LEGEND DESPERATELY

		} // closes loop through WPs

		std::string saveName = Form("efficiency_eta%.2f-%.2f.pdf", etaBinning[iEtaBin], etaBinning[iEtaBin+1]);
		// make plot with complimentary function
		effComparingWPsDraw(hEff, hDiv, saveName);

	} // closes loop through eta bins 
} // closes the function effComparingWPs

void PlottingVersionOne::effComparingWPsDraw(std::vector<TEfficiency> hEffD, std::vector<TH1F*> hDivD, std::string saveNameD)
{
    TCanvas* c=new TCanvas("c","c",650,600);
	
    for (size_t iEff=0; iEff<hEffD.size(); ++iEff){
	    hDivD[iEff]->Draw("same, P");
		hEffD[iEff].Draw("same");
	}
	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveNameD.c_str()));
	c->Close();
}






void PlottingVersionOne::effComparingEta()
{
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){

		std::vector<TEfficiency> hEff;
		std::vector<TH1F*> hDiv;

		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinning.size()-1; ++iEtaBin){

			TH1F * hDen = (TH1F*)f->Get( Form("effDenominator_eta%.2f-%.2f", etaBinning[iEtaBin], etaBinning[iEtaBin+1]) );
			TH1F * hNum = (TH1F*)f->Get( Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1]) );		
			hEff.push_back(TEfficiency(*hNum,*hDen));
			hNum->Divide(hDen);
			hDiv.push_back(hNum);

			// SETUP HOW YOU WOULD LIKE THE PLOT			
			hEff[iEtaBin].SetLineColor(iEtaBin); // or maybe do same line style...
			hDiv[iEtaBin]->SetLineColor(iEtaBin);

			// WILL NEED A LEGEND DESPERATELY

		} // closes loop through eta bins 

		std::string saveName = Form("efficiency_%sDoubleBTagWP.pdf", doubleBtagWPname[iWP].c_str());
		// make plot with complimentary function
		effComparingEtaDraw(hEff, hDiv, saveName);

	} // closes loop through WPs
} // closes the function effComparingEta

void PlottingVersionOne::effComparingEtaDraw(std::vector<TEfficiency> hEffD, std::vector<TH1F*> hDivD, std::string saveNameD)
{
    TCanvas* c=new TCanvas("c","c",650,600);
	
    for (size_t iEff=0; iEff<hEffD.size(); ++iEff){
	    hDivD[iEff]->Draw("same, P");
		hEffD[iEff].Draw("same");
	}
	c->SaveAs(Form("%s%s", outputDirectory.c_str(), saveNameD.c_str()));
	c->Close();
}


#endif