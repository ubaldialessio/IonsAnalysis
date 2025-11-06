#include "binning.h"
#include "utils.h"

#include <TLatex.h>

const int nHistPerFile = 1;
std::vector<TString> fluxName = {"rateMulti"};

std::vector<TH1D  *> getFluxHist(std::vector<TFile *> v);
void print(std::vector<TH1D*> fluxes, std::vector<unsigned int> chNum, TString timePeriod);
TString NumberWithCommas(int value);

int main(int argc, char **argv) {
    if (argc < 3) {
		printf("Usage: \n");
		printf("%s <time> <charge1> <charge2> ... <charge n> \n", argv[0]);
		return 1;
	}
	TString timePeriod = argv[1];
	std::vector<unsigned int> chNum;
	auto files  = getFilesFromUnfoldingFactor(argc, argv, chNum, timePeriod,"");
    auto fluxes = getFluxHist(files);
    print(fluxes, chNum, timePeriod);
}
std::vector<TH1D *> getFluxHist(std::vector<TFile *> v) {
	std::vector<TH1D*> FluxHist {};
    for (const auto &file : v) 
        for (const auto& name : fluxName) {
			TH1D* hist = (TH1D *)(file->Get(name.Data()) );
			if (hist) FluxHist.push_back(hist);
		}
	return FluxHist;
}
void print(std::vector<TH1D*> fluxes, std::vector<unsigned int> chNum, TString timePeriod) {
    TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
    setLogon();
	style->cd();
	gStyle->SetOptStat(0);
	TString prefix = getPrefixFile(chNum)+Form("_rates_%s",timePeriod.Data());
    for (int j=0; j<nHistPerFile; j++) {
        TCanvas *b = new TCanvas("","",2048,1280);
        if (j==0) b->SaveAs( Form("../output/%s.pdf[", prefix.Data()) ,"RECREATE");	
		auto a = new TLegend(0.10,0.9,0.9,0.65);
		a->SetNColumns(2);
		a->SetBorderSize(0);
		a->SetFillStyle(0);   
		a->SetTextFont(62);
        for (int i=0; i<chNum.size(); i++) {
			b->Update();
            
            /*TH2D *h1 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.001, 1.15);			
			setTitle(h1,"", "R (GV)", "Rate R^{2.7} (s^{-1} GV^{1.7})", chNum[i]);			
			formatAxis(h1, chNum.size());		
			formatTitle(h1,chNum.size());			
			adjustZoomY(h1,"",chNum[i]);			
			adjustZoomX(h1,"");			
			h1->Draw("SAME");*/	
			gPad->SetLogy();
			TString fluxTitle  = getIonName(chNum[i])+" flux";
            if (i==0) setTitle(fluxes[j+nHistPerFile*i], "", "R (GV)", "Rate R^{2.7} (s^{-1} GV^{1.7})", chNum[i]);	
            else setTitle(fluxes[j+nHistPerFile*i], "", "R (GV)", "", chNum[i]);
			fluxes[j+nHistPerFile*i]->SetTitle("");
			formatAxis(fluxes[j+nHistPerFile*i], 1);
			formatTitle(fluxes[j+nHistPerFile*i],1);
			getColor(fluxes[j+nHistPerFile*i],chNum[i], "dat");	
			adjustZoomX(fluxes[j+nHistPerFile*i],"rate");
			adjustZoomY(fluxes[j+nHistPerFile*i],"rate",chNum[i]);
			formatMarkerSize(fluxes[j+nHistPerFile*i],1);
            fluxes[j+nHistPerFile*i]->GetYaxis()->SetLabelFont(62);
            fluxes[j+nHistPerFile*i]->GetXaxis()->SetLabelFont(62);  
            fluxes[j+nHistPerFile*i]->GetYaxis()->SetTitleFont(62);
            fluxes[j+nHistPerFile*i]->GetXaxis()->SetTitleFont(62);
			fluxes[j+nHistPerFile*i]->Draw("SAME");
			a->AddEntry(fluxes[j+nHistPerFile*i],Form("%s, %s raw counts",getIonName(chNum[i]).Data(),NumberWithCommas(getNucleiCount(chNum[i],timePeriod)).Data()));
			a->SetTextSize(0.031);
			gPad->SetGridy();
            gPad->SetLogx();
        }
		a->Draw();
        b->SaveAs( Form("../output/%s.pdf", prefix.Data()));
        if(j==nHistPerFile-1) b->SaveAs( Form("../output/%s.pdf]", prefix.Data()));
    }
}

TString NumberWithCommas(int value) {
    TString num = TString::Format("%d", value);
    int len = num.Length();
    for (int i = len - 3; i > 0; i -= 3)
        num.Insert(i, ".");
    return num;
}