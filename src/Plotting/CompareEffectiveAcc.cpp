#include "binning.h"
#include "utils.h"

#include <TLatex.h>

const int nHistPerFile = 1;
std::vector<TString> fluxName = {"final_acc"};

std::vector<TH1D  *> getFluxHist(std::vector<TFile *> v);
void print(std::vector<TH1D*> fluxes, std::vector<unsigned int> chNum, TString timePeriod);

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
    print(fluxes, chNum, timePeriod );
}

std::vector<TH1D *> getFluxHist(std::vector<TFile *> v) {
	std::vector<TH1D*> FluxHist {};
    for (const auto &file : v) 
        for (const auto& name : fluxName) {
			TH1D* hist = (TH1D *)(file->Get(name.Data()) );
			if (hist) {
                FluxHist.push_back(hist);
            }
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
	TString prefix = getPrefixFile(chNum)+"_effAcc";
    for (int j=0; j<nHistPerFile; j++) {
        TCanvas *b = new TCanvas("","",2048,1280);
        if (j==0) b->SaveAs( Form("../output/%s_%s.pdf[", prefix.Data(), timePeriod.Data()) ,"RECREATE");	
		auto a = new TLegend(0.70,0.9,0.9,0.75);
        for (int i=0; i<chNum.size(); i++) {
			b->Update();
            
            TH2D *h1 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.001, 0.077);			
			setTitle(h1,"", "R (GV)", "Effective acc (m^{2}sr)", chNum[i]);			
			formatAxis(h1, chNum.size());		
			formatTitle(h1,chNum.size());			
			adjustZoomY(h1,"",chNum[i]);			
			adjustZoomX(h1,"rate");			
			h1->Draw("SAME");

			TString fluxTitle  = getIonName(chNum[i])+"";
            if (i==0) setTitle(fluxes[j+nHistPerFile*i], "", "R (GV)", "Effective acc (m^{2}sr)", chNum[i]);	
            if (i==1) setTitle(fluxes[j+nHistPerFile*i], "", "R (GV)", "", chNum[i]);	
            if (i==2) setTitle(fluxes[j+nHistPerFile*i], "", "R (GV)", "", chNum[i]);	
			formatAxis(fluxes[j+nHistPerFile*i], chNum.size());
			formatTitle(fluxes[j+nHistPerFile*i],chNum.size());
			getColor(fluxes[j+nHistPerFile*i],chNum[i], "dat");	
			adjustZoomX(fluxes[j+nHistPerFile*i],"Mc acceptance");
			adjustZoomY(fluxes[j+nHistPerFile*i],"Mc acceptance",chNum[i]);
            fluxes[j+nHistPerFile*i]->GetYaxis()->SetLabelFont(62);
            fluxes[j+nHistPerFile*i]->GetXaxis()->SetLabelFont(62);  
            fluxes[j+nHistPerFile*i]->GetYaxis()->SetTitleFont(62);
            fluxes[j+nHistPerFile*i]->GetXaxis()->SetTitleFont(62);
			fluxes[j+nHistPerFile*i]->Draw("SAME");
			a->AddEntry(fluxes[j+nHistPerFile*i],getIonName(chNum[i]));
			gPad->SetGridy();
            gPad->SetLogx();
        }
		a->Draw();
        b->SaveAs( Form("../output/%s_%s.pdf", prefix.Data(), timePeriod.Data()));
        if(j==nHistPerFile-1) b->SaveAs( Form("../output/%s_%s.pdf]", prefix.Data(), timePeriod.Data()));
    }
}