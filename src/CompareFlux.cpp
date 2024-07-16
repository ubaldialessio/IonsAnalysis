#include "binning.h"
#include "utils.h"

#include <TLatex.h>

const int nHistPerFile = 2;
std::vector<TString> fluxName = {"flux","fluxMultiplied"};

std::vector<TFile *> getFiles(int argc, char **argv, std::vector<unsigned int> &chNum);
std::vector<TH1D  *> getFluxHist(std::vector<TFile *> v);
void print(std::vector<TH1D*> fluxes, std::vector<unsigned int> chNum);

int main(int argc, char **argv) {
    if (argc < 2) {
		printf("Usage: \n");
		printf("%s <charge1> <charge2> ... <charge n> \n", argv[0]);
		return 1;
	}
	std::vector<unsigned int> chNum;
	auto files  = getFiles(argc, argv, chNum);
    auto fluxes = getFluxHist(files);
    print(fluxes, chNum );
}

std::vector<TFile *> getFiles(int argc, char **argv, std::vector<unsigned int> &chNum) {
	std::vector<TFile*> files;
	for (int i=1; i<argc; i++) { 
        TString ch = argv[i];
		chNum.push_back(atoi(argv[i]));
		TFile *file = new TFile( Form("../output/Flux/%s.root",getIonName(atoi(ch)).Data()) );
		if (!file || file->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",argv[i]); }
		files.push_back(file);
	}
	return files;
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
void print(std::vector<TH1D*> fluxes, std::vector<unsigned int> chNum) {
    TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	style->cd();
	gStyle->SetOptStat(0);
	TString prefix = getPrefixFile(chNum);
    for (int j=0; j<nHistPerFile; j++) {
        TCanvas *b = new TCanvas("","",2048,1280);
        if (j==0) b->SaveAs( Form("../output/%s.pdf[", prefix.Data()) ,"RECREATE");
        b->Divide(chNum.size(),1e-15,1e-15);		
        for (int i=0; i<chNum.size(); i++) {
			b->cd(i+1);
			b->Update();
			TString fluxTitle  = getIonName(chNum[i])+" flux";
            setTitle(fluxes[j+nHistPerFile*i], fluxTitle, "R (GV)", "", chNum[i]);
            if (j==0) gPad->SetLogy();
			formatAxis(fluxes[j+nHistPerFile*i], chNum.size());
			formatTitle(fluxes[j+nHistPerFile*i],chNum.size());
			getColor(fluxes[j+nHistPerFile*i],chNum[i], "");
			fluxes[j+nHistPerFile*i]->Draw();
			gPad->SetGridy();
            gPad->SetLogx();
        }
        b->SaveAs( Form("../output/%s.pdf", prefix.Data()));
        if(j==nHistPerFile-1) b->SaveAs( Form("../output/%s.pdf]", prefix.Data()));
    }
}