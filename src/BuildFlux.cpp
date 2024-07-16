#include "definition.h"
#include "binning.h"
#include "utils.h"

TH1D* LoadCounts(unsigned int charge) {
    //Getting counts
    TString fileName= "../IonsSelected/"+getIonPath(charge)+"/"+getIonPath(charge)+"_output.root";
    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
        printf("Errore nell'aprire il file\n");
    }
    TH1D *counts = (TH1D*)file->Get("counts");
    std::cout << fileName << std::endl;
    if (!counts) printf("Errore \n");

    //Getting contamination BelowL1
    TString contName= "../Fragmentation/BelowL1/Fractions/contamination_"+getIonName(charge)+".root";
    TFile *contFile = TFile::Open(contName.Data());
    if (!contFile || contFile->IsZombie()) {
        printf("Errore nell'aprire il file per la purity\n");
    } else {
        auto purity = (TGraph*)contFile->Get("f_vs_r");
        //Removing contamination BelowL1 from counts
        auto hcounts = new TH1D("Counts", "R (GV)", nRbins-1, Rbins);
        for (int ibin = 0; ibin < counts->GetNbinsX(); ibin++) {
        hcounts->SetBinContent(ibin + 1, counts->GetBinContent(ibin + 1) 
                                           *(1 - purity->Eval(counts->GetBinCenter(ibin + 1)))  );
        }
        return hcounts;
    }
    return counts;
}
TH1D *LoadLT(unsigned int charge) {
    //Getting livetime
    TString fileName= "../IonsSelected/"+getIonPath(charge)+"/"+getIonPath(charge)+"_output.root";
    TFile *file = TFile::Open(fileName.Data());
    TH1D *lvt = (TH1D*)file->Get("lvt_25");
    return lvt;
} 
TH1D *LoadAcc(unsigned int charge) {
    //Getting Final Acceptance
    TString fileName= "../IonsSelected/"+getIonPath(charge)+"/"+getIonPath(charge)+"_output.root";
    TFile *file = TFile::Open(fileName.Data());
    TH1D *acc= (TH1D*)file->Get("final_acc_7");
    return acc;
}
TH1D *LoadRBinWidth() {
    auto hwidth = (TH1D*)hist_rig->Clone();
	for(int i=1; i<=hist_rig->GetNbinsX(); ++i) {
		hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
		hwidth->SetBinError(i, 0);
	}
    return hwidth;
}
TH1D *BuildFlux(unsigned int charge) {
    auto counts = LoadCounts(charge);
    auto lvt    = LoadLT(charge);
    auto acc    = LoadAcc(charge);
    auto rWidth = LoadRBinWidth();
    counts->Divide(lvt);
    counts->Divide(acc);
    counts->Divide(rWidth);
    return counts;
}
void MultiplyByXPower(TH1 *hh, double power) {
  for (unsigned int ibin = 0; ibin < hh->GetNbinsX(); ibin++) {
    hh->SetBinContent(ibin + 1, hh->GetBinContent(ibin + 1) * pow(hh->GetBinCenter(ibin + 1), power));
    hh->SetBinError(ibin + 1, hh->GetBinError(ibin + 1) * pow(hh->GetBinCenter(ibin + 1), power));
  }
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: \n");
	    printf("%s <charge> \n", argv[0]);
	    return 1;
	}
    unsigned int charge=atoi(argv[1]);
    TH1D *flux  = BuildFlux(charge);
    TH1D *fluxMultiplied= (TH1D *)flux->Clone();
    MultiplyByXPower(fluxMultiplied,2.70);
    TFile *f = new TFile( Form("../output/Flux/%s.root",getIonName(charge).Data() ) ,"RECREATE");
    f->WriteTObject(flux,"flux");
    f->WriteTObject(fluxMultiplied,"fluxMultiplied");
    f->Close();
    return 1;
}
