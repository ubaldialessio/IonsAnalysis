#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"
#include "SplineUtility.h"

TH1D* LoadCounts(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting counts
    TString fileName= Form("../IonsSelected/"+getIonPath(charge)+"/Counts/%scounts.root",inp_sec_track.Data());
    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
        printf("Errore nell'aprire il file\n");
    }
    int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
    auto counts = (TH1D *)((TH2D *)file->Get("IL1/rigidity"))->ProjectionY("r", start_bin, stop_bin);
    if (!counts) printf("Errore \n");

    //Getting contamination BelowL1
    TString contName= "../Fragmentation/BelowL1/Fractions/contamination_"+getIonName(charge)+".root";
    TFile *contFile = TFile::Open(contName.Data());
    //Getting contamination AboveL1
    TString above = Form("../IonsSelected/"+getIonPath(charge)+"/Fragments/Result/%sfragment_acc7.root",inp_sec_track.Data());
    std::cout << "Using this TOI: " << above << std::endl;
    TFile *aboveFile = TFile::Open(above.Data());
    if (!contFile || contFile->IsZombie()) {
        printf("Errore nell'aprire il file per la purity\n");
    } if (!aboveFile || aboveFile->IsZombie()) {
        printf("Errore nell'aprire aboveL1\n");
    } else {
        auto purity = (TF1*)contFile->Get("purityFit");
        auto fragm  = (TH1D*)aboveFile->Get("fraction_spline_tot");
        //Removing contamination BelowL1 and AboveL1 from counts
        auto hcounts = new TH1D("Counts", "R (GV)", nRbins_HighZ-1, Rbins_HighZ);
        for (int ibin = 0; ibin < hcounts->GetNbinsX(); ibin++) {
            double N = counts->GetBinContent(ibin + 1);
            double P = purity->Eval(counts->GetBinCenter(ibin+1));
            double F = fragm->GetBinContent(ibin + 1);
            hcounts->SetBinContent(ibin + 1, N *(P) - N*F );
            double sigma_N = counts->GetBinError(ibin + 1);
            //double sigma_P = purity->GetBinError(ibin + 1);
            double sigma_F = fragm->GetBinError(ibin + 1);
            double sigma_hcounts = sqrt(pow((1 - P) * (1 - F) * sigma_N, 2) +
                                        //pow(N * (1 - F) * sigma_P, 2) +
                                        pow(N * (1 - P) * sigma_F, 2));
            hcounts->SetBinError(ibin + 1, sigma_hcounts);
        }
        return hcounts;
    }
    return counts;
}
TH1D *LoadCountsWithBackground(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting counts
    TString fileName= Form("../IonsSelected/"+getIonPath(charge)+"/Counts/%scounts.root",inp_sec_track.Data());
    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
        printf("Errore nell'aprire il file\n");
    }
    int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
    auto counts = (TH1D *)((TH2D *)file->Get("IL1/rigidity"))->ProjectionY("r", start_bin, stop_bin);
    if (!counts) printf("Errore \n");
    return counts;   
}
TH1D *LoadLT(unsigned int charge, TString timePeriod) {
    //Getting livetime
    TString fileName= "../IonsSelected/"+getIonPath(charge)+"/Livetime/livetime.root";
    TFile *file = TFile::Open(fileName.Data());
    int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
    TH1D *lvt = (TH1D *)((TH2D *)file->Get("lvt_25"))->ProjectionY("q", start_bin, stop_bin);
    return lvt;
} 
TH1D *LoadUnfoldedAcceptance(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting Final Acceptance
    TString fileName = "../IonsSelected/"+getIonPath(charge)+"/Unfolding/"+getIonPath(charge)+
                       Form("_UnfoldingFactor%s_%s.root",inp_sec_track.Data(),timePeriod.Data());
    TFile *file = TFile::Open(fileName.Data());
    TH1D *acc = (TH1D*)file->Get("final_acc");
    return acc;
}
TH1D *LoadRawAcceptance(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting Final Acceptance
    TString fileName= Form("../IonsSelected/"+getIonPath(charge)+"/RawAcceptance/rawacc.root");
    TFile *file = TFile::Open(fileName.Data());
    auto mc_pass_gen = (TH1D*)file->Get("mc_pass_gen");
    auto mc_samp = (TH1D*)file->Get("mc_samp");
    mc_pass_gen->Rebin(9);
    mc_samp->Rebin(9);
    auto acc = divide(mc_pass_gen, mc_samp, "efficiency");
    acc->Scale(TMath::Pi() * 3.9 * 3.9);
    auto acc_range= SetFitLimits(SplineUtility::Efficiency::Acc,charge);
    auto acc_knots =SetKnots(SplineUtility::Efficiency::Acc,charge);
    auto spline_acc = autospline(acc, acc_range.first, acc_range.second,acc_knots.first,acc_knots.second);
    return spline_acc;
}
TH1D *LoadRBinWidth() {
    auto hwidth = (TH1D*)hist_rig_highZ->Clone();
	for(int i=0; i<hist_rig_highZ->GetNbinsX(); ++i) {
		hwidth->SetBinContent(i+1, hwidth->GetBinWidth(i+1));
		hwidth->SetBinError(i+1, 0);
	}
    return hwidth;
}
TH1D *LoadUnfoldingFactor(unsigned int charge, TString timePeriod) {
        TString fileName = "../IonsSelected/"+getIonPath(charge)+"/Unfolding/"+getIonPath(charge)+
                       Form("_UnfoldingFactor_%s.root",timePeriod.Data());
    TFile *file = TFile::Open(fileName.Data());
    TH1D *u = (TH1D*)file->Get("unf_factor");
    return autospline(u,2.15,1000,5,20);
}
TH1D *BuildFlux(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    auto counts = LoadCounts(charge, timePeriod,inp_sec_track);
    auto lvt    = LoadLT(charge, timePeriod);    
    auto acc    = LoadUnfoldedAcceptance(charge, timePeriod,inp_sec_track);
    auto rWidth = LoadRBinWidth();
    for (int ibin = 1; ibin <= counts->GetNbinsX(); ibin++) {
        double N = counts->GetBinContent(ibin);
        double LT = lvt->GetBinContent(ibin);
        double A = acc->GetBinContent(ibin);
        double dR = acc->GetBinWidth(ibin);
        double sigma_N = counts->GetBinError(ibin);
        double sigma_LT = lvt->GetBinError(ibin);
        double sigma_A = acc->GetBinError(ibin);
        double sigma_dR = rWidth->GetBinError(ibin);
        double flux = (N / (LT * A * dR));
        double sigma_flux = flux * sqrt(
            pow(sigma_N / N, 2) +
            pow(sigma_LT / LT, 2) +
            pow(sigma_A / A, 2) +
            pow(sigma_dR / dR, 2)
        );
        counts->SetBinContent(ibin, flux);
        counts->SetBinError(ibin, sigma_flux);
    }
    return counts;
}
TH1D *BuildFluxWithBackground(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    auto counts = LoadCountsWithBackground(charge, timePeriod,inp_sec_track);
    auto lvt    = LoadLT(charge, timePeriod);    
    auto acc    = LoadUnfoldedAcceptance(charge, timePeriod,inp_sec_track);
    auto rWidth = LoadRBinWidth();
    for (int ibin = 1; ibin <= counts->GetNbinsX(); ibin++) {
        double N = counts->GetBinContent(ibin);
        double LT = lvt->GetBinContent(ibin);
        double A = acc->GetBinContent(ibin);
        double dR = rWidth->GetBinContent(ibin);
        double sigma_N = counts->GetBinError(ibin);
        double sigma_LT = lvt->GetBinError(ibin);
        double sigma_A = acc->GetBinError(ibin);
        double sigma_dR = rWidth->GetBinError(ibin);
        double flux = (N / (LT * A * dR));
        double sigma_flux = flux * sqrt(
            pow(sigma_N / N, 2) +
            pow(sigma_LT / LT, 2) +
            pow(sigma_A / A, 2) +
            pow(sigma_dR / dR, 2)
        );

        counts->SetBinContent(ibin, flux);
        counts->SetBinError(ibin, sigma_flux);
    }
    return counts;
}
TH1D *BuildRawFlux(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    auto counts = LoadCountsWithBackground(charge, timePeriod,inp_sec_track);
    auto lvt    = LoadLT(charge, timePeriod);    
    auto acc    = LoadRawAcceptance(charge, timePeriod,inp_sec_track);
    auto rWidth = LoadRBinWidth();
    for (int ibin = 1; ibin <= counts->GetNbinsX(); ibin++) {
        double N = counts->GetBinContent(ibin);
        double LT = lvt->GetBinContent(ibin);
        double A = acc->GetBinContent(ibin);
        double dR = rWidth->GetBinContent(ibin);
        double sigma_N = counts->GetBinError(ibin);
        double sigma_LT = lvt->GetBinError(ibin);
        double sigma_A = acc->GetBinError(ibin);
        double sigma_dR = rWidth->GetBinError(ibin);
        double flux = (N / (LT * A * dR));
        double sigma_flux = flux * sqrt(
            pow(sigma_N / N, 2) +
            pow(sigma_LT / LT, 2) +
            pow(sigma_A / A, 2) +
            pow(sigma_dR / dR, 2)
        );

        counts->SetBinContent(ibin, flux);
        counts->SetBinError(ibin, sigma_flux);
    }
    return counts;
}
TH1D *LoadTotalCorrection(unsigned int charge, TString timePeriod) {
    TString fileName = Form("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/%s_TotalCorrections%s_%s.root",getIonPath(charge).Data(),"",timePeriod.Data());
    TFile *f = TFile::Open(fileName.Data());
    auto totCorr = (TH1D*)f->Get("final_damc_tot");
    if (totCorr) std::cout << "retrieving tot corr\n";
    return totCorr;
}
TF1 *LoadPurity(unsigned int charge) {
    TString contName= "../Fragmentation/BelowL1/Fractions/contamination_"+getIonName(charge)+".root";
    TFile *contFile = TFile::Open(contName.Data());
    auto purity = (TF1*)contFile->Get("purityFit");
    return purity;
}
TH1D *LoadToi(unsigned int charge,TString inp_sec_track) {
    TString above = Form("../IonsSelected/"+getIonPath(charge)+"/Fragments/Result/%sfragment.root",inp_sec_track.Data());
    std::cout << "Using this TOI: " << above << std::endl;
    TFile *aboveFile = TFile::Open(above.Data());
    auto fragm  = (TH1D*)aboveFile->Get("fraction_spline");
    return fragm;
}
TH1D *LoadL1ChargeCut(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    TString fileName = "../IonsSelected/"+getIonPath(charge)+"/Unfolding/"+getIonPath(charge)+
                       Form("_UnfoldingFactor%s_%s.root",inp_sec_track.Data(),timePeriod.Data());
    TFile *file = TFile::Open(fileName.Data());
    TH1D *l1ch = (TH1D*)file->Get("final_l1ch");
    return l1ch;
}


int main(int argc, char *argv[]) {
	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <time> <sec_track> \n", argv[0]);
	    return 1;
	}
    unsigned int charge=atoi(argv[1]);
    TString timePeriod = argv[2];
    TString sec_track = argv[3];
    TString inp_sec_track;
    sec_track == "y" ? inp_sec_track = "_sec_track_" : "_";
    TH1D *flux  = BuildFlux(charge,timePeriod,inp_sec_track);
    TH1D* flux_with_background = BuildFluxWithBackground(charge,timePeriod,inp_sec_track);
    TH1D* raw_flux = BuildRawFlux(charge,timePeriod,inp_sec_track);
    TH1D *fluxMultiplied= (TH1D *)flux->Clone();
    TH1D* fluxMultiplied_with_background = (TH1D*)flux_with_background->Clone();
    MultiplyByXPower(fluxMultiplied_with_background,2.70);
    MultiplyByXPower(fluxMultiplied,2.70);
    TString out = StreamUtility::getOutputDir(charge,argv[0],getIonName(charge)+
                  Form("_Flux%s_%s.root",inp_sec_track.Data(),timePeriod.Data()) );
    TFile *f = new TFile(out.Data(),"RECREATE");
    f->WriteTObject(flux,"flux");
    f->WriteTObject(fluxMultiplied,"fluxMultiplied");
    f->WriteTObject(flux_with_background,"flux_with_background");
    f->WriteTObject(fluxMultiplied_with_background,"fluxMultiplied_with_background");
    f->WriteTObject(raw_flux,"raw_flux");
    auto unf_fac  = LoadUnfoldingFactor(charge,timePeriod);
    auto tot_corr = LoadTotalCorrection(charge,timePeriod);
    auto purity   = LoadPurity(charge);
    auto toi      = LoadToi(charge,inp_sec_track);
    auto final_acc = LoadUnfoldedAcceptance(charge,timePeriod,inp_sec_track);
    auto l1ch = LoadL1ChargeCut(charge,timePeriod,inp_sec_track);
    auto counts = LoadCountsWithBackground(charge,timePeriod,inp_sec_track);
    auto time = LoadLT(charge,timePeriod);
    auto raw_acc = LoadRawAcceptance(charge,timePeriod,inp_sec_track);
    f->WriteTObject(unf_fac,"unf_fac");
    f->WriteTObject(tot_corr,"tot_corr");
    f->WriteTObject(purity,"purity");
    f->WriteTObject(toi,"toi");
    f->WriteTObject(final_acc,"final_acc");
    f->WriteTObject(l1ch,"l1ch");
    f->WriteTObject(counts,"counts");
    f->WriteTObject(time,"time");
    f->WriteTObject(raw_acc,"raw_acc");
    f->Close();
    return 1;
}