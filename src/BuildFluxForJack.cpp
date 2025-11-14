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
TH1D *LoadRawAcceptance(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting Final Acceptance
    TString fileName= Form("../IonsSelected/"+getIonPath(charge)+"/Passed/passed.root");
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
TH1D *LoadUnfoldedAcceptance(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting Final Acceptance
    TString fileName = "../IonsSelected/"+getIonPath(charge)+"/Unfolding/"+getIonPath(charge)+
                       Form("_UnfoldingFactor%s_%s.root",inp_sec_track.Data(),timePeriod.Data());
    TFile *file = TFile::Open(fileName.Data());
    TH1D *acc = (TH1D*)file->Get("final_acc");
    return acc;
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
TH1D *LoadCorrectionNoTrigger(unsigned int charge, TString timePeriod, TH1D* rawacc) {
    auto corrAcc = (TH1D*)rawacc->Clone("corrAcc");
    TString ionPath = getIonPath(charge);
	TFile *file = new TFile("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+Form("_UnfoldingFactor_%s.root",timePeriod.Data()) );
    for (auto i : fit) {
        if (i=="final_damc_tr") continue;
        auto hist = (TH1D*)file->Get(i.Data());
        printf("Retrieving ISS/MC corrections : %s \n",i.Data());
        if (hist) corrAcc->Multiply(hist);
    }
    return corrAcc;
}
TH1D *LoadTriggerCorrection(unsigned int charge, TString timePeriod) {
    TString ionPath = getIonPath(charge);
	TFile *file = new TFile("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+Form("_UnfoldingFactor_%s.root",timePeriod.Data()) );
    for (auto i : fit) {
        if (i=="final_damc_tr") {
            auto hist = (TH1D*)file->Get(i.Data());
            printf("Retrieving ISS/MC corrections : %s \n",i.Data());
            return hist;
        }
    }
    return nullptr; 
}
TF1 *LoadPurity(unsigned int charge) {
    TString contName= "../Fragmentation/BelowL1/Fractions/contamination_"+getIonName(charge)+".root";
    TFile *contFile = TFile::Open(contName.Data());
    auto purity = (TF1*)contFile->Get("purityFit");
    return purity;
}
TH1D *LoadToi(unsigned int charge,TString inp_sec_track, TString ACC, TString nucleus) {
    TString above = Form("../IonsSelected/"+getIonPath(charge)+"/Fragments/Result/%sfragment_%s.root",
                         inp_sec_track.Data(), ACC.Data() );
    TFile *aboveFile = TFile::Open(above.Data());
    if (nucleus=="fraction") {
        auto fragm  = (TH1D*)aboveFile->Get( Form("fraction_spline_%s",ACC.Data()) );
        return fragm;
    } if (nucleus=="total") {
        auto fragm  = (TH1D*)aboveFile->Get( Form("fraction_spline_tot%s",ACC.Data()) );
        return fragm;
    } else {
        auto fragm  = (TH1D*)aboveFile->Get( Form("single_spline_%s_%s",nucleus.Data(),ACC.Data()) );
        return fragm; 
    }
}
TH1D *LoadToiAcceptance(unsigned int charge,TString inp_sec_track, TString ACC, TString nucleus) {
    TString above = Form("../IonsSelected/"+getIonPath(charge)+"/Fragments/Result/%sfragment_%s.root",
                         inp_sec_track.Data(), ACC.Data() );
    TFile *aboveFile = TFile::Open(above.Data());
    auto fragm  = (TH1D*)aboveFile->Get( Form("spl_%s_%s",nucleus.Data(),ACC.Data()) );
    return fragm; 
}
TH1D *LoadL1ChargeCut(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    TString fileName = "../IonsSelected/"+getIonPath(charge)+"/Unfolding/"+getIonPath(charge)+
                       Form("_UnfoldingFactor%s_%s.root",inp_sec_track.Data(),timePeriod.Data());
    TFile *file = TFile::Open(fileName.Data());
    TH1D *l1ch = (TH1D*)file->Get("final_l1ch");
    return l1ch;
}
TTree *LoadTree(unsigned int charge, TString inp_sec_track, TString ACC, TString nucleus) {
    TString above = Form("./../IonsSelected/%s/Fragments/%s.B1236.root",
                         getIonPath(charge).Data(), nucleus.Data());
    std::cout << "Getting tree from file " << above.Data() << std::endl;

    TFile *aboveFile = TFile::Open(above.Data(), "READ");
    if (!aboveFile || aboveFile->IsZombie()) {
        std::cerr << "Error opening file " << above << std::endl;
        return nullptr;
    }

    TTree *sourceTree = (TTree*)aboveFile->Get(Form("events_%s", ACC.Data()));
    if (!sourceTree) {
        std::cerr << "Error: cannot find TTree events_" << ACC << std::endl;
        aboveFile->Close();
        delete aboveFile;
        return nullptr;
    }

    // ✅ Carichiamo almeno un evento (forza i basket)
    sourceTree->GetEntry(0);

    // ✅ Copiamo tutto in memoria
    TTree *ttree = sourceTree->CopyTree("");
    ttree->SetDirectory(nullptr);

    aboveFile->Close();
    delete aboveFile;

    return ttree;
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
    sec_track == "y" ? inp_sec_track = "_sec_track_" : "";
    TString out = StreamUtility::getOutputDir(charge,argv[0],getIonName(charge)+
                  Form("_FluxForJack%s_%s.root",inp_sec_track.Data(),timePeriod.Data()) );
    TFile *f = new TFile(out.Data(),"RECREATE");

    auto flux  = BuildFlux(charge,timePeriod,inp_sec_track);
    f->WriteTObject(flux, Form("Z%d_flux",charge) );

    auto time = LoadLT(charge,timePeriod);
    f->WriteTObject(time,"exposure_time");

    auto event_counts = LoadCountsWithBackground(charge,timePeriod,inp_sec_track);
    f->WriteTObject(event_counts, Form("Z%d_raw_events",charge) );

    auto raw_acc = LoadRawAcceptance(charge,timePeriod,inp_sec_track);
    f->WriteTObject(raw_acc, Form("Z%d_acceptance_raw",charge) );

    auto corr_acc = LoadCorrectionNoTrigger(charge, timePeriod,raw_acc);
    f->WriteTObject(corr_acc, Form("Z%d_acceptance_corrected",charge) );

    auto trigger = LoadTriggerCorrection(charge, timePeriod);
    f->WriteTObject(trigger, Form("Z%d_trigger_efficiency",charge) );

    auto toi_acc4      = LoadToi(charge,inp_sec_track,"acc4","fraction");
    f->WriteTObject(toi_acc4, Form("Z%d_al1_total_acc4",charge) );

    auto toi_acc7      = LoadToi(charge,inp_sec_track,"acc7","fraction");
    f->WriteTObject(toi_acc7, Form("Z%d_al1_total_acc7",charge) );

    auto toi_S_acc4      = LoadToi(charge,inp_sec_track,"acc4","S");
    f->WriteTObject(toi_S_acc4, Form("Z%d_al1_sulfur_acc4",charge) );

    auto toi_S_acc7      = LoadToi(charge,inp_sec_track,"acc7","S");
    f->WriteTObject(toi_S_acc7, Form("Z%d_al1_sulfur_acc7",charge) );

    auto toi_Fe_acc4      = LoadToi(charge,inp_sec_track,"acc4","Fe");
    f->WriteTObject(toi_Fe_acc4, Form("Z%d_al1_iron_acc4",charge) );

    auto toi_Fe_acc7      = LoadToi(charge,inp_sec_track,"acc7","Fe");
    f->WriteTObject(toi_Fe_acc7, Form("Z%d_al1_iron_acc7",charge) );
     
    auto toi_S_acc_acc4 = LoadToiAcceptance(charge,inp_sec_track,"acc4","S");
    f->WriteTObject(toi_S_acc_acc4, Form("Z16_to_Z15_acceptance_acc4") );

    auto toi_S_acc_acc7 = LoadToiAcceptance(charge,inp_sec_track,"acc7","S");
    f->WriteTObject(toi_S_acc_acc7, Form("Z16_to_Z15_acceptance_acc7") );

    auto toi_Fe_acc_acc4 = LoadToiAcceptance(charge,inp_sec_track,"acc4","Fe");
    f->WriteTObject(toi_Fe_acc_acc4, Form("Z26_to_Z15_acceptance_acc4") );

    auto toi_Fe_acc_acc7 = LoadToiAcceptance(charge,inp_sec_track,"acc7","Fe");
    f->WriteTObject(toi_Fe_acc_acc7, Form("Z26_to_Z15_acceptance_acc7") );

    auto purity   = LoadPurity(charge);
    f->WriteTObject(purity, Form("Z%d_bl1_purity",charge) );

    auto l1ch = LoadL1ChargeCut(charge,timePeriod,inp_sec_track);
    f->WriteTObject(l1ch, Form("Z%d_l1q_efficiency",charge) );

    auto unf_fac  = LoadUnfoldingFactor(charge,timePeriod);
    f->WriteTObject(unf_fac, Form("Z%d_unfolding_factor",charge) );

    f->Close();

    //event and count numbers
    TString outNumbers_26 = Form("./../IonsSelected/%s/Flux/TOI_26_to_15.root",getIonPath(charge).Data());
    TString outNumbers_16 = Form("./../IonsSelected/%s/Flux/TOI_16_to_15.root",getIonPath(charge).Data());
    TFile *f_26 = new TFile(outNumbers_26.Data(),"RECREATE");
    TFile *f_16 = new TFile(outNumbers_16.Data(),"RECREATE");

    auto tree_S_acc4 = LoadTree(charge,inp_sec_track,"acc4","S");
    f_16->WriteTObject(tree_S_acc4, Form("events_acc4") );
    auto tree_S_acc7 = LoadTree(charge,inp_sec_track,"acc7","S");
    f_16->WriteTObject(tree_S_acc7, Form("events_acc7") );

    auto tree_Fe_acc4 = LoadTree(charge,inp_sec_track,"acc4","Fe");
    f_26->WriteTObject(tree_Fe_acc4, Form("events_acc4") );
    auto tree_Fe_acc7 = LoadTree(charge,inp_sec_track,"acc7","Fe");
    f_26->WriteTObject(tree_Fe_acc7, Form("events_acc7") );

    f_16->Close();
    f_26->Close();
    return 1;
}