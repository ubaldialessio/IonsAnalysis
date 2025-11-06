#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"
#include "SplineUtility.h"

TH1D* getAcceptance(unsigned int charge);
TString McFileName(unsigned int charge);
std::pair<int,int> lvtBasedOnQYan(unsigned int charge);

int main(int argc, char **argv) {
    TH1::SetDefaultSumw2();
	if (argc < 2) {
		printf("Usage: \n");
		printf("%s <charge> <time> \n", argv[0]);
		return 1;
	}
    unsigned int charge=atoi(argv[1]);
    TString timePeriod = argv[2];
    auto acc = getAcceptance(charge+1);
    auto acc_knots =SetKnots(SplineUtility::Efficiency::Acc,charge+1);
    auto spline_acc = autospline(acc,2.15,2000,acc_knots.first,acc_knots.second);
    auto flux = RebinHistogramWeighted(GetQYanFlux(charge+1));
    MultiplyByXPower(flux,-2.7);
    auto hwidth = (TH1D*)hist_rig_highZ->Clone();
	for(int i=1; i<=hist_rig_highZ->GetNbinsX(); ++i) {
		hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
		hwidth->SetBinError(i, 0);
	}
    //Getting livetime
    TString fileName= "../IonsSelected/Si/Livetime/livetime.root";
    TFile *file = TFile::Open(fileName.Data());
    int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
    auto lvt = (TH1D *)((TH2D *)file->Get("lvt_25"))->ProjectionY("q", start_bin, stop_bin);
    auto nBkg = (TH1D*)spline_acc->Clone();
    nBkg->Multiply(lvt);
    nBkg->Multiply(hwidth);
    nBkg->Multiply(flux);

    //for the nucleus under study, get the counts from data
    TFile *file2 = TFile::Open("../IonsSelected/"+getIonPath(charge)+"/Counts/counts.root");
    auto counts = (TH1D *)((TH2D *)file2->Get("rigidity"))->ProjectionY("counts", start_bin, stop_bin);

    TH1D* p = (TH1D*)hist_rig_highZ->Clone();
    for (int i=1; i<= p->GetNbinsX(); i++) {
        p->SetBinContent(i, 1-(nBkg->GetBinContent(i) )/(counts->GetBinContent(i)) );
    }

    TCanvas *b = new TCanvas("","",2048,1280);
    b->SaveAs( Form("../output/%s_mcPurity_%s.pdf[",getIonPath(charge).Data() ,timePeriod.Data()) ,"RECREATE");	
    p->GetXaxis()->SetRangeUser(1.5,1000);
    b->SetLogx();	
    p->Draw();
    b->SaveAs( Form("../output/%s_mcPurity_%s.pdf", getIonPath(charge).Data() ,timePeriod.Data()) ,"RECREATE");		
    b->SetLogy();
    acc->Draw();
    spline_acc->Draw("SAME");
    b->SaveAs( Form("../output/%s_mcPurity_%s.pdf", getIonPath(charge).Data() ,timePeriod.Data()) ,"RECREATE");		
    nBkg->Draw();
    b->SaveAs( Form("../output/%s_mcPurity_%s.pdf", getIonPath(charge).Data() ,timePeriod.Data()) ,"RECREATE");		
    counts->Draw();
    b->SaveAs( Form("../output/%s_mcPurity_%s.pdf", getIonPath(charge).Data() ,timePeriod.Data()) ,"RECREATE");	
    b->SaveAs( Form("../output/%s_mcPurity_%s.pdf]",getIonPath(charge).Data() ,timePeriod.Data()) ,"RECREATE");		
}

TH1D* getAcceptance(unsigned int charge) {
    TString path = "../IonsSelected/P/Fragments/June/"+McFileName(charge)+".root";
    TFile *file = TFile::Open(path.Data(), "READ");
    TH1D *mc_pass = (TH1D*)file->Get("mc_pass");
    TH1D *mc_samp = (TH1D*)file->Get("mc_samp");
    mc_pass->Sumw2();
    mc_samp->Sumw2();
    mc_pass->Rebin(9);
    mc_samp->Rebin(9);
    auto h = divide(mc_pass,mc_samp,"efficiency");
    h->Scale(TMath::Pi() * 3.9 * 3.9,"nosw2");
    return h;
}
TString McFileName(unsigned int charge) {
    TString a;
    switch(charge) {
        case 14:
            a = "Si.B1236";
            break;
        case 15:
            a = "P.B1236";
            break;
        case 16:
            a = "S.B1236";
            break;
        case 17:
            a = "Cl.B1236";
            break;
        case 18:
            a = "Ar.B1236";
            break;
        case 19:
            a = "K.B1236";
            break;
        case 20:
            a = "Ca.B1236";
            break;
        case 26:
            a = "Fe.B1236";
            break;
    }
    return a;
}
std::pair<int,int> lvtBasedOnQYan(unsigned int charge) {
    switch(charge) {
        case 14:
        return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1620259200) ); //10y
        break;
        case 15:
        return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1654087477) ); //11y
        break;              
        case 16:
        return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1620259200) ); //10y
        break;
        case 17:
        return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1654087477) ); //11y
        break;
        case 18:
        return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1654087477) ); //11y
        break;
        case 19:
        return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1654087477) ); //11y
        break;
        case 20:
        return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1654087477) ); //11y
        break;
        case 26:
        return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1620259200) ); //10y
        break;
    }
    return std::make_pair(-1,-1);
}