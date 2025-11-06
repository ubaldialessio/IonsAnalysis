#include "binning.h"
#include "utils.h"
#include "SplineUtility.h"

std::vector<TString> acceptances = {"acc", "spline_acc"};

TString getEffName(int j) {
    switch (j) {
        case 0: return "Mc acceptance";
        case 1: return "Spline mc acceptance";
        default: return "Unknown";
    }
}

std::vector<TH1D*> getAcceptanceMC(unsigned int charge) {
    std::vector<TH1D*> v;
    TString ionPath = getIonPath(charge);
    TFile *file = new TFile("../IonsSelected/" + ionPath + "/Passed/passed.root");

    auto mc_pass_gen = (TH1D*)file->Get("mc_pass_gen");
    auto mc_samp = (TH1D*)file->Get("mc_samp");
    mc_pass_gen->Rebin(9);
    mc_samp->Rebin(9);

    auto acc = divide(mc_pass_gen, mc_samp, "efficiency");
    acc->Scale(TMath::Pi() * 3.9 * 3.9);

    auto acc_range = SetFitLimits(SplineUtility::Efficiency::Acc, charge);
    auto acc_knot  = SetKnots(SplineUtility::Efficiency::Acc, charge);
    auto spline_acc = autospline(acc, 2.15, acc_range.second, acc_knot.first , acc_knot.second);
    //acc->Rebin(9);
    //acc->Scale(1./9);

    TFile *ff = new TFile("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/"+ ionPath+"_mcAcc.root","RECREATE");
    ff->WriteTObject(spline_acc,"spline_acc");
    ff->Close();

    v.push_back(acc);
    v.push_back(spline_acc);

    return v;
}

void setCorrectionColor(TH1D* h, int i) {
    if (i == 0) {
        h->SetMarkerColor(kRed-4);
        h->SetLineColor(kRed-4);
        }
    if (i == 1) {
         h->SetLineColor(kBlack);
         h->SetMarkerColor(kBlack);
    }
}

void print(std::vector<std::vector<TH1D*>> allAccs, std::vector<unsigned int> charges) {
    TStyle *style = effHistStyle(0);
    gROOT->SetStyle("MyStyle");
    gROOT->ForceStyle();
    setLogon();
    style->cd();
    gStyle->SetOptStat(0);

    TString prefix = getPrefixFile(charges) + "_mcAcc";

    TCanvas *c = new TCanvas("c", "", 2048, 1280);
    c->SaveAs(Form("../output/%s.pdf[", prefix.Data()), "RECREATE");
    TLegend *leg = new TLegend(0.35,0.9,0.9,0.65);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(4);

    for (int i = 0; i < charges.size(); ++i) {
        c->cd(i + 1);
        c->Update();
        gPad->SetLogx();
        auto ch = charges[i];
        /*TH2D *frame = new TH2D("", "", nRbins - 1, Rbins, 200, 0, 0.06);
        setTitle(frame, "", "R (GV)", "Acceptance (m^{2}sr)", ch);
        formatAxis(frame, charges.size());
        formatTitle(frame, charges.size());
        adjustZoomY(frame, "Mc acceptance", ch);
        adjustZoomX(frame, "Mc acceptance");
        frame->Draw();*/
        
        for (int j = 0; j < allAccs[i].size(); ++j) {
            TString effName = getEffName(j);
            auto hist = allAccs[i][j];
            if (j == 0) leg->AddEntry(hist, getIonName(ch));
            if (j == 1) {
                const char* input_cstr = hist->GetTitle();
                TString input(input_cstr);
                TString result;
                int chi2Pos = input.Index("Chi2/ndf:");
                if (chi2Pos != kNPOS) {
                    TString chi2Part = input(chi2Pos + 9, input.Length() - chi2Pos - 9);
                    int commaPos = chi2Part.Index(",");
                    TString chi2ValueStr = TString(chi2Part(0, commaPos));
                    chi2ValueStr.Strip(TString::kBoth);
                    double chi2Value = chi2ValueStr.Atof();
                    result.Form("Chi2/ndf = %.2f", chi2Value);
                }
                leg->AddEntry(hist, result, "L");
            }
            setTitle(hist, effName, "R (GV)", "", ch);		
            formatAxis(hist, 1);			
            formatTitle(hist,1);		
			adjustZoomX(hist,"Mc acceptance");
			adjustZoomY(hist,"Mc acceptance",ch);
            formatMarkerSize(hist,1);
            hist->GetYaxis()->SetTitle("Acceptance (m^{2} sr)");
            hist->GetYaxis()->SetTitleOffset(1.1);
            hist->GetYaxis()->SetLabelSize(0.045);
            hist->SetTitle(effName);
            hist->SetFillStyle(0);
            hist->SetLineWidth(2);
            //getColor(hist,ch, "dat");	
            if (j == 0) {
                getColor(hist, ch, "dat");
                hist->Draw("hist LP SAME");
            } if (j==1) {
                getColor(hist, ch, "mc");	
                hist->Draw("hist L SAME");
            }
        }
        leg->SetTextSize(0.022);
        leg->Draw("SAME");
    }

    c->SaveAs(Form("../output/%s.pdf", prefix.Data()));
    c->SaveAs(Form("../output/%s.pdf]", prefix.Data()));
}

// --------------------------- MAIN ----------------------------
int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: \n");
        printf("%s <charge1> <charge2> ... <charge n> \n", argv[0]);
        return 1;
    }

    std::vector<unsigned int> charges;
    std::vector<std::vector<TH1D*>> allAccs;

    for (int i = 1; i < argc; ++i) {
        unsigned int charge = atoi(argv[i]);
        charges.push_back(charge);
        allAccs.push_back(getAcceptanceMC(charge));
    }

    print(allAccs, charges);
    return 0;
}
