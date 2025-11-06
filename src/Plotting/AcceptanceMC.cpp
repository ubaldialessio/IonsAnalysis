#include "binning.h"
#include "utils.h"
#include "SplineUtility.h"

std::vector<TString> acceptances = {"acc","spline_acc"};
TString getEffName(int j) {
	TString A;
	switch(j) {
		case 0:
		A="Mc acceptance"; 
		break;
		case 1:
		A="Spline mc acceptance";
		break;
	}
	return A;
}
std::vector<TH1D*> getAcceptanceMC(unsigned int charge) {
    std::vector<TH1D*> v= {};
    TString ionPath = getIonPath(charge);
	TFile *file = new TFile("../IonsSelected/"+ionPath+"/RawAcceptance/"+"rawacc.root");
    auto mc_pass_gen = (TH1D*)file->Get("mc_pass_gen");
    auto mc_samp = (TH1D*)file->Get("mc_samp");
    mc_pass_gen->Rebin(9);
    mc_samp->Rebin(9);
    auto acc = divide(mc_pass_gen, mc_samp, "efficiency");
    acc->Scale(TMath::Pi() * 3.9 * 3.9);
    auto acc_range= SetFitLimits(SplineUtility::Efficiency::Acc,charge);
    auto spline_acc = autospline(acc, acc_range.first, acc_range.second,5,20);
    v.push_back(acc);
    v.push_back(spline_acc);
    return v;
}
void setCorrectionColor(TH1D* h, int i) {
    if (i == 0) h->SetLineColor(kRed);
    if (i == 1) h->SetLineColor(kBlue);
}

void print(std::vector<TH1D*> v, unsigned int charge) {
    TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
    std::vector<unsigned int> chNum;
    chNum.push_back(charge);
    TString prefix = getPrefixFile(chNum);
    TCanvas *b = new TCanvas("","",2048,1280);
    b->SaveAs( Form("../output/%s_mcAcc.pdf[", prefix.Data()) ,"RECREATE");
    auto a = new TLegend(0.70,0.9,0.9,0.75);
    b->cd();		
    b->Update();		
	TH2D *h = new TH2D("", "", nRbins-1, Rbins, 200, 0., 0.06);						
	setTitle(h, "", "R (GV)", "Mc acceptance", charge);			
	formatAxis(h, chNum.size());			
	formatTitle(h,chNum.size());			
    adjustZoomY(h, "Mc acceptance",charge);
    adjustZoomX(h, "Mc acceptance");
	h->Draw("");

    for (int i=0; i<v.size(); i++) {
        TString effName = getEffName(i);
        if (i==0)  a->AddEntry(v[i],effName,"L");
        if (i==1){ 
            const char* input_cstr = v[i]->GetTitle();
            TString input(input_cstr);  // Convert to TString
            TString result;  // Declare your output string
            int chi2Pos = input.Index("Chi2/ndf:");
            if (chi2Pos != kNPOS) {
                TString chi2Part = input(chi2Pos + 9, input.Length() - chi2Pos - 9);
                int commaPos = chi2Part.Index(",");
                TString chi2ValueStr = TString(chi2Part(0, commaPos));  // Convert TSubString to TString
                chi2ValueStr.Strip(TString::kBoth);  // Strip whitespace

                double chi2Value = chi2ValueStr.Atof();
                result.Form("Chi2/ndf: %.2f", chi2Value);  // Format
            }
            a->AddEntry(v[i],"Spline: "+result,"L");
        }
        setTitle(v[i], effName, "R (GV)", "Mc acceptance", chNum[i]);		
	    formatAxis(v[i], chNum.size());			
		formatTitle(v[i],chNum.size());		
        getColor(v[i], chNum[i], "dat");

        v[i]->SetTitle(effName);
        v[i]->SetFillStyle(0);
        v[i]->SetLineWidth(3);
        // Assign custom colors
        setCorrectionColor(v[i],i);	
        if (i==0) {
            v[i]->SetMarkerColor(kRed);
            v[i]->Draw("hist LP SAME");
        }
        if (i==1)
            v[i]->Draw("hist L SAME");
    }
    b->Update();
    a->Draw("SAME");
    a->SetTextSize(0.02);
    gPad->SetLogx();	
    b->SaveAs( Form("../output/%s_mcAcc.pdf", prefix.Data()) );
    b->SaveAs( Form("../output/%s_mcAcc.pdf]", prefix.Data()) );
}

//-------------------------MAIN----------------------------
int main(int argc, char **argv) {
	if (argc < 2) {
		printf("Usage: \n");
		printf("%s <charge> \n", argv[0]);
		return 1;
	}
    unsigned int charge=atoi(argv[1]);
    auto acc = getAcceptanceMC(charge);
    print(acc, charge);
}