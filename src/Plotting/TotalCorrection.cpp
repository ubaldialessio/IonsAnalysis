#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

std::vector<TString> corr = {"final_damc_l1","final_damc_tf","final_damc_tr","final_damc_tk",
                             "final_damc_l1u","final_damc_tkch","final_daq","final_l1ch","final_damc_tot"};
TString getEffName(int j) {
	TString A;
	switch(j) {
		case 0:
            A="L1 pickup"; 
            break;
		case 1:
            A="Tof";
            break;
		case 2:
            A="Trigger";
            break;
		case 3:
            A="Track";
            break;
        case 4:
            A="L1 detection";
            break;
        case 5:
            A="Track charge";
            break;
        case 6:
            A="Daq";
            break;
        case 7:
            A="L1 charge";
            break;
        case 8:
            A ="Total";
            break;
	}
	return A;
}
std::vector<TH1D*> getCorrections(unsigned int charge, TString timePeriod) {
    std::vector<TH1D*> v= {};
    TString ionPath = getIonPath(charge);
	TFile *file = new TFile("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+Form("_UnfoldingFactor_%s.root",timePeriod.Data()) );
    for (auto i : corr) {
        //per il Fosforo non prendere il total
        //if (charge == 15 && i == "final_damc_tot") continue;
        auto hist = (TH1D*)file->Get(i.Data());
        printf("Retrieving histogram : %s \n",i.Data());
        if (hist) v.push_back(hist);
    }
    //produci il total per il fosforo
    /*if (charge == 15) {
        auto tot = (TH1D*)v[0]->Clone();
        for (int j=1; j<v.size(); j++)
            tot->Multiply(v[j]);
        v.push_back(tot);
    }*/
    return v;
}

void print(std::vector<TH1D*> v, unsigned int charge, TString timePeriod) {
    TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
    std::vector<unsigned int> chNum;
    chNum.push_back(charge);
    TString prefix = getPrefixFile(chNum);
    TCanvas *b = new TCanvas("","",2048,1280);
    b->SaveAs( Form("../output/%s_totCorr_%s.pdf[", prefix.Data(), timePeriod.Data()) ,"RECREATE");
    auto a = new TLegend(0.15,0.9,0.6,0.60);
    b->cd();		
    b->Update();		
	TH2D *h = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 2.15);						
	setTitle(h, "", "R (GV)", "Data/Mc corrections", charge);			
	formatAxis(h, chNum.size());			
	formatTitle(h,chNum.size());			
    adjustZoomY(h, "Total Data/Mc correction overlap", charge);
    h->GetXaxis()->SetRangeUser(1.92,1000);
    h->GetYaxis()->SetNdivisions(18);
	h->Draw("");
    for (int i=0; i<v.size(); i++) {
        TString effName = getEffName(i);
        a->AddEntry(v[i],effName,"L");
        setTitle(v[i], effName, "R (GV)", "Data/Mc ratio", chNum[i]);		
	    formatAxis(v[i], chNum.size());			
		formatTitle(v[i],chNum.size());		
        getColor(v[i], chNum[i], "dat");	
        v[i]->SetTitle(effName);
        v[i]->SetFillStyle(0);
        v[i]->SetLineWidth(3.8);
        // Assign custom colors
        setFriendColor(v[i],i);	
        printf("Drawing histogram : %s \n",effName.Data());
        v[i]->Draw("hist L SAME");
    }
    a->Draw("SAME");
    a->SetTextSize(0.038);
    a->SetTextFont(62);
    a->SetBorderSize(0);
    a->SetFillStyle(0);
    gPad->SetLogx();	
    b->SaveAs( Form("../output/%s_totCorr_%s.pdf", prefix.Data(), timePeriod.Data()) );
    b->SaveAs( Form("../output/%s_totCorr_%s.pdf]", prefix.Data(), timePeriod.Data()) );
}

//-------------------------MAIN----------------------------
int main(int argc, char **argv) {
	if (argc < 3) {
		printf("Usage: \n");
		printf("%s <time> <charge> \n", argv[0]);
		return 1;
	}
    TString timePeriod = argv[1];
    unsigned int charge=atoi(argv[2]);
    auto corrections = getCorrections(charge, timePeriod);
    print(corrections, charge, timePeriod);
    TString out = Form("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/%s_TotalCorrections%s_%s.root",
                        getIonPath(charge).Data(),"",timePeriod.Data());
    TFile *f = new TFile(out.Data(),"RECREATE");
    std::cout << out << std::endl;
    f->cd();
    for (int i=0; i<corrections.size(); i++) {
        std::cout << "Writing object " << corr[i].Data() << std::endl;
        f->WriteTObject(corrections[i], corr[i].Data());
    }
    f->Write();
    f->Close();
}