#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "SplineUtility.h"


std::vector<TString> finalNames= {"damc_l1","damc_tf","damc_tr","damc_tk","damc_l1u","damc_tkch","",
    "final_damc_l1","final_damc_tf","final_damc_tr","final_damc_tk","final_damc_l1u","final_damc_tkch","final_daq"};
std::vector<TString> correctionNames= {"damc_l1","damc_tf","damc_tr","damc_tk","damc_l1u","damc_tkch"};
std::vector<TString> splineNames= {"final_damc_l1","final_damc_tf","final_damc_tr","final_damc_tk","final_damc_l1u","final_damc_tkch","final_daq"};

std::vector<TH1D*> Interpolate(std::vector<std::pair<TH1D*,TH1D*>> v, TString timePeriod,TString inp_sec_track);
std::vector<TH1D*> InterpolateWeightedAverage(std::vector<std::pair<TH1D*,TH1D*>> v, TString timePeriod);
std::vector<std::pair<TH1D*,TH1D*>> getCorrections(unsigned int charge1, unsigned int charge2, TString timePeriod,TString inp_sec_track);
void appendToOutput(std::vector<TH1D*> v, unsigned int charge, TString timePeriod, TString inp_sec_track);

int main(int argc, char *argv[]) {
	if (argc < 4) {
		printf("Usage: \n");
	    printf("%s <time> <charge 1> <charge 2> <sec_track> \n", argv[0]);
	    return 1;
	}
    TString timePeriod = argv[1];
    unsigned int charge1=atoi(argv[2]), charge2=atoi(argv[3]);
    TString sec_track = argv[4];
    TString inp_sec_track;
    sec_track == "y" ? inp_sec_track = "_sec_track_" : "";
    auto corrections = getCorrections(charge1,charge2,timePeriod,inp_sec_track);
    auto v = Interpolate(corrections,timePeriod,inp_sec_track);
    appendToOutput(v,15,timePeriod,inp_sec_track);
    return 1;
}

std::vector<std::pair<TH1D*,TH1D*>> getCorrections(unsigned int charge1, unsigned int charge2, TString timePeriod,TString inp_sec_track) {
    std::vector<std::pair<TH1D*,TH1D*>> v = {};
    TString ionPath1 = getIonPath(charge1);
    TString ionPath2 = getIonPath(charge2);
    TFile *f1 = new TFile("../IonsSelected/"+ionPath1+"/Unfolding/"+ionPath1+
                          Form("_UnfoldingFactor_%s.root",timePeriod.Data()));
    TFile *f2 = new TFile("../IonsSelected/"+ionPath2+"/Unfolding/"+ionPath2+
                          Form("_UnfoldingFactor_%s.root",timePeriod.Data()));
    for (const auto name : splineNames) {
        auto graph1 = (TH1D*)f1->Get(name.Data());
        auto graph2 = (TH1D*)f2->Get(name.Data());
        if (graph1 && graph2)  {
            printf("--- Retrieving histogram %s  \n",name.Data());
            auto p = std::make_pair(graph1,graph2);
            p.first->SetName(name);
            p.second->SetName(name);
            v.push_back(p);
        }
    }
    return v;
}
std::vector<TH1D*> Interpolate(std::vector<std::pair<TH1D*,TH1D*>> v, TString timePeriod,TString inp_sec_track) {
    std::vector<TH1D*> vv;
    for (const auto p : v) { //loop on the elements of the vector of the corrections
        TH1D *g = (TH1D*)hist_rig_highZ->Clone();
        for (int i = 1; i <= g->GetNbinsX(); ++i) {
            double x = g->GetBinCenter(i);

            int bin1 = p.first->FindBin(x);
            int bin2 = p.second->FindBin(x);

            double val1 = p.first->GetBinContent(bin1);
            double val2 = p.second->GetBinContent(bin2);
            double err1 = p.first->GetBinError(bin1);
            double err2 = p.second->GetBinError(bin2);

            double mean = (val1 + val2) / 2;
            double err = 0.5 * TMath::Sqrt(err1 * err1 + err2 * err2);

            g->SetBinContent(i, mean);
            g->SetBinError(i, err2);  //putting errors of Sulfur
        }
        vv.push_back(g);
    }
    //Producing also the total correction 
    auto final_damc_tot = new TH1D();
    final_damc_tot = (TH1D*)vv[0]->Clone();
    final_damc_tot->Multiply(vv[1]);
    final_damc_tot->Multiply(vv[2]);
    final_damc_tot->Multiply(vv[3]);
    final_damc_tot->Multiply(vv[4]);
    final_damc_tot->Multiply(vv[5]);
    final_damc_tot->Multiply(vv[6]);
    TString f = Form("../IonsSelected/P/Unfolding/P_UnfoldingFactor_%s%s.root",inp_sec_track.Data(),timePeriod.Data());
    auto ff = new TFile(f);
    auto spline_acc= (TH1D*)ff->Get("spline_acc");
    auto final_acc = (TH1D*)spline_acc->Clone();
    final_acc->Multiply(final_damc_tot);
    final_acc->SetName("final_acc");
    final_damc_tot->SetName("final_damc_tot");
    vv.push_back(vv[0]);
    vv.push_back(vv[1]);
    vv.push_back(vv[2]);
    vv.push_back(vv[3]);
    vv.push_back(vv[4]);
    vv.push_back(vv[5]);
    vv.push_back(vv[6]);
    vv.push_back(final_damc_tot);
    vv.push_back(final_acc);
    return vv;
}
std::vector<TH1D*> InterpolateWeightedAverage(std::vector<std::pair<TH1D*,TH1D*>> v, TString timePeriod) {
    std::vector<TH1D*> vv;  
    for (const auto p : v) { //loop on the elements of the vector of the corrections
        TH1D *g = (TH1D*)hist_rig_highZ->Clone();
        for (int i=1; i<=g->GetNbinsX(); ++i) { //loop on the common points
            double x_center = g->GetBinCenter(i);
            double value1 = p.first->Interpolate(x_center);
            double value2 = p.second->Interpolate(x_center);
            int bin1 = p.first->FindBin(x_center);
            int bin2 = p.second->FindBin(x_center);
            double error1 = p.first->GetBinError(bin1);
            double error2 = p.second->GetBinError(bin2);
            /*double value1 = p.first->GetBinContent(i);
            double error1 = p.first->GetBinError(i);
            double value2 = p.second->GetBinContent(i);
            double error2 = p.second->GetBinError(i);*/
            // Skip bins with zero or invalid errors
            if (error1 == 0 || error2 == 0) {
                std::cerr << "Error in bin " << i << " is zero or negative!" << std::endl;
                continue;
            }
            // Compute the weighted average for this bin
            double weight1 = 1.0 / (error1 * error1);
            double weight2 = 1.0 / (error2 * error2);
            double weightedAverage = (value1 * weight1 + value2 * weight2) / (weight1 + weight2);
            double weightedError = sqrt(1.0 / (weight1 + weight2));

            g->SetBinContent(i,weightedAverage);
            g->SetBinError(i,error2);
        }
        vv.push_back(g);
    }
    //Producing also the total correction 
    auto final_damc_tot = new TH1D();
    final_damc_tot = (TH1D*)vv[0]->Clone();
    final_damc_tot->Multiply(vv[1]);
    final_damc_tot->Multiply(vv[2]);
    final_damc_tot->Multiply(vv[3]);
    final_damc_tot->Multiply(vv[4]);
    final_damc_tot->Multiply(vv[5]);
    TString f = Form("../IonsSelected/P/Unfolding/P_UnfoldingFactor_%s.root",timePeriod.Data());
    auto ff = new TFile(f);
    auto spline_acc= (TH1D*)ff->Get("spline_acc");
    auto final_acc = (TH1D*)spline_acc->Clone();
    final_acc->Multiply(final_damc_tot);
    final_acc->SetName("final_acc");
    final_damc_tot->SetName("final_damc_tot");
    vv.push_back(vv[0]);
    vv.push_back(vv[1]);
    vv.push_back(vv[2]);
    vv.push_back(vv[3]);
    vv.push_back(vv[4]);
    vv.push_back(vv[5]);
    vv.push_back(final_damc_tot);
    vv.push_back(final_acc);
    return vv;
}
void appendToOutput(std::vector<TH1D*> v, unsigned int charge,TString timePeriod, TString inp_sec_track) {
    TString ionPath = getIonPath(charge);
    TString outPath = "../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+
                       Form("_UnfoldingFactor%s_%s.root",inp_sec_track.Data(),timePeriod.Data());
    TString intPath = "../IonsSelected/"+ionPath+"/Interpolation/"+ionPath+
                       Form("_interpolation%s_%s.root",inp_sec_track.Data(),timePeriod.Data());                   
    TString name;
    TFile *file = new TFile(outPath.Data(), "UPDATE");
    file->cd();
    for (int ii=0; ii<finalNames.size(); ii++) {
        v[ii]->SetName(finalNames[ii].Data() );
        name = finalNames[ii].Data();
    }
    for (const auto i : v) {
        printf("Writing histogram %s \n",i->GetName());
        file->Delete( Form("%s;*",i->GetName()) );
        if (i) i->Write();
    }
    file->Close();

    TFile *file1 = new TFile(intPath.Data(), "UPDATE");
    file1->cd();
    for (int ii=0; ii<finalNames.size(); ii++) {
        v[ii]->SetName(finalNames[ii].Data() );
        name = finalNames[ii].Data();
    }
    for (const auto i : v) {
        printf("Writing histogram %s \n",i->GetName());
        file1->Delete( Form("%s;*",i->GetName()) );
        if (i) i->Write();
    }
    file1->Close();

}