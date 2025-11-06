#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"
#include "SplineUtility.h"

std::vector<unsigned int> allowedZ = {14, 15, 16, 17, 18, 19, 20, 26};

std::vector<TH1D*> LoadLT(unsigned int charge, TString timePeriod);
TH1D *LoadRBinWidth();
std::vector<TH1D*> GetParentCounts(unsigned int charge,TString timePeriod);
std::map<std::string,std::vector<TH1D*>> GetParentAcceptances(unsigned int charge, std::vector<TString> &Parents,
                                        TString sec_track, TString l1_cut);
std::map<std::string,std::vector<TH1D*>> SplineParentAcceptances(std::map<std::string,std::vector<TH1D*>> acc,
                                        std::vector<TString> Parents, unsigned int charge);
std::vector<TH1D*> GetGalpropFluxes(unsigned int charge);
std::vector<TH1D*> GetQyanFluxes(unsigned int charge);
std::map<std::string,std::vector<TH1D*>>  Contamination(std::map<std::string,std::vector<TH1D*>> spline_acc,
                                                        std::vector<TH1D*> fluxes, 
                                 std::vector<TH1D*> lvt, TH1D* rwidth,std::vector<TString> Parents);
std::map<std::string,std::vector<TH1D*>> SingleContribution(std::vector<TH1D*> counts, TH1D* den);
std::map<std::string,TH1D*> GetSum(std::map<std::string,std::vector<TH1D*>>  v, int start);
void Save(const std::map<std::string, std::vector<TH1D*>> &m, TFile *f, TString ACC) {
    if (!f || !f->IsOpen()) {
        std::cerr << "Save(): invalid or closed TFile!" << std::endl;
        return;
    }

    // Look for the matching key in the map
    auto it = m.find(ACC.Data());
    if (it == m.end()) {
        std::cerr << "Save(): no entry found for key \"" << ACC << "\"" << std::endl;
        return;
    }

    // Write only histograms from the matching vector
    const auto &vec = it->second;
    for (const auto *obj : vec) {
        if (!obj) continue;
        f->WriteTObject(obj, Form("%s_%s", obj->GetName(), ACC.Data()));
    }
}
template <typename T>
void Save(const std::vector<T*> &v, TFile *f, const std::string &tag = "") {
    static_assert(std::is_base_of<TObject, T>::value, "T must inherit from TObject");

    if (!f || !f->IsOpen()) {
        std::cerr << "Save(): invalid or closed TFile!" << std::endl;
        return;
    }

    for (const auto *obj : v) {
        if (!obj) continue;
        if (tag.empty())
            f->WriteTObject(obj, obj->GetName());
        else
            f->WriteTObject(obj, Form("%s_%s", obj->GetName(), tag.c_str()));
    }
}
void PrintImg(std::vector<TH1D*> acc, std::vector<TH1D*> spl, std::vector<TH1D*> cont, TH1D *sum, unsigned int charge,
             std::vector<TH1D*> SingleContribution, TH1D* fraction, TString sec_track, TString l1_cut, TH1D* fraction_spline, TString ACC);
TString McFileName(unsigned int charge);
std::pair<int,int> lvtBasedOnQYan(unsigned int charge);
TH1D *buildDummyFlux(unsigned int charge, TString timePeriod, TString inp_sec_track);
std::map<std::string,TH1D*> WeightForToi(unsigned int charge, TString timePeriod, TString inp_sec_track);
TH1D* TotalWeightedToi(std::map<std::string,TH1D*> w, TH1D* f_acc4, TH1D* f_acc7);

int main(int argc, char **argv) {
    TH1::SetDefaultSumw2();
	if (argc < 5) {
		printf("Usage: \n");
		printf("%s <charge> -print <y/n> -sec track <y/n> -l1Cut <y/n> <time> \n", argv[0]);
		return 1;
	}
    unsigned int charge=atoi(argv[1]);
    TString print = argv[2];
    TString sec_track = argv[3];
    TString l1_cut = argv[4];
    TString timePeriod = argv[5];
    TString inp_sec_track;
    sec_track == "y" ? inp_sec_track = "_sec_track_" : "";
    TString out_acc3 = StreamUtility::getOutputDir(charge,argv[0], Form("%sfragment_acc4.root",inp_sec_track.Data()) );
    TString out_acc7 = StreamUtility::getOutputDir(charge,argv[0], Form("%sfragment_acc7.root",inp_sec_track.Data()) );
    TFile *f_acc4 = new TFile(out_acc3.Data(),"RECREATE");
    TFile *f_acc7 = new TFile(out_acc7.Data(),"RECREATE");
    std::vector<TString> Parents;
    auto ParentAcceptances = GetParentAcceptances(charge,Parents,sec_track,l1_cut);
    auto Counts = GetParentCounts(charge,timePeriod);
    auto SplineAcceptances = SplineParentAcceptances(ParentAcceptances,Parents,charge);
    //auto GalpropFluxes = GetGalpropFluxes(charge); //orignal
    auto QYanFluxes = GetQyanFluxes(charge);
    auto Contaminations = Contamination(SplineAcceptances,QYanFluxes,LoadLT(charge,timePeriod),LoadRBinWidth(),Parents);
    // auto den = GetSum(Contaminations); original
    auto den = Counts[0];
    auto num = GetSum(Contaminations,0);
    Save(ParentAcceptances,f_acc4,"acc4");
    Save(SplineAcceptances,f_acc4,"acc4");
    Save(Contaminations,f_acc4,"acc4");
    Save(ParentAcceptances,f_acc7,"acc7");
    Save(SplineAcceptances,f_acc7,"acc7");
    Save(Contaminations,f_acc7,"acc7");
    //Save(Counts,f);
    f_acc4->WriteTObject(num["acc4"],"num_acc4");
    f_acc4->WriteTObject(den,"den");
    f_acc7->WriteTObject(num["acc7"],"num_acc7");
    f_acc7->WriteTObject(den,"den");
    //Evaluation of the contamination fraction for each single period 
    auto fraction_acc4 = divide(num["acc4"],den,"ratio");
    auto fraction_acc7 = divide(num["acc7"],den,"ratio");
    auto weights = WeightForToi(charge,timePeriod,inp_sec_track);
    auto total_fraction = TotalWeightedToi(weights,fraction_acc4,fraction_acc7);
    //Evaluat single contributions
    //Contaminations.erase(Contaminations.begin());
    auto single_acc4 = SingleContribution(Contaminations["acc4"],den);
    auto single_acc7 = SingleContribution(Contaminations["acc7"],den);
    f_acc4->WriteTObject(fraction_acc4,"fraction_acc4");
    f_acc4->WriteTObject(total_fraction,"total_fraction");
    f_acc7->WriteTObject(fraction_acc7,"fraction_acc7");
    f_acc7->WriteTObject(total_fraction,"total_fraction");
    Save(single_acc4["point"],f_acc4,"acc4");
    Save(single_acc7["point"],f_acc7,"acc7");
    Save(single_acc4["spline"],f_acc4,"acc4");
    Save(single_acc7["spline"],f_acc7,"acc7");
    auto fraction_spline_acc4 = autospline(fraction_acc4,0.8,1000,2,10);
    auto fraction_spline_acc7 = autospline(fraction_acc7,0.8,1000,2,10);
    auto fraction_spline_tot  = autospline(total_fraction,0.8,1000,2,10);
    f_acc4->WriteTObject(fraction_spline_acc4,"fraction_spline_acc4");
    f_acc4->WriteTObject(fraction_spline_tot,"fraction_spline_tot");
    f_acc7->WriteTObject(fraction_spline_acc7,"fraction_spline_acc7");
    f_acc7->WriteTObject(fraction_spline_tot,"fraction_spline_tot");
    f_acc4->Close();
    f_acc7->Close();
    if (print=="y") {
        PrintImg(ParentAcceptances["acc4"],SplineAcceptances["acc4"],Contaminations["acc4"],
                 den,charge,single_acc4["point"],fraction_acc4,sec_track,l1_cut,fraction_spline_acc4,"acc4");

        PrintImg(ParentAcceptances["acc7"],SplineAcceptances["acc7"],Contaminations["acc7"],
                 den,charge,single_acc7["point"],fraction_acc7,sec_track,l1_cut,fraction_spline_acc7,"acc7");

        PrintImg(ParentAcceptances["acc7"],SplineAcceptances["acc7"],Contaminations["acc7"],
                 den,charge,single_acc7["point"],total_fraction,sec_track,l1_cut,fraction_spline_tot,"tot");
    }
    return 1;
}
std::vector<TH1D*> LoadLT(unsigned int charge, TString timePeriod) {
    std::vector<TH1D*> v;
    for (auto k : allowedZ) {
        if (k > charge) {
            //Getting livetime
            TString fileName= "../IonsSelected/Si/Livetime/livetime.root";
            TFile *file = TFile::Open(fileName.Data());
            int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
            auto lvt = (TH1D *)((TH2D *)file->Get("lvt_25"))->ProjectionY("lvt", start_bin, stop_bin);
            v.push_back(lvt);
        }
    }
    return v;
} 
TH1D *LoadRBinWidth() {
    auto hwidth = (TH1D*)hist_rig_highZ->Clone();
	for(int i=1; i<=hist_rig_highZ->GetNbinsX(); ++i) {
		hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
		hwidth->SetBinError(i, 0);
	}
    return hwidth;
}
std::vector<TH1D*> GetParentCounts(unsigned int charge,TString timePeriod) {
    std::vector<TH1D*> v;
    int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
    TString path = "../IonsSelected/"+getIonPath(charge)+"/Fragments/";
    //for the nucleus under study, get the counts from data
    TFile *file = TFile::Open("../IonsSelected/"+getIonPath(charge)+"/Counts/counts.root");
    auto counts = (TH1D *)((TH2D *)file->Get("IL1/rigidity"))->ProjectionY("rr", start_bin, stop_bin);
    //subtract the purity
    TString contName= "../Fragmentation/BelowL1/Fractions/contamination_"+getIonName(charge)+".root";
    TFile *contFile = TFile::Open(contName.Data());
    auto purity = (TF1*)contFile->Get("purityFit");
    auto hcounts = new TH1D("cc", "R (GV)", nRbins_HighZ-1, Rbins_HighZ);
    for (int ibin = 0; ibin < counts->GetNbinsX(); ibin++) {
        hcounts->SetBinContent(ibin + 1, counts->GetBinContent(ibin + 1));
                                           
    }
    v.push_back(hcounts);
    
    for (auto i : allowedZ) {
            if (i<=charge) continue;
            std::cout << "Getting counts for charge : " << i << std::endl;
            auto a = path+McFileName(i)+".root";
            TFile *file = TFile::Open(a.Data(), "READ");
            if (!file||file->IsZombie()) {
                delete file;
                continue;
            }
            TH1D *count = (TH1D*)file->Get("hRig_IL1");
            count->SetTitle("counts_"+McFileName(i));
            count->SetName("counts_"+McFileName(i));
            v.push_back(count);
    }
    return v; 
}
std::map<std::string,std::vector<TH1D*>> GetParentAcceptances(unsigned int charge, std::vector<TString> &Parents, 
                                                              TString sec_track, TString l1_cut) {
    std::map<std::string,std::vector<TH1D*>> v;
    TString path = "../IonsSelected/"+getIonPath(charge)+"/Fragments/";
    for (auto i : allowedZ) {
        if (i <= charge) continue;
        auto a = path+McFileName(i)+".root";
        if (sec_track=="y") a = path+"/WithSecTrack/"+McFileName(i)+".root";
        if (l1_cut=="y" )   a = path+"/WithL1Cut/"+McFileName(i)+".root";
        std::cout << "Using this fragment files : " << a << std::endl;
        TFile *file = TFile::Open(a.Data(), "READ");
        if (!file||file->IsZombie()) {
            delete file;
            continue;
        }
        Parents.push_back(McFileName(i)(0,McFileName(i).First('.')) );

        //Retrieving acceptance acc4 and acc7
        TH1D *mc_pass_acc4 = (TH1D*)file->Get("mc_pass_acc4");
        TH1D *mc_pass_acc7 = (TH1D*)file->Get("mc_pass_acc7");
        TH1D *mc_samp = (TH1D*)file->Get("mc_samp");
        mc_pass_acc4->Sumw2();
        mc_pass_acc7->Sumw2();
        mc_samp->Sumw2();
        mc_pass_acc4->Rebin(9);
        mc_pass_acc7->Rebin(9);
        mc_samp->Rebin(9);

        auto acc_4 = divide(mc_pass_acc4,mc_samp,"efficiency");
        auto acc_7 = divide(mc_pass_acc7,mc_samp,"efficiency");
        acc_4->Scale(TMath::Pi() * 3.9 * 3.9,"nosw2");
        acc_4->SetTitle("acc_"+McFileName(i));
        acc_4->SetName("acc_"+McFileName(i));

        acc_7->Scale(TMath::Pi() * 3.9 * 3.9,"nosw2");
        acc_7->SetTitle("acc_"+McFileName(i));
        acc_7->SetName("acc_"+McFileName(i));
        v["acc4"].push_back(acc_4);
        v["acc7"].push_back(acc_7);
    }
    return v;
}
std::map<std::string,std::vector<TH1D*>> SplineParentAcceptances(std::map<std::string,std::vector<TH1D*>> acc, 
                                           std::vector<TString> Parents, unsigned int charge) {
    std::map<std::string,std::vector<TH1D*>> v;
    //Normalize to the nuclei under study
    /*for (int i=0; i<acc.size(); i++) {
        acc[i]->Divide(acc[0]);
    }*/
    auto acc_knots =SetKnots(SplineUtility::Efficiency::Acc,charge);
    auto acc_range =SetFitLimits(SplineUtility::Efficiency::Acc,charge);
    for (int i=0; i<acc["acc4"].size(); i++) {
        auto h_acc4 = autospline(acc["acc4"][i],acc_range.first,acc_range.second,acc_knots.first,acc_knots.second);
        auto h_acc7 = autospline(acc["acc7"][i],acc_range.first,acc_range.second,acc_knots.first,acc_knots.second);
        h_acc4->SetTitle("spl_"+Parents[i]);
        h_acc4->SetName("spl_"+Parents[i]);
        h_acc7->SetTitle("spl_"+Parents[i]);
        h_acc7->SetName("spl_"+Parents[i]);
        v["acc4"].push_back(h_acc4);
        v["acc7"].push_back(h_acc7);
    }
    return v;
}
std::vector<TH1D*> GetGalpropFluxes(unsigned int charge) {
    std::vector<TH1D*> v;
    TString path = "../../IonsAnalysis/Published/Galprop/";
    for (auto k : allowedZ) {
        TString a = path+getIonPath(k)+".txt";
        std::ifstream infile(a.Data());
        if (!infile) {
            std::cerr << "Error: Could not open "<< a << std::endl;
            continue;
        }
        std::vector<double> rigidity, intensity;
        double r, i;
        while (infile >> r >> i) {
            rigidity.push_back(r);
            intensity.push_back(i);
        }
        int nPoints = rigidity.size();
        if (nPoints == 0) {
            std::cerr << "Error: No data read from the file." << std::endl;
        }
        TGraph *graph = new TGraph(nPoints, rigidity.data(), intensity.data());
        auto h = new TH1D ("", "", nRbins_HighZ-1, Rbins_HighZ);
        for (int i = 0; i < h->GetNbinsX(); i++) {
            double binCenter = h->GetBinCenter(i);
            double value = graph->Eval(binCenter);
            h->SetBinContent(i, value);
        }
        //MultiplyByXPower(h,2.7);
        h->SetTitle("flux_"+McFileName(i));
        h->SetName("flux_"+McFileName(i));
        v.push_back(h);
    }
    return v;
}
std::vector<TH1D*> GetQyanFluxes(unsigned int charge) {
    std::vector<TH1D*> v;
    for (auto k : allowedZ) {
        if (k > charge) {
            if (k == 21 || k == 22 || k == 23 || k == 24 || k == 25) {
                // per Sc, Ti, V usa il flusso dummy
                v.push_back(buildDummyFlux(k,"13.5y",""));
            } else {
                auto hh = GetQYanFlux(k);
                MultiplyByXPower(hh,-2.7);
                auto hhh = RebinHistogramWeighted(hh);
                v.push_back(hhh);
            }
        }
    }
    return v;
}
std::map<std::string,std::vector<TH1D*>> Contamination(std::map<std::string,std::vector<TH1D*>> spline_acc,std::vector<TH1D*> fluxes, 
                                std::vector<TH1D*> lvt, TH1D* rwidth,std::vector<TString> Parents) {
    TH1::SetDefaultSumw2();
    std::map<std::string,std::vector<TH1D*>> v;
    for (int i=0; i<spline_acc["acc4"].size(); i++) {

        spline_acc["acc4"][i]->Sumw2();
        TH1D *h_acc4 = (TH1D*)spline_acc["acc4"][i]->Clone( Form("h_acc4_%d",i) );
        h_acc4->Multiply(lvt[i]);
        h_acc4->Multiply(rwidth);
        h_acc4->Multiply(fluxes[i]);
        h_acc4->SetName(Parents[i]);
        h_acc4 ->SetTitle(Parents[i]);

        spline_acc["acc7"][i]->Sumw2();
        TH1D *h_acc7 = (TH1D*)spline_acc["acc7"][i]->Clone( Form("h_acc7_%d",i) );
        h_acc7->Multiply(lvt[i]);
        h_acc7->Multiply(rwidth);
        h_acc7->Multiply(fluxes[i]);
        h_acc7->SetName(Parents[i]);
        h_acc7 ->SetTitle(Parents[i]);

        v["acc4"].push_back(h_acc4);
        v["acc7"].push_back(h_acc7);
    }
    return v;
}
std::map<std::string,std::vector<TH1D*>> SingleContribution(std::vector<TH1D*> counts, TH1D* den) {
    std::map<std::string,std::vector<TH1D*>> v;
    for (int i=0; i<counts.size(); i++) {
            TH1D* hh = divide(counts[i],den,"ratio");
            hh->SetName( Form("single_%s",counts[i]->GetName()));
            hh->SetTitle( Form("single_%s",counts[i]->GetName()));
            auto hh_spline = autospline(hh,0.8,1000,2,10);
            hh_spline->SetName( Form("single_spline_%s",counts[i]->GetName()));
            hh_spline->SetTitle( Form("single_spline_%s",counts[i]->GetName()));
            v["spline"].push_back(hh_spline);
            v["point"].push_back(hh);
    }
    return v;
}
std::map<std::string,TH1D*> GetSum(std::map<std::string,std::vector<TH1D*>>  v, int start) {
    std::map<std::string,TH1D*> g;
    TH1D *h_acc4 = new TH1D("","",nRbins_HighZ-1,Rbins_HighZ);
    TH1D *h_acc7 = new TH1D("","",nRbins_HighZ-1,Rbins_HighZ);
    for (int i=0; i<v["acc4"].size(); i++) {
        if (i>=start) {
            h_acc4->Add(v["acc4"][i]);
            h_acc7->Add(v["acc7"][i]);
        }
    } 
    g["acc4"] = h_acc4;
    g["acc7"] = h_acc7;
    return g;
}
void PrintImg(std::vector<TH1D*> acc, std::vector<TH1D*> spl, std::vector<TH1D*> cont, TH1D *sum, unsigned int charge,
              std::vector<TH1D*> SingleContribution, TH1D* fraction, TString sec_track, TString l1_cut, TH1D* fraction_spline, TString ACC) {

    // --- safety debug ---
    std::cerr << "PrintImg debug: acc=" << acc.size()
              << " spl=" << spl.size()
              << " cont=" << cont.size()
              << " SingleContribution=" << SingleContribution.size()
              << "\n";
    if (!sum) std::cerr << "Warning: sum is nullptr\n";
    if (!fraction) std::cerr << "Warning: fraction is nullptr\n";
    if (!fraction_spline) std::cerr << "Warning: fraction_spline is nullptr\n";

    // style + canvas
    TStyle *style = effHistStyle(0);
    gROOT->SetStyle("MyStyle");
    gROOT->ForceStyle();
    setLogon();
    TCanvas *b = new TCanvas("","",2048,1280);
    b->SetLogy();
    b->SetLogx();

    // open pdf (attenzione all'uso corretto delle parentesi [ e ])
    TString pdfBase = Form("../IonsSelected/%s/Fragments/Result/img_%s.pdf", getIonPath(charge).Data(), ACC.Data());
    if (sec_track != "y") b->SaveAs(pdfBase + "[");
    else b->SaveAs(Form("../IonsSelected/%s/Fragments/Result/img_sectrack_%s.pdf[", getIonPath(charge).Data(), ACC.Data()));

    // --- acceptances ---
    auto a = new TLegend(0.60,0.9,0.9,0.65);
    a->SetNColumns(2);
    TH2D *h = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.5*1e-8, 0.1);
    setTitle(h, "", "R (GV)", "Acceptance (m^{2} sr)", charge);
    formatAxis(h, 1);
    formatTitle(h,1);
    formatZoomY(h, "Spline mc accepatnce");
    h->Draw();

    int col_index = 0;
    for (auto *iph : acc) {
        if (!iph) { std::cerr << "Warning: acc element nullptr\n"; continue; }
        TString t = iph->GetTitle();
        a->AddEntry((TObject*)iph, t(t.Index("_")+1,t.Length()).Data());
        iph->SetMarkerSize(2.2);
        iph->SetMarkerStyle(20);
        iph->SetLineWidth(1.8);
        iph->GetXaxis()->SetRangeUser(2.15,999);
        setFriendColor(const_cast<TH1D*>(iph), col_index); // se setFriendColor prende TH1D*
        iph->Draw("hist PE SAME");
        ++col_index;
    }

    col_index = 0;
    for (auto *sph : spl) {
        if (!sph) { std::cerr << "Warning: spl element nullptr\n"; continue; }
        sph->SetMarkerSize(2.2);
        sph->SetMarkerStyle(20);
        sph->SetLineWidth(1.8);
        sph->GetXaxis()->SetRangeUser(2.15,1000);
        setFriendColor(const_cast<TH1D*>(sph), col_index);
        sph->DrawCopy("hist L SAME");
        ++col_index;
    }

    a->SetBorderSize(0);
    a->SetFillStyle(0);
    a->SetTextSize(0.028);
    a->Draw("SAME");

    if (sec_track != "y") b->SaveAs(pdfBase);
    else b->SaveAs(Form("../IonsSelected/%s/Fragments/Result/img_sectrack_%s.pdf", getIonPath(charge).Data(), ACC.Data()));

    if (!sum) { std::cerr << "Error: sum histogram missing, skipping sum plot\n"; }
    else {
        // All corrections + sum
        TH2D *h1 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.001, 1e6);
        setTitle(h1, "", "R (GV)", "Counts", charge);
        formatAxis(h1, 1);
        formatTitle(h1,1);
        formatZoomY(h1, "Spline mc accepatnce");
        h1->Draw("");
        // draw cont
        auto c = new TLegend(0.60,0.9,0.9,0.65);
        c->SetBorderSize(0); c->SetFillStyle(0); c->SetNColumns(2);
        col_index = 0;
        for (auto *ct : cont) {
            if (!ct) { std::cerr << "Warning: cont element nullptr\n"; ++col_index; continue; }
            TString t = ct->GetTitle();
            c->AddEntry((TObject*)ct, t(t.Index("_")+1,t.Length()).Data());
            ct->SetMarkerSize(2.2);
            ct->SetMarkerStyle(20);
            setFriendColor(const_cast<TH1D*>(ct), col_index);
            ct->Draw("hist P SAME");
            ++col_index;
        }
        c->SetTextSize(0.028);
        c->Draw("SAME");

        // draw sum
        sum->SetMarkerSize(2);
        sum->SetMarkerStyle(20);
        sum->Draw("hist P SAME");
    }

    if (sec_track != "y") b->SaveAs(pdfBase);
    else b->SaveAs(Form("../IonsSelected/%s/Fragments/Result/img_sectrack_%s.pdf", getIonPath(charge).Data(), ACC.Data()));

    // --- Single contributions and fractions ---
    double ymax = 0.1; // fallback
    switch(charge) {
        case 14: ymax = 0.014; break;
        case 15: ymax = 0.15;  break;
        case 16: ymax = 0.07;  break;
        case 18: ymax = 0.09;  break;
        case 20: ymax = 0.09;  break;
        default: break;
    }

    TH2D *h2 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0., ymax);
    setTitle(h2, "", "R (GV)", "Fractions", charge);
    h2->SetTitle( Form("%s_%s",h2->GetTitle(),ACC.Data() ));
    formatAxis(h2, 1);
    formatTitle(h2,1);
    formatZoomY(h2, "Spline mc accepatnce");
    getColor(h2,charge,"dat");
    h2->GetXaxis()->SetRangeUser(2.15,1000);
    h2->Draw("");
    auto c1 = new TLegend(0.70,0.9,0.9,0.55);
    c1->SetBorderSize(0); c1->SetFillStyle(0); c1->SetNColumns(2);
    b->SetLogy(0);

    if (fraction) {
        fraction->GetYaxis()->SetRangeUser(1.,ymax);
        fraction->GetXaxis()->SetRangeUser(2.15,1000);
        fraction->Draw("SAME");
        fraction->SetMarkerSize(2.2);
        fraction->SetMarkerStyle(20);
        fraction->SetMarkerColor(kBlack);
        fraction->SetLineColor(kBlack);
        c1->AddEntry(fraction,"Total");
    }

    col_index = 0;
    for (size_t i = 0; i < SingleContribution.size(); ++i) {
        TH1D *sc = SingleContribution[i];
        if (!sc) { std::cerr << "Warning: SingleContribution["<<i<<"] is nullptr\n"; ++col_index; continue; }
        TString t = sc->GetTitle();
        c1->AddEntry(sc, t(t.Index("_")+1,t.Length()).Data());
        sc->SetMarkerSize(2.2);
        sc->SetMarkerStyle(20);
        sc->SetLineWidth(1.8);
        sc->GetXaxis()->SetRangeUser(2.15,1000);
        setFriendColor(sc, col_index);
        sc->DrawCopy("SAME");
        ++col_index;
    }

    if (fraction_spline) {
        getColor(fraction_spline,999,"mc");
        fraction_spline->SetMarkerSize(2.2);
        fraction_spline->SetMarkerStyle(20);
        fraction_spline->SetLineWidth(1.8);
        fraction_spline->GetXaxis()->SetRangeUser(2.15,1000);
        fraction_spline->SetFillStyle(0);
        fraction_spline->DrawCopy("hist L SAME");
        fraction_spline->SetFillStyle(1001);
        fraction_spline->SetMarkerSize(0);
        fraction_spline->Draw("E3 L SAME");
    }

    c1->SetTextSize(0.038);
    c1->Draw("SAME");

    // close pdf (use the closing bracket)
    if (sec_track != "y") b->SaveAs(pdfBase);
    else b->SaveAs(Form("../IonsSelected/%s/Fragments/Result/img_sectrack_%s.pdf", getIonPath(charge).Data(), ACC.Data()));

    if (sec_track != "y") b->SaveAs(pdfBase + "]");
    else b->SaveAs(Form("../IonsSelected/%s/Fragments/Result/img_sectrack_%s.pdf]", getIonPath(charge).Data(), ACC.Data()));

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
        case 21:
            a = "Sc.B1236";
            break;
        case 22:
            a = "Ti.B1236";
            break;
        case 23:
            a = "V.B1236";
            break;
        case 24:
            a = "Cr.B1236";
            break;
        case 25:
            a = "Mn.B1236";
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
TH1D *buildDummyFlux(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting counts
    TString fileName= Form("../IonsSelected/"+getIonPath(charge)+"/Counts/%scounts.root",inp_sec_track.Data());
    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
        printf("Errore nell'aprire il file\n");
    }
    int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
    auto counts = (TH1D *)((TH2D *)file->Get("IL1/rigidity"))->ProjectionY("r", start_bin, stop_bin);
    if (!counts) printf("Errore \n");
    //Getting livetime
    TString fileName1= "../IonsSelected/"+getIonPath(charge)+"/Livetime/livetime.root";
    TFile *file1 = TFile::Open(fileName1.Data());
    TH1D *lvt = (TH1D *)((TH2D *)file1->Get("lvt_25"))->ProjectionY("q", start_bin, stop_bin);
    //Getting Final Raw acceptance
    TString fileName2= Form("../IonsSelected/"+getIonPath(charge)+"/RawAcceptance/rawacc.root");
    TFile *file2 = TFile::Open(fileName2.Data());
    auto mc_pass_gen = (TH1D*)file2->Get("mc_pass_gen");
    auto mc_samp = (TH1D*)file2->Get("mc_samp");
    mc_pass_gen->Rebin(9);
    mc_samp->Rebin(9);
    auto acc = divide(mc_pass_gen, mc_samp, "efficiency");
    acc->Scale(TMath::Pi() * 3.9 * 3.9);
    auto acc_range= SetFitLimits(SplineUtility::Efficiency::Acc,charge);
    auto acc_knots =SetKnots(SplineUtility::Efficiency::Acc,charge);
    auto spline_acc = autospline(acc, acc_range.first, acc_range.second,acc_knots.first,acc_knots.second);
    //Bin Width
    auto hwidth = (TH1D*)hist_rig_highZ->Clone();
    for(int i=1; i<=hist_rig_highZ->GetNbinsX(); ++i) {
      hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
      hwidth->SetBinError(i, 0);
    }
    //Raw Flux
    for (int ibin = 1; ibin <= counts->GetNbinsX(); ibin++) {
        double N = counts->GetBinContent(ibin);
        double LT = lvt->GetBinContent(ibin);
        double A = spline_acc->GetBinContent(ibin);
        double dR = hwidth->GetBinContent(ibin);
        double sigma_N = counts->GetBinError(ibin);
        double sigma_LT = lvt->GetBinError(ibin);
        double sigma_A = spline_acc->GetBinError(ibin);
        double sigma_dR = hwidth->GetBinError(ibin);
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
std::map<std::string,TH1D*> WeightForToi(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting livetime
    TString fileName= "../IonsSelected/"+getIonPath(charge)+"/Livetime/livetime.root";
    TFile *file = TFile::Open(fileName.Data());

    //produce the amount of livetime we were in ACC4 and ACC7
    int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
    TH1D* lvt_acc4 = (TH1D *)((TH2D *)file->Get("lvt_25"))->ProjectionY("lvt_acc4", hist_time->FindBin(0),hist_time->FindBin(1456503197));
    TH1D* lvt_acc7 = (TH1D *)((TH2D *)file->Get("lvt_25"))->ProjectionY("lvt_acc7", hist_time->FindBin(1456503197), stop_bin);

    //produce the two weights
    TH1D* w_acc4 = (TH1D*)hist_rig_highZ->Clone("w_acc4");
    TH1D* w_acc7 = (TH1D*)hist_rig_highZ->Clone("w_acc7");
    for (int i = 1; i <= hist_rig_highZ->GetNbinsX(); ++i) {
        double t4  = lvt_acc4->GetBinContent(i);
        double t7  = lvt_acc7->GetBinContent(i);
        double e4  = lvt_acc4->GetBinError(i);
        double e7  = lvt_acc7->GetBinError(i);
        double tot = t4 + t7;

        if (tot <= 0) {
            w_acc4->SetBinContent(i, 0);
            w_acc4->SetBinError(i, 0);
            w_acc7->SetBinContent(i, 0);
            w_acc7->SetBinError(i, 0);
            continue;
        }

        double w4 = t4 / tot;
        double w7 = t7 / tot;

        double err_w4 = sqrt( (t7*t7*e4*e4 + t4*t4*e7*e7) ) / pow(tot, 2);
        double err_w7 = err_w4; // perché w7 = 1 - w4

        w_acc4->SetBinContent(i, w4);
        w_acc4->SetBinError(i, err_w4);
        w_acc7->SetBinContent(i, w7);
        w_acc7->SetBinError(i, err_w7);
    }
    std::map<std::string,TH1D*> g;
    g["acc4"] = w_acc4;
    g["acc7"] = w_acc7;
    return g;   
}
TH1D* TotalWeightedToi(std::map<std::string,TH1D*> w, TH1D* f_acc4, TH1D* f_acc7) {
    TH1D* tot = (TH1D*)hist_rig_highZ->Clone("tot");

    for (int i = 1; i <= hist_rig_highZ->GetNbinsX(); ++i) {
        double w4 = w["acc4"]->GetBinContent(i);
        double w7 = w["acc7"]->GetBinContent(i);
        double ew4 = w["acc4"]->GetBinError(i);
        double ew7 = w["acc7"]->GetBinError(i);

        double toi4 = f_acc4->GetBinContent(i);
        double toi7 = f_acc7->GetBinContent(i);
        double etoi4 = f_acc4->GetBinError(i);
        double etoi7 = f_acc7->GetBinError(i);

        double val = w4 * toi4 + w7 * toi7;

        // Propagazione completa dell’errore
        double err2 =
            pow(toi4 * ew4, 2) +
            pow(toi7 * ew7, 2) +
            pow(w4 * etoi4, 2) +
            pow(w7 * etoi7, 2);

        double err = sqrt(err2);

        tot->SetBinContent(i, val);
        tot->SetBinError(i, err);
    }

    return tot;
}

