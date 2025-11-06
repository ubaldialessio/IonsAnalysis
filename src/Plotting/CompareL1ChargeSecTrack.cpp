#include "binning.h"
#include "utils.h"

std::map<std::string, std::vector<TH1D*>> getTemplates(unsigned int charge, TString cut);
std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> loadTemplates(unsigned int charge, TString cut);
void print(unsigned int charge,
           std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> sec_track_templates, 
           std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> templates,
           TString cut);

int main(int argc, char *argv[]) {
    if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> -cut\n",argv[0]);
        printf("-cut can be: sec_track/l1_cut\n");
	    return 1;
	}
    unsigned int charge=atoi(argv[1]);
    TString cut = argv[2];
    auto templates_sec_track = loadTemplates(charge,cut);
    auto templates_no_sec_track = loadTemplates(charge,"");
    print(charge, templates_sec_track,templates_no_sec_track,cut);
}

std::map<std::string, std::vector<TH1D*>> getTemplates(unsigned int charge, TString cut) {
    std::map<std::string, std::vector<TH1D*>> templates;
    // Costruisci il nome della lista basato sulla carica
    TString listName1 = "L1TemplateList_" + TString::Format("%d", charge);
    TString listName2 = "L2TemplateList_" + TString::Format("%d", charge);
    // Apri il file ROOT
    TString path1 = "../Fragmentation/BelowL1/";
    TString nucleus = getIonName(charge);
    TString fileName = Form("%s%s/%s%s.root",path1.Data(),nucleus.Data(),cut.Data(),nucleus.Data());
    std::cout << fileName << std::endl;
    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
        printf("Errore nell'aprire il file\n");
        return templates;
    }
    // Recupera le liste
    TList *list1 = dynamic_cast<TList*>(file->Get(listName1.Data()));
    TList *list2 = dynamic_cast<TList*>(file->Get(listName2.Data()));
    if (!list1 || !list2) {
        printf("Errore nel recuperare le liste\n");
        return templates;
    }
    // Get L1 Templates
    for (auto *obj : *list1) {
        TH1D *hist1 = dynamic_cast<TH1D*>(obj);
        if (hist1) {
            templates["L1"].push_back(hist1);
        } else {
            printf("Errore nel recuperare un istogramma da L1TemplateList per carica %u\n", charge);
        }
    }
    // Get L2 Templates
    for (auto *obj : *list2) {
        TH1D *hist2 = dynamic_cast<TH1D*>(obj);
        if (hist2) {
            templates["L2"].push_back(hist2);
        } else {
            printf("Errore nel recuperare un istogramma da L2TemplateList per carica %u\n", charge);
        }
    }
    return templates;
}
std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> loadTemplates(unsigned int charge, TString cut) {
    std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> templatesByCharge;
    //Original with charge, charge+1, charge+2
    /*for (unsigned int i = 0; i <= n; i++) {
        unsigned int currentCharge = charge + i;
        templatesByCharge[currentCharge] = getTemplates(currentCharge); 
    }*/
    //Now including also charge -1 
    std::vector<unsigned int> allowedZ = {charge};
    for (auto k : allowedZ)
        templatesByCharge[k] = getTemplates(k,cut); 
    return templatesByCharge;
}
void print(unsigned int charge,
           std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> cut_templates, 
           std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> templates,
           TString cut) {
    TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();

    const int nP_rigBins = 13; 
	const double P_rigBins[nP_rigBins] = {0.8,3.64,5.37,7.76,11.00,15.3,21.1,38.9,70.00,150.00,300.00,500.00,1000.00};	

    auto templ        = templates[charge]["L1"];
    auto secTrackTemplate= cut_templates[charge]["L1"];
    auto nBins = templ.size();
    auto ion = getIonName(charge).Data();
    TString output = Form("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/%s_templateComparison_%s.pdf",ion,cut.Data());
    TCanvas* c = new TCanvas("c", "", 1024,640);
    c->SetLogy();
    for (int i=0; i < nBins ; ++i) {
        if (i==0) c->SaveAs( Form("%s[",output.Data()) );

        TH1D* temp      =templ.at(i);
        TH1D* temp_track=secTrackTemplate.at(i);
        formatAxis(temp,1);
        setTitle(temp, "", "Q_{L1}", "Counts", charge);
        formatTitle(temp,1);
        formatMarkerSize(temp,1);
        formatMarkerSize(temp_track,1);
        temp_track->SetMarkerColor(kRed);
        temp->SetMarkerColor(kBlack);
        temp_track->SetMarkerStyle(20);
        temp->SetMarkerStyle(20);
        temp->GetXaxis()->SetRangeUser(charge-20,charge+20);
        temp_track->GetXaxis()->SetRangeUser(charge-12,charge+12);

        auto leg = new TLegend(0.73, 0.75, 0.90, 0.9);
        leg->AddEntry(temp,"Default");
        leg->AddEntry(temp_track, Form("%s",cut.Data()) );
        leg->AddEntry((TObject*)nullptr, Form("[%.1f,%.1f] GV",P_rigBins[i],P_rigBins[i+1]),"");
        leg->SetTextSize(0.033);
        temp->Draw();
        temp_track->Draw("SAME");
        leg->Draw("SAME");
        c->SaveAs( Form("%s",output.Data()) );

        if (i== nBins-1) c->SaveAs( Form("%s]",output.Data()) );
    }

    //ratio 
    TString output2 = Form("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/%s_ratioTemplateComparison_%s.pdf",ion,cut.Data());
    TCanvas* c1 = new TCanvas("c", "", 1024,640);
    c1->SetGridy();
    for (int i=0; i < nBins ; ++i) {
        if (i==0) c1->SaveAs( Form("%s[",output2.Data()) );

        TH1D* temp      =templ.at(i);
        TH1D* temp_track=secTrackTemplate.at(i);
        temp_track->Divide(temp);
        auto ratio = (TH1D*)temp_track->Clone("ratio");
        formatAxis(ratio,1);
        setTitle(ratio, "", "Q_{L1}", "Ratio", charge);
        formatTitle(ratio,1);
        formatMarkerSize(ratio,1);
        ratio->SetMarkerColor(kBlack);
        ratio->SetMarkerStyle(20);
        ratio->GetXaxis()->SetRangeUser(charge-20,charge+20);
        ratio->GetYaxis()->SetRangeUser(0.,1.2);

        auto leg = new TLegend(0.53, 0.78, 0.90, 0.9);
        leg->AddEntry(ratio,Form("ratio %s/default",cut.Data()) );
        leg->AddEntry((TObject*)nullptr, Form("[%.1f,%.1f] GV",P_rigBins[i],P_rigBins[i+1]),"");
        leg->SetTextSize(0.033);
        ratio->Draw("hist P");
        leg->Draw("SAME");
        c1->SaveAs( Form("%s",output2.Data()) );

        if (i== nBins-1) c1->SaveAs( Form("%s]",output2.Data()) );
    }
}