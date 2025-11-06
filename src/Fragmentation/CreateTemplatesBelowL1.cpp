#include "definition.h"
#include "binning.h"
#include "utils.h"

void sum(int charge);
std::pair< std::vector<TH1D*>,std::vector<TH1D*> > getTemplates(unsigned int charge);

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <charge>" << std::endl;
        return 1;
    }
    int charge = std::stoi(argv[1]);
    sum(charge);
    return 0;
}

void sum(int charge) {
    std::pair< std::vector<TH1D*>,std::vector<TH1D*> > z1 = getTemplates(charge);

    //Crea il file .root di output
    TString nucleus = getIonName(charge);
    TString path = "../Fragmentation/BelowL1/Templates/";
    TString froot = path+nucleus+".root";
    auto outfile = new TFile(froot, "RECREATE");

    // Crea la canvas per gli istogrammi sommati e bin x bin
    TCanvas *canvasSum = new TCanvas("canvasSum", "Summed Histograms", 800, 600);
    canvasSum->SetLogy();
    // Crea un file PDF per l'output
    TString out1 = path+nucleus+"_Bin.pdf[";
    TString out2 = path+nucleus+"_Bin.pdf";
    TString out3 = path+nucleus+"_Bin.pdf]";
    TCanvas *canvasBin = new TCanvas("canvasBin", "Bin", 800, 600);
    canvasBin->SetLogy();
    canvasBin->cd();
    canvasBin->Print(out1.Data() );

    // Somma tutti gli istogrammi in L1
    TH1D *sumL1 = nullptr;
    TH1D *sumL2 = nullptr;
    TString hist1_name;
    TString hist2_name;
    int binL2z1 = z1.second.size();
    for (int i = 0; i < binL2z1; i++) {
        // Costruisci i nomi degli istogrammi
        /*if (i <= 9) {
            hist1_name = "L1Template_00" + std::to_string(i);
            hist2_name = "L2Template_00" + std::to_string(i);
        }
        else {
            hist1_name = "L1Template_0" + std::to_string(i);
            hist2_name = "L2Template_0" + std::to_string(i);
        }
        // Recupera gli istogrammi dalle liste
        TH1D *hist1 = (TH1D*)list1->FindObject(hist1_name.Data());
        TH1D *hist2 = (TH1D*)list2->FindObject(hist2_name.Data());
        if (!hist1 || !hist2) {
            printf("Errore nel recuperare gli istogrammi %s o %s\n", hist1_name.Data(), hist2_name.Data());
            continue;
        }*/
         // Setta i marker style e colore
        z1.first.at(i)->SetMarkerStyle(20);
        z1.first.at(i)->SetMarkerColor(kRed);
        z1.first.at(i)->SetLineColor(kRed);
        z1.second.at(i)->SetMarkerStyle(20);
        z1.second.at(i)->SetMarkerColor(kBlue);
        z1.second.at(i)->SetLineColor(kBlue);
        // Imposta il range dell'asse X
        std::pair<int,int> range = getRange(charge);
        z1.first.at(i)->GetXaxis()->SetRangeUser(range.first, range.second);
        z1.second.at(i)->GetXaxis()->SetRangeUser(range.first, range.second);
        // Disegna gli istogrammi sommati
        z1.first.at(i)->SetTitle(nucleus.Data() );
        z1.first.at(i)->GetXaxis()->SetTitle("Q");
        z1.first.at(i)->GetYaxis()->SetTitle("Normalized Entries");
        z1.second.at(i)->SetTitle(nucleus.Data() );
        z1.second.at(i)->GetXaxis()->SetTitle("Q");
        z1.second.at(i)->GetYaxis()->SetTitle("Normalized Entries");
        z1.first.at(i)->SetStats(kFALSE);
        z1.second.at(i)->SetStats(kFALSE);
        //Normalizzali
        //z1.first.at(i)->Scale(1.0/z1.first.at(i)->GetEntries() );
        //z1.second.at(i)->Scale(1.0/z1.second.at(i)->GetEntries() );
        // Aggiungi una legenda
        TLegend *legend1 = new TLegend(0.7, 0.8, 0.9, 0.9);
        legend1->AddEntry(z1.first.at(i), "L1", "lp");
        legend1->AddEntry(z1.second.at(i), "L2", "lp");
        z1.first.at(i)->Draw();
        z1.second.at(i)->Draw("SAME");
        legend1->Draw();
        // Aggiungi la pagina al PDF Bin
        canvasBin->Print(out2.Data() );
        delete legend1;

        // Somma gli istogrammi
        if (sumL1 == nullptr) {
            sumL1 = (TH1D*)z1.first.at(i)->Clone("sumL1");
        } else {
            sumL1->Add(z1.first.at(i));
        }

        if (sumL2 == nullptr) {
            sumL2 = (TH1D*)z1.second.at(i)->Clone("sumL2");
        } else {
            sumL2->Add(z1.second.at(i));
        }
    }

    if (!sumL1 || !sumL2) {
        printf("Errore: le somme degli istogrammi non sono state create correttamente.\n");
        return;
    }

    // Setta i marker style e colore
    sumL1->SetMarkerStyle(20);
    sumL1->SetMarkerColor(kRed);
    sumL1->SetLineColor(kRed);
    sumL2->SetMarkerStyle(20);
    sumL2->SetMarkerColor(kBlue);
    sumL2->SetLineColor(kBlue);
    // Imposta il range dell'asse X
    std::pair<int,int> range = getRange(charge);
    sumL1->GetXaxis()->SetRangeUser(range.first, range.second);
    sumL2->GetXaxis()->SetRangeUser(range.first, range.second);
    // Disegna gli istogrammi sommati
    sumL1->SetTitle(nucleus.Data() );
    sumL1->GetXaxis()->SetTitle("Q");
    sumL1->GetYaxis()->SetTitle("Entries");
    sumL2->SetTitle(nucleus.Data() );
    sumL2->GetXaxis()->SetTitle("Q");
    sumL2->GetYaxis()->SetTitle("Entries");
    sumL1->SetStats(kFALSE);
    sumL2->SetStats(kFALSE);
    //Normalizzali
    //sumL1->Scale(1.0/sumL1->GetEntries() );
    //sumL2->Scale(1.0/sumL2->GetEntries() );
    canvasSum->cd();
    sumL1->Draw("E");
    sumL2->Draw("E SAME");

    // Aggiungi una legenda
    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(sumL1, "L1", "lp");
    legend->AddEntry(sumL2, "L2", "lp");
    legend->Draw();

    // Aggiungi la pagina al PDF
    TString out4 = path+nucleus+"_Summed.pdf";
    TString out5 = path+nucleus+"_Summed.pdf]";
    canvasSum->Print(out4.Data() );

    // Chiudi i pdf
    canvasSum->Print(out5.Data() );
    canvasBin->cd();
    canvasBin->Print(out3.Data() );

    //Salva il file .root
    outfile->WriteTObject(sumL1, "sumL1");
	outfile->WriteTObject(sumL2, "sumL2");
    outfile->Close();

    // Pulizia
    delete legend;
    delete canvasBin;
    delete canvasSum;

}

std::pair< std::vector<TH1D*>,std::vector<TH1D*> > getTemplates(unsigned int charge) {
    std::pair< std::vector<TH1D*>,std::vector<TH1D*> > p {};
    // Costruisci il nome della lista basato sulla carica
    TString listName1 = "L1TemplateList_" + TString::Format("%d", charge);
    TString listName2 = "L2TemplateList_" + TString::Format("%d", charge);
    // Apri il file ROOT
    TString path1 = "../Fragmentation/BelowL1/";
    TString nucleus = getIonName(charge);
    TString fileName = path1+nucleus+"/"+nucleus+".root";
    TFile *file = TFile::Open(fileName.Data() );
    if (!file || file->IsZombie()) printf("Errore nell'aprire il file\n");
    // Recupera le liste
    TList *list1 = (TList*)file->Get(listName1.Data() );
    TList *list2 = (TList*)file->Get(listName2.Data() );
    if (!list1 || !list2) printf("Errore nel recuperare le liste\n");
    TString hist1_name;
    TString hist2_name;
    // Get L1 Template
    for (auto *obj : *list1) {
      TH1D *hist1 = dynamic_cast<TH1D*>(obj);
        if (hist1) {
          p.first.push_back(hist1);
        } else {
          printf("Errore nel recuperare un istogramma da L1TemplateList per carica %u \n", charge);
        }
    }
    // Get L2 Template
    for (auto *obj : *list2) {
        TH1D *hist2 = dynamic_cast<TH1D*>(obj);
        if (hist2) {
            p.second.push_back(hist2);
        } else {
            printf("Errore nel recuperare un istogramma da L2TemplateList per carica %u \n", charge);
        }
    }
    return p;
}
