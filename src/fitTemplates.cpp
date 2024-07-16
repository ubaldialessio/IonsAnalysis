#include "definition.h"
#include "binning.h"
#include "utils.h"

//Roofit headers
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooHistConstraint.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooCmdArg.h"
#include "RooArgSet.h"

//additional header
#include <string_view>


using namespace RooFit;

const int nRigbins_LightIons = 24;
const double Rigbins_LightIons[nRigbins_LightIons] = {0.8,  1.16, 1.51, 1.92, 2.40,  2.97, 3.64, 4.43,
	                                                  5.37, 6.47, 7.76, 9.26, 11.00, 13.0, 15.3, 18.0,
	                                                  21.1, 24.7, 28.8, 33.5, 38.9,  45.1, 52.2, 5000};

void fitHistograms(std::pair< std::vector<TH1D*>,std::vector<TH1D*> > z1,
    std::pair< std::vector<TH1D*>,std::vector<TH1D*> > z2, unsigned int charge);

std::pair< std::vector<TH1D*>,std::vector<TH1D*> > getTemplates(unsigned int charge);
std::pair<double,double> getCoefBoron(int i );
std::pair<double,double> getCoefCarbon(int i );

TString getIonName(unsigned int charge);
std::pair<int,int> getRange(unsigned int charge );


int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: \n");
	    printf("%s <charge> \n", argv[0]);
	    return 1;
	}
    unsigned int charge=atoi(argv[1]);
    std::pair< std::vector<TH1D*>,std::vector<TH1D*> > z1 = getTemplates(charge);
    std::pair< std::vector<TH1D*>,std::vector<TH1D*> > z2 = getTemplates(charge+1);
    fitHistograms(z1,z2,charge);
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
    if (!file || file->IsZombie()) {
        printf("Errore nell'aprire il file\n");
    }
    // Recupera le liste
    TList *list1 = (TList*)file->Get(listName1.Data() );
    TList *list2 = (TList*)file->Get(listName2.Data() );
    if (!list1 || !list2) {
        printf("Errore nel recuperare le liste\n");
    }
    TString hist1_name;
    TString hist2_name;
    for (int i = 0; i < 23; i++) {
        // Costruisci i nomi degli istogrammi
        if (i <= 9) {
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
        }
        p.first.push_back(hist1);
        p.second.push_back(hist2);
    }
    return p;
}

void fitHistograms(std::pair< std::vector<TH1D*>,std::vector<TH1D*> > z1,
    std::pair< std::vector<TH1D*>,std::vector<TH1D*> > z2, unsigned int charge) {
    // Crea un file PDF per l'output
    TString path1 = "../Fragmentation/BelowL1/Fractions/";
    TString nucleus = getIonName(charge);
    TString out1 = path1+nucleus+"_purityFit.pdf[";
    TString out2 = path1+nucleus+"_purityFit.pdf";
    TString out3 = path1+nucleus+"_purityFit.pdf]";
    TCanvas *canvasBin = new TCanvas("canvasBin", "Bin", 800, 600);
    canvasBin->SetLogy();
    canvasBin->cd();
    canvasBin->Print(out1.Data() );

    auto f_vs_r = new TGraph();
    f_vs_r->GetXaxis()->SetTitle("R (GV)");
    f_vs_r->GetYaxis()->SetTitle("Fraction");

    for (int i=0; i<23; i++) {
        // Definisci la variabile osservabile
        RooRealVar x("x", "Q", 3, 10);
        // Crea RooDataHist da TH1D ---- first=L1, second=L2
        TH1D* hist1 = z1.first[i];
        TH1D* hist2 = z1.second[i];
        TH1D* hist3 = z2.second[i];

        RooDataHist dataHistSumL1_Boron("dataHistSumL1_Boron", "dataHistSumL1_Boron", x, Import(*hist1) ); //z1.first[i]
        RooDataHist dataHistSumL2_Boron("dataHistSumL2_Boron", "dataHistSumL2_Boron", x, Import(*hist2) ); //z1.second[i]
        RooDataHist dataHistSumL2_Carbon("dataHistSumL2_Carbon", "dataHistSumL2_Carbon", x, Import(*hist3) ); //z2.second[i]
        // Crea RooHistPdf da RooDataHist
        RooHistPdf pdfSumL2_Boron("pdfSumL2_Boron", "pdfSumL2_Boron", x, dataHistSumL2_Boron );
        RooHistPdf pdfSumL2_Carbon("pdfSumL2_Carbon", "pdfSumL2_Carbon", x, dataHistSumL2_Carbon );
        // Normalizza i coefficienti in base al numero di entries degli istogrammi
        double L1_Boron   = z1.first.at(i)->Integral(0,z1.first.at(i)->FindBin(5.8),"width");
        double L2_Boron   = z1.second.at(i)->Integral(0,z1.second.at(i)->FindBin(6),"width");
        double L1_Carbon  = z2.first.at(i)->Integral(0,z2.first.at(i)->FindBin(7),"width");
        double L2_Carbon  = z2.second.at(i)->Integral(0,z2.second.at(i)->FindBin(7),"width");
        double nTot = L2_Boron+L2_Carbon;
        // Coefficienti normalizzati
        double frac;
        if (i>=0 && i<=2) frac = 0.3;
        if (i==3) frac = 0.38;
        if (i==4) frac = 0.44;
        if (i>=5 && i<=10) frac = 0.5;
        if (i>=11 && i<=23) frac = 0.45;
        std::pair<double,double> b = getCoefBoron(i);
        std::pair<double,double> c = getCoefCarbon(i);
        RooRealVar coefBoron("coefBoron", "Coefficient for Boron", L1_Boron,0.,L1_Boron); //L1_Boron,0.,L1_Boron*0.3
        RooRealVar coefCarbon("coefCarbon", "Coefficient for Carbon", (z1.first.at(i)->Integral("width") - L1_Boron)*frac, 0. , (z1.first.at(i)->Integral("width") - L1_Boron)*frac); //(0.28*L1_Boron-L2_Boron)*1.02,0.5,L1_Boron-L2_Carbon
        // Crea il modello somma
        RooRealSumPdf model("model", "model", RooArgList(pdfSumL2_Boron, pdfSumL2_Carbon), RooArgList(coefBoron,coefCarbon) );
        // Fitta il modello ai dati
        RooFitResult* r = model.fitTo(dataHistSumL1_Boron, RooFit::Minimizer("Minuit2", "Migrad"), Timer(true), Save(), PrintLevel(3)) ;
        RooPlot* frame = x.frame();
        frame->GetYaxis()->SetTitle("Events");
        TString title = hist1->GetTitle();
        frame->SetTitle(nucleus+" @ "+title );
        dataHistSumL1_Boron.plotOn(frame);
        model.plotOn(frame,Components(pdfSumL2_Boron), LineColor(kBlue));
        model.plotOn(frame,Components(pdfSumL2_Carbon), LineColor(kRed));
        model.plotOn(frame, LineColor(kGreen));    
        frame->Draw();
        canvasBin->Print(out2.Data() );

        //integral della RooAddPdf
        //auto pdf= model.pdf("pdf");
        x.setRange("range",0,charge+0.8);
        auto integralBoron = model.createIntegral(x,RooFit::NormSet(x), RooFit::Range("range"));
        auto integralCarbon= pdfSumL2_Carbon.createIntegral(x,RooFit::NormSet(x),RooFit::Range("range"));

        const RooArgList& params = r->floatParsFinal();
        RooRealVar* param_b = (RooRealVar*)params.find("coefBoron");
        RooRealVar* param_c = (RooRealVar*)params.find("coefCarbon");
        double bor = param_b->getVal();
        double car = param_c->getVal();
        double den = bor*z1.second.at(i)->Integral(0, z1.second.at(i)->FindBin(charge+0.8)) + car*z2.second.at(i)->Integral(0, z2.second.at(i)->FindBin(charge+0.8)); //conteggi fino al taglio
        double num = car*z2.second.at(i)->Integral(0, z2.second.at(i)->FindBin(charge+0.8)); //carbonii
        double h = (num/den);
        f_vs_r->SetPoint(f_vs_r->GetN(),Rigbins_LightIons[i], h);
    }   
    canvasBin->SetLogx();
    canvasBin->SetLogy(0);
    TAxis *axisY = f_vs_r->GetYaxis();
    TAxis *axisX = f_vs_r->GetXaxis();
    //axisX->SetRangeUser(3,6000);
    axisY->SetRangeUser(0,0.07);
    f_vs_r->Draw("ACP");
    canvasBin->Print(out2.Data() );
    canvasBin->Print(out3.Data() );
    TString fileOut = "../Fragmentation/BelowL1/Fractions/contamination_"+nucleus+".root";
    TFile *output = new TFile(fileOut.Data(), "RECREATE");
    output->WriteTObject(f_vs_r,"f_vs_r");
}

std::pair<double,double> getCoefBoron(int i ) {
  std::pair<double,double> p {};
  switch (i) {
  case 0:
    p.first = 0.1;
    p.second = 0.2;
    break;
  case 1:
    p.first = 8;
    p.second = 9;
    break;
  case 2:
    p.first = 8;
    p.second = 9;
    break;
  case 3:
    p.first = 6;
    p.second = 8;
    break;
  case 16:
    p.first = 3;
    p.second = 5;
    break;
  default:
    p.first = 7;
    p.second = 9;
  }
  return p;
}

std::pair<double,double> getCoefCarbon(int i ) {
  std::pair<double,double> p {};
  switch (i) {  
  case 0:
    p.first = 0.1;
    p.second = 0.7;
    break;
  case 1:
    p.first = 0.1;
    p.second = 0.8;
    break;
  case 2:
    p.first = 0.1;
    p.second = 0.9;
    break;
  case 3:
    p.first = 0.1;
    p.second = 0.7;
    break;
  case 4:
    p.first = 0.1;
    p.second = 0.7;
    break;
  case 5:
    p.first = 0.1;
    p.second = 0.5;
    break;
  case 6:
    p.first = 0.1;
    p.second = 0.4;
    break;
  case 7:
    p.first = 0.1;
    p.second = 0.6;
    break;
  case 8:
    p.first = 0.1;
    p.second = 0.4;
    break;
  case 9:
    p.first = 0.1;
    p.second = 0.35;
    break;
  case 10:
    p.first = 0.1;
    p.second = 0.25;
    break;
  case 11:
    p.first = 0.1;
    p.second = 0.4;
  case 12:
    p.first = 0.1;
    p.second = 0.25;
    break;
  case 13:
    p.first = 0.1;
    p.second = 0.25;
    break;
  case 14:
    p.first = 0.1;
    p.second = 0.15;
    break;
  case 15:
    p.first = 0.01;
    p.second = 0.25;
    break;
  case 16:
    p.first = 0.01;
    p.second = 0.08;
    break;
  case 17:
    p.first = 0.01;
    p.second = 0.25;
    break;
  case 18:
    p.first = 0.01;
    p.second = 0.25;
    break;
  
  default:
    p.first = 0.001;
    p.second = 0.25;
  }
  return p;
}