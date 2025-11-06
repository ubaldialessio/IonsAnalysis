#include "binning.h"
#include "utils.h"

const int nPlot = 10;



// Define the tuple structure with (plot, title, xLabel, yLabel)
using PlotInfo = std::tuple<TString, TString, TString, TString>;

// Populate the data in a vector of tuples
std::vector<PlotInfo> plotData = {
    {"damc_l1", "L1 Data/Mc", "R (GV)", ""},
    {"damc_tf", "Tof Data/Mc", "R (GV)", ""},
    {"damc_tr", "Trigger Data/Mc", "R (GV)", ""},
    {"damc_tk", "Track Data/Mc", "R (GV)", ""},
    {"l1_da", "Data L1", "R (GV)", "#varepsilon"},
    {"tf_da", "Data Tof", "R (GV)", "#varepsilon"},
    {"tr_da", "Data Trigger", "R (GV)", "#varepsilon"},
    {"tk_da", "Data Track", "R (GV)", "#varepsilon"},
    {"l1_mc", "Mc L1", "R (GV)", "#varepsilon"},
    {"tf_mc", "Mc Tof", "R (GV)", "#varepsilon"},
    {"tr_mc", "Mc Trigger", "R (GV)", "#varepsilon"},
    {"tk_mc", "Mc Track", "R (GV)", "#varepsilon"},
    {"acc", "Mc acceptance", "R (GV)", ""},
    {"spline_acc", "Spline mc acceptance", "R (GV)", ""},
    {"flux", "Flux", "R (GV)", "#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})"},
    {"fluxMultiplied", "Flux", "R (GV)", "#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})"},
    {"flux", "Flux", "R (GV)", "#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{-1})"},
    {"lvt_25", "Livetime", "R (GV)", "t (s)"},
    {"rate", "Rate", "R (GV)", "Rate (s^{-1} GV^{-1})"},
    {"counts", "Counts", "R (GV)", "Counts"},
    {"final_acc", "Total acceptance", "R (GV)", ""},
    {"final_damc_tot", "Data/Mc Tot", "R (GV)", ""},
    {"final_damc_l1", "L1 spline", "R (GV)", ""},
    {"final_damc_tf", "Tof spline", "R (GV)", ""},
    {"final_damc_tr", "Trigger spline", "R (GV)", ""},
    {"final_damc_tk", "Track spline", "R (GV)", ""},
    {"Total correction", "Total Data/Mc correction", "R (GV)", ""},
    {"unf_factor", "Unfolding Factor", "R(GV)", ""},
    {"sp_unf", "Spline Unfolding Factor", "R(GV)", ""},
    {"flux_model_0", "Flux fitting", "R(GV)", ""},
    {"pub_flux", "Published flux as model", "R(GV)", ""},
    {"rateMulti", "Rate Multiplied", "R (GV)", "Rate R^{2.7} (s^{-1} GV^{1.7})"}
};


// Array di plots che mi interessano
    TString plots[] = {
        "damc_l1","damc_tf",
        "damc_tr","damc_tk",
        "acc", "flux", 
        "lvt_25","rate","counts","final_acc"
    };
    TString title[] = {"L1 Data/Mc", "Tof Data/Mc", "Trigger Data/Mc","Track Data/Mc",
                        "Mc acceptance", "Flux", "Livetime", "Rate", "Counts", "Total acceptance"};
    TString xLabel[] = {"R (GV)","R (GV)","R (GV)","R (GV)",
                        "R (GV)","R (GV)","R (GV)","R (GV)","R (GV)", "R (GV)"};
    TString yLabel[] = {"","","","","","#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})", 
                        "t (s)", "Rate (s^{-1} GV^{-1})", "Counts"};


std::vector<TFile *> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum);
void addStuffForP(std::vector<TFile*> &v);
void printAll(std::vector<TFile *> f, std::vector<unsigned int> chNum);
void SaveAllPlots(std::vector<TFile *> f, std::vector<unsigned int> chNum);
template <typename T>
T* getObjectFromFile(TFile* file, const std::string& objName);
TFile CombineFiles();


int main(int argc, char **argv) {
	if (argc < 2) {
		printf("Usage: \n");
		printf("%s <charge1> <charge2> ... <charge n> \n", argv[0]);
		return 1;
	}
  std::vector<unsigned int> charges, chNum;
	TString opt = argv[argc];

	//Get the files
	auto files = getFiles(charges, argc, argv, chNum);
  //Add the single files for P (Z=15), namely acceptance and damc_corrections
  //if (std::find(charges.begin(), charges.end(), 15) == charges.end()) addStuffForP(files);

  //print all the histos
  setLogon();
  SaveAllPlots(files, chNum);

}

//Function definitions
std::vector<TFile *> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum) {
	std::vector<TFile*> files;
	for (int i=1; i<argc; i++) { 
		chNum.push_back(atoi(argv[i]));
		TString ionPath = getIonPath(atoi(argv[i]));
		TFile *file = new TFile("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+"_UnfoldingFactor.root");
		if (!file || file->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",charges[i-1]); }
		files.push_back(file);
	}
	return files;
}  
void addStuffForP(std::vector<TFile*> &v) {
  TFile *fileAcceptance  = new TFile("../IonsSelected/P/P_acceptance.root");
  TFile *fileCorrections = new TFile("../IonsSelected/P/P_damc_corrections.root");
  if (fileAcceptance) v.push_back(fileAcceptance);
  if (fileCorrections)v.push_back(fileCorrections);
}

template <typename T>
T* getObjectFromFile(TFile* file, const std::string& objName) {
    T* obj = dynamic_cast<T*>(file->Get(objName.c_str()));
    if (!obj) {
        std::cerr << "Error: Object " << objName << " not found or wrong type!" << std::endl;
        return nullptr;
    }
    return obj;
}

void SaveAllPlots(std::vector<TFile *> f, std::vector<unsigned int> chNum) {
  TCanvas *c = new TCanvas();
  gStyle->SetOptStat(0);
  c->SaveAs("../output/everything.pdf[");
  for (int j=0; j<f.size(); j++) {
    TKey *key;
    TIter next((TList *)f[j]->GetListOfKeys());
    while (key = (TKey *)next()) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        TString objName = key->GetName();  // Retrieve object name
        TString plotName,titleName,xName,yName;
        for (const auto& plot : plotData) { //search the correct tuple by name
            if (std::get<0>(plot) == objName) {
              plotName      = std::get<0>(plot);
              titleName = std::get<1>(plot);
              xName     = std::get<2>(plot);
              yName     = std::get<3>(plot);
            }
        }
              if (cl->InheritsFrom("TGraphAsymmErrors")) {
                  c->Clear();
                  TGraphAsymmErrors* k = getObjectFromFile<TGraphAsymmErrors>(f[j], key->GetName());
                  gPad->SetGridy();
                  formatAxis(k,1);
                  setTitle(k,titleName,xName,yName, chNum[j]);
                  getColor(k,chNum[j],"");
                  adjustLabelOffeset(k,0);
                  adjustZoomY(k,titleName, chNum[j]);
                  adjustZoomX(k,titleName);
                  formatAxis(k,1);
                  k->Draw();
                  gPad->SetLogx();
                  c->SaveAs("../output/everything.pdf");
              } else if (cl->InheritsFrom("TH1D")) {
                  TH1D* h = getObjectFromFile<TH1D>(f[j],key->GetName());
                  if (plotName == "pub_flux") gPad->SetLogy();
                  if (plotName == "rate" || plotName == "counts") gPad->SetLogy();
                  if (plotName == "Total correction") gPad->SetGridy();
                  gPad->SetGridy();
                  formatAxis(h,1);
                  setTitle(h,titleName,xName,yName, chNum[j]);
                  getColor(h,chNum[j],"");
                  adjustLabelOffeset(h,0);
                  adjustZoomY(h,titleName, chNum[j]);
                  adjustZoomX(h,titleName);
                  formatAxis(h,1);
                  h->Draw();
                  gPad->SetLogx();
                  c->SaveAs("../output/everything.pdf");
                  gPad->SetLogy(0);
              } else if (cl->InheritsFrom("TF1")) {
                  TF1* h = getObjectFromFile<TF1>(f[j],key->GetName());
                  if (plotName == "flux_model_0") gPad->SetLogy();
                  if (plotName == "rate" || plotName == "counts") gPad->SetLogy();
                  if (plotName == "Total correction") gPad->SetGridy();
                  gPad->SetGridy();
                  formatAxis(h,1);
                  setTitle(h,titleName,xName,yName, chNum[j]);
                  getColor(h,chNum[j],"");
                  adjustLabelOffeset(h,0);
                  adjustZoomY(h,titleName, chNum[j]);
                  adjustZoomX(h,titleName);
                  formatAxis(h,1);
                  h->Draw();
                  gPad->SetLogx();
                  c->SaveAs("../output/everything.pdf");
                  gPad->SetLogy(0);
              }
    }
  }
  c->SaveAs("../output/everything.pdf]");
}

void printAll(std::vector<TFile *> f, std::vector<unsigned int> chNum) {
  TCanvas *c = new TCanvas();
  gStyle->SetOptStat(0);
  c->SaveAs("../output/everything.pdf[");
  for (int j=0; j<f.size(); j++) {
    for (int i=0; i<nPlot; i++) {
      c->Update();
      if (plots[i].Contains("damc") ) {
        TGraphAsymmErrors *k = (TGraphAsymmErrors *)(f[j]->Get(plots[i].Data()) );
        if (k) {
          c->Clear();
          formatAxis(k,1);
          setTitle(k,title[i],xLabel[i],yLabel[i], chNum[j]);
          getColor(k,chNum[j],"");
          adjustLabelOffeset(k,0);
          adjustZoomY(k,title[i], chNum[j]);
          adjustZoomX(k,title[i]);
          formatAxis(k,1);
          k->Draw();
          gPad->SetLogx();
          c->SaveAs("../output/everything.pdf");
        }
      }
      else {
      	TH1D* h = (TH1D *)(f[j]->Get(plots[i].Data()) );
        if (h) {
          if (plots[i] == "flux") MultiplyByXPower(h,2.70);
          if (plots[i] == "rate" || plots[i] == "counts") gPad->SetLogy();
          formatAxis(h,1);
          setTitle(h,title[i],xLabel[i],yLabel[i], chNum[j]);
          getColor(h,chNum[j],"");
          adjustLabelOffeset(h,0);
          adjustZoomY(h,title[i],chNum[j]);
          adjustZoomX(h,title[i]);
          formatAxis(h,1);
          h->Draw();
          gPad->SetLogx();
          c->SaveAs("../output/everything.pdf");
          gPad->SetLogy(0);
        }
      }
		} //for
  }
  c->SaveAs("../output/everything.pdf]");
}

