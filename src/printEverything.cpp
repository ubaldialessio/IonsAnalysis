#include "binning.h"
#include "utils.h"

const int nPlot = 10;
// Array di plots che mi interessano
    TString plots[] = {
        "damc_l1_7","damc_tf_7",
        "damc_tr_7","damc_tk_7",
        "acc_7", "flux_7", 
        "lvt_25","rate","counts","final_acc_7"
    };
    TString title[] = {"L1 Data/Mc", "Tof Data/Mc", "Trigger Data/Mc","Track Data/Mc",
                        "Acceptance", "Flux", "Livetime", "Rate", "Counts", "Total acceptance"};
    TString xLabel[] = {"R (GV)","R (GV)","R (GV)","R (GV)",
                        "R (GV)","R (GV)","R (GV)","R (GV)","R (GV)", "R (GV)"};
    TString yLabel[] = {"","","","","","#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})", "t (s)", "Rate (s^{-1} GV^{-1})", "Counts",""};

std::vector<TFile *> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum);
void MultiplyByXPower(TH1 *hh, double power);
void printAll(std::vector<TFile *> f, std::vector<unsigned int> chNum);

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

  //print all the histos
  setLogon();
  printAll(files, chNum);

}

//Function definitions
std::vector<TFile *> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum) {
	std::vector<TFile*> files;
	for (int i=1; i<argc; i++) { 
		chNum.push_back(atoi(argv[i]));
		TString ionPath = getIonPath(atoi(argv[i]));
		TFile *file = new TFile("../IonsSelected/"+ionPath+"/"+ionPath+"_output.root" );
		if (!file || file->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",charges[i-1]); }
		files.push_back(file);
	}
	return files;
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
          adjustZoomY(k,title[i]);
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
          if (plots[i] == "flux_7") MultiplyByXPower(h,2.70);
          if (plots[i] == "rate" || plots[i] == "counts") gPad->SetLogy();
          formatAxis(h,1);
          setTitle(h,title[i],xLabel[i],yLabel[i], chNum[j]);
          getColor(h,chNum[j],"");
          adjustLabelOffeset(h,0);
          adjustZoomY(h,title[i]);
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


void MultiplyByXPower(TH1 *hh, double power) {
  for (unsigned int ibin = 0; ibin < hh->GetNbinsX(); ibin++) {
    hh->SetBinContent(ibin + 1, hh->GetBinContent(ibin + 1) * pow(hh->GetBinCenter(ibin + 1), power));
    hh->SetBinError(ibin + 1, hh->GetBinError(ibin + 1) * pow(hh->GetBinCenter(ibin + 1), power));
  }
}

