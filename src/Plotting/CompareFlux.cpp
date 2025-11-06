#include "binning.h"
#include "utils.h"

#include <TLatex.h>

const int nHistPerFile = 2;
std::vector<TString> fluxName = {"flux","fluxMultiplied"};

std::vector<TFile *> getFiles(int argc, char **argv, std::vector<unsigned int> &chNum);
std::vector<TH1D  *> getFluxHist(std::vector<TFile *> v);
void print(std::vector<TH1D*> fluxes, std::vector<unsigned int> chNum);

int main(int argc, char **argv) {
    if (argc < 2) {
		printf("Usage: \n");
		printf("%s <charge1> <charge2> ... <charge n> \n", argv[0]);
		return 1;
	}
	std::vector<unsigned int> chNum;
	auto files  = getFiles(argc, argv, chNum);
    auto fluxes = getFluxHist(files);
    print(fluxes, chNum );
}

std::vector<TFile *> getFiles(int argc, char **argv, std::vector<unsigned int> &chNum) {
	std::vector<TFile*> files;
	for (int i=1; i<argc; i++) { 
        TString ch = argv[i];
		chNum.push_back(atoi(argv[i]));
		TFile *file = new TFile( Form("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/%s/Flux/%s.root",
			getIonName(atoi(ch)).Data() ,getIonName(atoi(ch)).Data() ));
		if (!file || file->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",argv[i]); }
		files.push_back(file);
	}
	return files;
}
std::vector<TH1D *> getFluxHist(std::vector<TFile *> v) {
	std::vector<TH1D*> FluxHist {};
    for (const auto &file : v) 
        for (const auto& name : fluxName) {
			TH1D* hist = (TH1D *)(file->Get(name.Data()) );
			if (hist) FluxHist.push_back(hist);
		}
	return FluxHist;
}
void print(std::vector<TH1D*> fluxes, std::vector<unsigned int> charge_number) {
	TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TString prefix = getPrefixFile(charge_number);
	for (int j=0; j<nHistPerFile; j++) { 
		TCanvas *b = new TCanvas("","",2048,1280);
		if (j==0) b->SaveAs( Form("../output/%s.pdf[", prefix.Data()) ,"RECREATE");	
		b->Divide(charge_number.size(),2.,0.,0.);		
		
		//aspect ratio
		int nx = charge_number.size();
		int ny = 2;
		float padW = 1.0 / nx;
		float padH = 1.0 / ny;
		// Loop over pads and set their size and margins
		for (int ix = 0; ix < nx; ++ix) {
			for (int iy = 0; iy < ny; ++iy) {
				int padNumber = ix + 1 + iy * nx;
				b->cd(padNumber);
				TPad *pad = (TPad*)gPad;
				
				float xlow = ix * padW;
				float ylow = 1.0 - (iy + 1) * padH;
				float xup;
				xup  = xlow + padW;
				if (ix==nx-1) xup  = xlow + 0.989*padW;
				float yup  = ylow + 0.989*padH;

				pad->SetPad(xlow, ylow, xup, yup);
			}
		}

		for (int i=0; i<charge_number.size(); i++) {	
		//----------top row
			auto a1 = formatLegend(charge_number.size());			
			b->cd(i+1);
			b->Update();
			TString title = getIonName(i);		
			TH2D *h1 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);			
			setTitle(h1, "", "R (GV)", "", charge_number[i]);	
			formatAxis(h1, charge_number.size());		
			adjustZoomY(h1,"",charge_number[i]);			
			adjustZoomX(h1,"s");			
			formatTitle(h1,charge_number.size());			
			h1->Draw();			
			
			
			if (i!=charge_number.size()-1) h1->GetXaxis()->SetTitle(""); //remove repeated label			
			a1->AddEntry(fluxes[j+nHistPerFile*i], "Data");
			setTitle(fluxes[j+nHistPerFile*i], "", "R (GV)", "", charge_number[i]);					
			formatAxis(fluxes[j+nHistPerFile*i], charge_number.size());						
			formatTitle(fluxes[j+nHistPerFile*i],charge_number.size());				
			getColor(fluxes[j+nHistPerFile*i], charge_number[i], "dat");						
			/*adjustZoomY(fluxes[j+nHistPerFile*i], effName,charge_number[i]);
			formatMarkerSize(fluxes[j+nHistPerFile*i],charge_number.size());		
			adjustZoomX(fluxes[j+nHistPerFile*i], effName);
			formatMarkerSize(h1,charge_number.size());*/
			
			fluxes[j+nHistPerFile*i]->Draw();		

			if (i!=charge_number.size()-1) fluxes[j+nHistPerFile*i]->GetXaxis()->SetTitle(""); //remove repeated label		
			gPad->SetLogx();		
			a1->Draw();

			/*if (i==0) {
				TLegend *y = new TLegend(0.1,0.,0.3,0.13);
				y->SetBorderSize(0);
 				y->SetFillStyle(0);  // Transparent background
				y->AddEntry((TObject*)nullptr, effName, "");
				y->SetTextFont(62);
				y->SetTextSize(0.06);
				y->Draw();
			}*/

		//----------bottom row
			
			auto a2 = formatLegendBottom(charge_number.size());		
			b->cd(i+1+charge_number.size());			
			b->Update();
			TString title2 = getIonName(i);		
			TH2D *h2 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);			
			if (i==charge_number.size() / 2) setTitle(h2, "", "R (GV)", "Data/MC", charge_number[i]);			
			h2->SetTitle("");
			formatAxis(h2, charge_number.size());			
			formatTitle(h2,charge_number.size());			
			adjustZoomY(h2,"flux",charge_number[i]-1);			
			adjustZoomX(h2,"flux");			
			h2->Draw("");	

			if (i!=charge_number.size()-1) h1->GetXaxis()->SetTitle(""); //remove repeated label			
			//a2->AddEntry(fluxes[j+nHistPerFile*i] , "Data/Mc");	
			//setTitle(fluxes[j+nHistPerFile*i], effName2, "R (GV)", "Data/Mc", charge_number[i]);		
			formatAxis(fluxes[j+nHistPerFile*i], charge_number.size()-1);			
			formatTitle(fluxes[j+nHistPerFile*i],charge_number.size());		
			getColor(fluxes[j+nHistPerFile*i], charge_number[i], "dat");			
			adjustZoomY(fluxes[j+nHistPerFile*i], "",charge_number[i]);
			formatMarkerSize(fluxes[j+nHistPerFile*i],charge_number.size());
			adjustZoomX(fluxes[j+nHistPerFile*i], "");

			if (i!=charge_number.size()-1) fluxes[j+nHistPerFile*i]->GetXaxis()->SetTitle(""); //remove repeated label

			gPad->SetLogx();		
			a2->Draw();
		}
		 b->SaveAs( Form("../output/%s.pdf", prefix.Data()));
		 if(j==nHistPerFile-1) b->SaveAs( Form("../output/%s.pdf]", prefix.Data()));
	}






    /*TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	style->cd();
	gStyle->SetOptStat(0);
	TString prefix = getPrefixFile(chNum);
    for (int j=0; j<nHistPerFile; j++) {
        TCanvas *b = new TCanvas("","",2048,1280);
        if (j==0) b->SaveAs( Form("../output/%s.pdf[", prefix.Data()) ,"RECREATE");
        b->Divide(chNum.size(),1e-15,1e-15);		
        for (int i=0; i<chNum.size(); i++) {
			b->cd(i+1);
			b->Update();
			TString fluxTitle  = getIonName(chNum[i])+" flux";
            setTitle(fluxes[j+nHistPerFile*i], fluxTitle, "R (GV)", "", chNum[i]);
            if (j==0) gPad->SetLogy();
			formatAxis(fluxes[j+nHistPerFile*i], chNum.size());
			formatTitle(fluxes[j+nHistPerFile*i],chNum.size());
			getColor(fluxes[j+nHistPerFile*i],chNum[i], "");
			fluxes[j+nHistPerFile*i]->Draw();
			gPad->SetGridy();
            gPad->SetLogx();
        }
        b->SaveAs( Form("../output/%s.pdf", prefix.Data()));
        if(j==nHistPerFile-1) b->SaveAs( Form("../output/%s.pdf]", prefix.Data()));
    }*/
}