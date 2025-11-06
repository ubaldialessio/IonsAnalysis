#include "binning.h"
#include "utils.h"

const int nEff = 6; //tof, track, trigger, l1, l1 unb
std::vector<TString> eff = {"damc_l1","damc_tf","damc_tr","damc_tk","damc_l1u","damc_tkch",
	"final_damc_l1","final_damc_tf","final_damc_tr","final_damc_tk","final_damc_l1u","final_damc_tkch"};					
TString plots[6] = {"damc_l1","damc_tf","damc_tr","damc_tk","damc_l1u","damc_tkch"};	

std::vector<TFile *> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum);
TString getEffName(int j);

template <typename T>
std::vector<T> getEffHist(std::vector<TFile *> v);

template <typename T>
void printSingleOption(std::vector<T> eff, int Ncharges , std::vector<unsigned int> chNum);

template <typename T>
void print(std::pair<std::vector<T>, std::vector<T>> eff, int Ncharges, std::vector<unsigned int> chNum);

template <typename T>
std::pair<std::vector<T>,std::vector<T> > getEffAndFitHist(std::vector<TFile *> v, TString opt);

//-------------------------MAIN----------------------------
int main(int argc, char **argv) {
	if (argc < 3) {
		printf("Usage: \n");
		printf("%s <charge1> <charge2> ... <charge n> <overlapFit/no> \n", argv[0]);
		return 1;
	}
    std::vector<unsigned int> charges, chNum;
	TString opt = argv[argc-1];

	//Get the files
	auto files = getFiles(charges, argc, argv, chNum);
	if (opt == "no") {
    	auto eff = getEffHist<TH1D*>(files);
		printSingleOption(eff, files.size(), chNum);
	}
	if (opt == "overlapFit") {
		auto eff = getEffAndFitHist<TH1D*>(files, opt);
		print(eff, files.size(), chNum);
	}
    setLogon(); 
}

//Function definitions
std::vector<TFile *> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum) {
	std::vector<TFile*> files;
	for (int i=1; i<argc-1; i++) { 
		chNum.push_back(atoi(argv[i]));
		TString ionPath = getIonPath(atoi(argv[i]));
		TFile *file = new TFile("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+"_UnfoldingFactor.root" );
		if (!file || file->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",argv[i]); }
		files.push_back(file);
	}
	return files;
}

template <typename T>
std::vector<T> getEffHist(std::vector<TFile *> v) {
	std::vector<T> effHist {};
		for (const auto &file : v) {
			for (const auto& name : eff) {
				auto hist = (T)file->Get(name.Data());
				if (hist) {
					effHist.push_back(hist);
				}
			}
		}
		return effHist;
}

template <typename T>
void printSingleOption(std::vector<T> eff, int Ncharges , std::vector<unsigned int> chNum) {
	TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	style->cd();
	TString prefix = getPrefixFile(chNum);
	gStyle->SetOptStat(0);
	for (int j=0; j<nEff; j++) {
         
		TCanvas *b = new TCanvas("","",2048,1280);
		if (j==0) b->SaveAs( Form("../output/%s_corr.pdf[", prefix.Data()) ,"RECREATE");
		b->Divide(chNum.size(),1,0,0);
		for (int i=0; i<Ncharges; i++) {
             
			b->cd(i+1);
			b->Update();
			TString title = getIonName(chNum[i]);
			TString effName = getEffName(j);
			TH2D *h = new TH2D("", "", nRbins - 1, Rbins, 100, 0.2, 1.15);
             
			setTitle(h, effName, "R (GV)", "#varepsilon", chNum[i]);
			formatAxis(h, chNum.size());
			formatTitle(h,chNum.size());
			adjustZoomY(h,effName, chNum[i]); //Contains also grid y option
			adjustZoomX(h,effName);
			h->Draw();
			setTitle(eff[j+nEff*i], effName, "R (GV)", "#varepsilon", chNum[i]);
			formatAxis(eff[ j+nEff*i], chNum.size());
			formatTitle(eff[ j+nEff*i],chNum.size());
			adjustZoomY(eff[ j+nEff*i], effName, chNum[i]);
			getColor(eff[ j+nEff*i], chNum[i], "");
			eff[ j+nEff*i]->Draw("P E3 SAME");
			if (i!=chNum.size()-1) eff[ j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			if (i!=chNum.size()-1) h->GetXaxis()->SetTitle(""); //remove repeated label
			gPad->SetLogx();
		}
		b->SaveAs( Form("../output/%s_corr.pdf", prefix.Data()) );
		if(j==nEff-1) b->SaveAs( Form("../output/%s_corr.pdf]", prefix.Data()) );
	}
}
template <typename T>
void print(std::pair<std::vector<T>, std::vector<T>> eff, int Ncharges, std::vector<unsigned int> chNum) {
	TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TString prefix = getPrefixFile(chNum);
	for (int j=0; j<nEff; j++) { 
		TCanvas *b = new TCanvas("","",2048,1280);
		if(j==0) b->SaveAs( Form("../output/%s_corr.pdf[", prefix.Data()) ,"RECREATE");		
		b->Divide(Ncharges,1,0,0);		
		for (int i=0; i<Ncharges; i++) {	
			auto a = formatLegend(chNum.size());			
			b->cd(i+1);			
			b->Update();
			TString title = getIonName(chNum[i]);		
			TString effName = getEffName(j);			
			TH2D *h = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);			
			setTitle(h, effName, "R (GV)", effName, chNum[i]);			
			formatAxis(h, chNum.size());			
			formatTitle(h,chNum.size());			
			adjustZoomY(h,effName,chNum[i]);			
			adjustZoomX(h,effName);			
			h->Draw("");	

			if (i!=chNum.size()-1) h->GetXaxis()->SetTitle(""); //remove repeated label			
			a->AddEntry(eff.first[j+nEff*i], "Data/Mc");	
			setTitle(eff.first[j+nEff*i], effName, "R (GV)", "Data/Mc", chNum[i]);		
			formatAxis(eff.first[j+nEff*i], chNum.size());			
			formatTitle(eff.first[j+nEff*i],chNum.size());		
			getColor(eff.first[j+nEff*i], chNum[i], "dat");			
			adjustZoomY(eff.first[j+nEff*i], effName,chNum[i]);
			formatMarkerSize(eff.first[j+nEff*i],chNum.size());
			//adjustZoomX(eff.first[j+nEff*i], effName);

			if (i!=chNum.size()-1) eff.first[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			a->AddEntry(eff.second[j+nEff*i], "Spline");
			setTitle(eff.second[j+nEff*i], effName, "R (GV)", "Data/Mc", chNum[i]);
			formatAxis(eff.second[j+nEff*i], chNum.size());
			formatTitle(eff.second[j+nEff*i],chNum.size());
			getColor(eff.second[j+nEff*i], chNum[i], "mc");
			adjustZoomY(eff.second[j+nEff*i], effName, chNum[i]);
			//formatMarkerSize(eff.second[j+nEff*i],chNum.size());

			eff.first[j+nEff*i]->Draw("P E1 SAME");		
			eff.second[j+nEff*i]->SetFillStyle(0);
			eff.second[j+nEff*i]->DrawCopy("hist L SAME");
			eff.second[j+nEff*i]->SetFillStyle(1001);
			eff.second[j+nEff*i]->SetMarkerSize(0);
			eff.second[j+nEff*i]->Draw("E3 L SAME");
			b->RedrawAxis();
			//adjustZoomX(eff.second[j+nEff*i], effName);

			if (i!=chNum.size()-1) eff.second[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			a->SetTextSize(0.035);
			gPad->SetLogx();		
			a->Draw();
		}
		b->SaveAs( Form("../output/%s_corr.pdf", prefix.Data()) );
		if(j==nEff-1 ) b->SaveAs( Form("../output/%s_corr.pdf]", prefix.Data()) );
	}
}
template <typename T>
std::pair<std::vector<T>, std::vector<T>> getEffAndFitHist(std::vector<TFile *> v, TString opt) {
	std::pair<std::vector<T>, std::vector<T>> effHist {};
    for (const auto &file : v) {
        for (const auto & name : eff) {
			auto hist = (T)file->Get(name.Data());
			if (hist) {
				if ( !name.Contains("final") ) {
					printf("Retrieving histogram : %s \n",name.Data());
					effHist.first.push_back(hist);
				} else if ( name.Contains("final") ) {
					printf("Retrieving histogram : %s \n",name.Data());
					effHist.second.push_back(hist);
				}
			}
		}
	}
	return effHist;
}

TString getEffName(int j) {
	TString A;
	switch(j) {
		case 0:
		A="L1 pickup Data/Mc";
		break;
		case 1:
		A="Tof Data/Mc";
		break;
		case 2:
		A="Trigger Data/Mc";
		break;
		case 3:
		A="Track Data/Mc";
		break;
		case 4:
		A="L1 detection Data/Mc";
		break;
		case 5:
		A="Track Charge Data/Mc";
		break;
	}
	return A;
}