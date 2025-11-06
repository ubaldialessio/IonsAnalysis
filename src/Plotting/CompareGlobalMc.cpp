#include "binning.h"
#include "utils.h"

#include <TLatex.h>

int nEff = 6; //tof, track, trigger, l1, l1unb, track charge
std::vector<TString> eff = {"l1_da","tf_da","tr_da","tk_da","l1u_da","tkch_da",
								"l1_mc","tf_mc","tr_mc","tk_mc","l1u_mc","tkch_mc"};		

std::vector<TString> da_eff = {"l1_da","tf_da","tr_da","tk_da","l1u_da","tkch_da"};

std::vector<TString> mc_eff = {"l1_mc","tf_mc","tr_mc","tk_mc","l1u_mc","tkch_mc"};

TString plots[16] = {"l1_da","tf_da","tr_da","tk_da","l1u_da","tkch_da",
								"l1_mc","tf_mc","tr_mc","tk_mc","l1u_mc","tkch_mc"};

//Function declarations
template <typename T>
std::pair<std::vector<T>,std::vector<T> > getEffHist(std::pair<std::vector<TFile*>, std::vector<TFile*>>v, TString opt);

std::pair<std::vector<TFile*>, std::vector<TFile*>> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum);

template <typename T>
void print(std::pair<std::vector<T>, std::vector<T> > eff, TString opt, int Ncharges, std::vector<unsigned int> chNum);

TString getEffName(int j);

template <typename T>
void printSingleOption(std::pair<std::vector<T>,std::vector<T> > eff, TString opt, int Ncharges , std::vector<unsigned int> chNum);


int main(int argc, char **argv) {
	if (argc < 3) {
		printf("Usage: \n");
		printf("%s <charge1> <charge2> ... <charge n> <dat/mc/both> \n", argv[0]);
		return 1;
	}
	std::vector<unsigned int> charges, chNum;
	TString opt = argv[argc-1];
	//Get the files
	auto files = getFiles(charges, argc, argv, chNum);

	//Get the histograms (data,mc)
	
	auto eff = getEffHist<TH1D*>(files, opt);

	//print the histos on the canvas (only works for "both" option)
	if (opt=="both") { print(eff, opt, files.first.size(), chNum); }

	//print with the given option
	else { setLogon(); printSingleOption(eff, opt, files.first.size(), chNum); }
	return 1;
}

//Function definitions
std::pair<std::vector<TFile*>,std::vector<TFile*> >  getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum) {
	std::pair<std::vector<TFile*>,std::vector<TFile*> > files {};
	for (int i=1; i<argc-1; i++) { 
		chNum.push_back(atoi(argv[i]));
		TString ionPath = getIonPath(atoi(argv[i]));
		TFile *file1 = new TFile("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+"_UnfoldingFactor_13.5y.root" );
        TFile *file2 = new TFile("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+"_UnfoldingFactor_global_13.5y.root" );
		if (!file1 || file1->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",argv[i]); }
		if (!file2 || file2->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",argv[i]); }
		files.first.push_back(file1);
        files.second.push_back(file2);
	}
	return files;
}
template <typename T>
std::pair<std::vector<T>, std::vector<T>> getEffHist(std::pair<std::vector<TFile*>, std::vector<TFile*>> v, TString opt) {
    std::pair<std::vector<T>, std::vector<T>> effHist {};

    for (size_t i = 0; i < v.first.size(); ++i) {
        TFile* file1 = v.first[i];
        TFile* file2 = v.second[i];

        for (const auto& name : mc_eff) {
            if ((opt == "mc" || opt == "both") && std::find(mc_eff.begin(), mc_eff.end(), name) != mc_eff.end()) {
                if (file1) {
                    auto hist1 = (T)file1->Get(name.Data());
                    if (hist1) {
                        printf("Retrieving histogram from file1: %s\n", name.Data());
                        effHist.first.push_back(hist1);
                    }
                }
                if (file2) {
                    auto hist2 = (T)file2->Get(name.Data());
                    if (hist2) {
                        printf("Retrieving histogram from file2: %s\n", name.Data());
                        effHist.second.push_back(hist2);
                    }
                }
            }
        }
    }

    return effHist;
}


template <typename T>
void print(std::pair<std::vector<T>, std::vector<T>> eff, TString opt, int Ncharges, std::vector<unsigned int> chNum) {
	TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TString prefix = getPrefixFile(chNum);
	for (int j=0; j<nEff; j++) { 
		TCanvas *b = new TCanvas("","",2048,1280);
		if(j==0) b->SaveAs( Form("../output/%s_globalMc.pdf[", prefix.Data()) ,"RECREATE");		
		b->Divide(Ncharges,1,0,0);		
		for (int i=0; i<Ncharges; i++) {	
			auto a = formatLegend(chNum.size());			
			b->cd(i+1);			
			b->Update();
			TString title = getIonName(chNum[i]);		
			TString effName = getEffName(j);			
			TH2D *h = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);			
			setTitle(h, effName, "R (GV)", effName+" #bf{#varepsilon}", chNum[i]);			
			formatAxis(h, chNum.size());		
			formatTitle(h,chNum.size());			
			adjustZoomY(h,effName,chNum[i]);			
			adjustZoomX(h,effName);			

			h->Draw();			
			eff.first[j+nEff*i]->Draw("P E1 SAME");		
			if (i!=chNum.size()-1) h->GetXaxis()->SetTitle(""); //remove repeated label			
			a->AddEntry(eff.first[j+nEff*i], "Mc single");			
			setTitle(eff.first[j+nEff*i], effName, "R (GV)", "#varepsilon", chNum[i]);		
			formatAxis(eff.first[j+nEff*i], chNum.size());			
			formatTitle(eff.first[j+nEff*i],chNum.size());		
			getColor(eff.first[j+nEff*i], chNum[i], "dat");			
			adjustZoomY(eff.first[j+nEff*i], effName,chNum[i]);
			//adjustZoomX(eff.first[j+nEff*i], effName);
			formatMarkerSize(eff.first[j+nEff*i],chNum.size());

			eff.second[j+nEff*i]->Draw("P E1 SAME");
			if (i!=chNum.size()-1) eff.first[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			a->AddEntry(eff.second[j+nEff*i], "Mc global");
			setTitle(eff.second[j+nEff*i], effName, "R (GV)", "#varepsilon", chNum[i]);
			formatAxis(eff.second[j+nEff*i], chNum.size());
			formatTitle(eff.second[j+nEff*i],chNum.size());
			getColor(eff.second[j+nEff*i], chNum[i], "mc");
			adjustZoomY(eff.second[j+nEff*i], effName, chNum[i]);
			//adjustZoomX(eff.second[j+nEff*i], effName);
			formatMarkerSize(eff.second[j+nEff*i],chNum.size());

			if (i!=chNum.size()-1) eff.second[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			a->SetTextSize(0.045);
			gPad->SetLogx();		
			a->Draw();
		}
		b->SaveAs( Form("../output/%s_globalMc.pdf", prefix.Data()) );
		if(j==nEff-1 ) b->SaveAs( Form("../output/%s_globalMc.pdf]", prefix.Data()) );
	}
}
TString getEffName(int j) {
	TString A;
	switch(j) {
		case 0:
		A="L1 pickup";
		break;
		case 1:
		A="UTof";
		break;
		case 2:
		A="Trigger";
		break;
		case 3:
		A="Inner tracker";
		break;
		case 4:
		A="L1 detect";
		break;
		case 5:
		A="Track Charge";
		break;
	}
	return A;
}

template <typename T>
void printSingleOption(std::pair<std::vector<T>,std::vector<T> > eff, TString opt, int Ncharges , std::vector<unsigned int> chNum) {
	TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	style->cd();
	TString prefix = getPrefixFile(chNum)+"_"+opt;
	gStyle->SetOptStat(0);
	for (int j=0; j<nEff; j++) {
		TCanvas *b = new TCanvas("","",2048,1280);
		if (j==0) b->SaveAs( Form("../output/%s_eff.pdf[", prefix.Data()) ,"RECREATE");
		b->Divide(chNum.size(),1,0,0);
		for (int i=0; i<Ncharges; i++) {
			auto a = formatLegend(chNum.size());
			b->cd(i+1);
			b->Update();
			TString title = getIonName(chNum[i]);
			TString effName = getEffName(j);
			TH2D *h = new TH2D("", "", nRbins - 1, Rbins, 100, 0.2, 1.15);
			setTitle(h, effName, "R (GV)",effName+" #varepsilon", chNum[i]);
			formatAxis(h, chNum.size());
			formatTitle(h,chNum.size());
			formatZoomY(h,effName); //Contains also grid y option
			formatZoomX(h,effName);
			h->Draw();
			if (opt=="dat") {
				a->AddEntry(eff.first[j+nEff*i], "Data");
				setTitle(eff.first[j+nEff*i], effName, "R (GV)", effName+" #varepsilon", chNum[i]);
				formatAxis(eff.first[j+nEff*i], chNum.size());
				formatTitle(eff.first[j+nEff*i],chNum.size());
				formatZoomY(eff.first[j+nEff*i], effName);
				getColor(eff.first[j+nEff*i], chNum[i], opt);
				eff.first[j+nEff*i]->Draw("P E3 SAME");
				if (i!=chNum.size()-1) eff.first[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			}
			if (opt == "mc") {
				a->AddEntry(eff.second[j+nEff*i], "Mc");
				setTitle(eff.second[j+nEff*i], effName, "R (GV)", effName+" #varepsilon", chNum[i]);
				formatAxis(eff.second[j+nEff*i], chNum.size());
				formatTitle(eff.second[j+nEff*i],chNum.size());
				formatZoomY(eff.second[j+nEff*i], effName);
				getColor(eff.second[j+nEff*i], chNum[i], opt);
				eff.second[j+nEff*i]->Draw("P E3 SAME");
				if (i!=chNum.size()-1) eff.second[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			}
			if (i!=chNum.size()-1) h->GetXaxis()->SetTitle(""); //remove repeated label
			a->SetTextSize(0.045);
			a->Draw();
			gPad->SetLogx();
		}
		b->SaveAs( Form("../output/%s_eff.pdf", prefix.Data()) );
		if(j==nEff-1) b->SaveAs( Form("../output/%s_eff.pdf]", prefix.Data()) );
	}
}
