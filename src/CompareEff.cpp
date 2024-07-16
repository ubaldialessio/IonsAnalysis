#include "binning.h"
#include "utils.h"

#include <TLatex.h>

const int nEff = 4; //tof, track, trigger, l1
std::vector<TString> eff = {"dataEff_l1","dataEff_tf","dataEff_tr","dataEff_tk",
								"mcEff_l1","mcEff_tf","mcEff_tr","mcEff_tk"};								
TString plots[8] = {"dataEff_l1","dataEff_tf","dataEff_tr","dataEff_tk",
								"mcEff_l1","mcEff_tf","mcEff_tr","mcEff_tk"};

//Function declarations
template <typename T>
std::pair<std::vector<T>,std::vector<T> > getEffHist(std::vector<TFile *> v, TString opt);

std::vector<TFile *> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum);

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
	
	auto eff = getEffHist<TGraphAsymmErrors* >(files, opt);

	//print the histos on the canvas (only works for "both" option)
	if (opt=="both") { print(eff, opt, files.size(), chNum); }

	//print with the given option
	else { setLogon(); printSingleOption(eff, opt, files.size(), chNum); }
	return 1;
}

//Function definitions
std::vector<TFile *> getFiles(std::vector<unsigned int> charges,int argc, char **argv, std::vector<unsigned int> &chNum) {
	std::vector<TFile*> files;
	for (int i=1; i<argc-1; i++) { 
		chNum.push_back(atoi(argv[i]));
		TString ionPath = getIonPath(atoi(argv[i]));
		TFile *file = new TFile("../IonsSelected/"+ionPath+"/"+ionPath+"_output.root" );
		if (!file || file->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",argv[i]); }
		files.push_back(file);
	}
	return files;
}
template <typename T>
std::pair<std::vector<T>, std::vector<T>> getEffHist(std::vector<TFile *> v, TString opt) {
	std::pair<std::vector<T>, std::vector<T>> effHist {};
    for (const auto &file : v) {
        for (const auto& name : eff) {
			auto hist = (T)file->Get(name.Data());
			if (hist) {
				if ( (opt == "dat" || opt == "both") && name.Contains("dataEff")) {
					effHist.first.push_back(hist);
				} else if ( (opt == "mc" || opt == "both") && name.Contains("mcEff")) {
					effHist.second.push_back(hist);
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
		if(j==0) b->SaveAs( Form("../output/%s_eff.pdf[", prefix.Data()) ,"RECREATE");
		b->Divide(Ncharges,1,0,0);
		for (int i=0; i<Ncharges; i++) {
			auto a = formatLegend(chNum.size());
			b->cd(i+1);
			b->Update();
			TString title = getIonName(chNum[i]);
			TString effName = getEffName(j);
			TH2D *h = new TH2D("", "", nRbins - 1, Rbins, 100, 0.2, 1.15);
			setTitle(h, effName, "R (GV)", "#varepsilon", chNum[i]);
			formatAxis(h, chNum.size());
			formatTitle(h,chNum.size());
			formatZoomY(h,effName);
			formatZoomX(h,effName);
			h->Draw("AXIS");
			eff.first[j+nEff*i]->Draw("P E3 SAME");

			if (i!=chNum.size()-1) h->GetXaxis()->SetTitle(""); //remove repeated label
			a->AddEntry(eff.first[j+nEff*i], "Data");
			setTitle(eff.first[j+nEff*i], effName, "R (GV)", "#varepsilon", chNum[i]);
			formatAxis(eff.first[j+nEff*i], chNum.size());
			formatTitle(eff.first[j+nEff*i],chNum.size());
			getColor(eff.first[j+nEff*i], chNum[i], "dat");
			formatZoomY(eff.first[j+nEff*i], effName);
			formatZoomX(eff.first[j+nEff*i], effName);

			eff.second[j+nEff*i]->Draw("P E3 SAME");
			if (i!=chNum.size()-1) eff.first[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			a->AddEntry(eff.second[j+nEff*i], "Mc");
			setTitle(eff.second[j+nEff*i], effName, "R (GV)", "#varepsilon", chNum[i]);
			formatAxis(eff.second[j+nEff*i], chNum.size());
			formatTitle(eff.second[j+nEff*i],chNum.size());
			getColor(eff.second[j+nEff*i], chNum[i], "mc");
			formatZoomY(eff.second[j+nEff*i], effName);
			formatZoomX(eff.second[j+nEff*i], effName);

			if (i!=chNum.size()-1) eff.second[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			a->SetTextSize(0.045);
			gPad->SetLogx();		
			a->Draw();
		}
		b->SaveAs( Form("../output/%s_eff.pdf", prefix.Data()) );
		if(j==nEff-1) b->SaveAs( Form("../output/%s_eff.pdf]", prefix.Data()) );
	}
}
TString getEffName(int j) {
	TString A;
	switch(j) {
		case 0:
		A="L1";
		break;
		case 1:
		A="Tof";
		break;
		case 2:
		A="Trigger";
		break;
		case 3:
		A="Track";
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
			setTitle(h, effName, "R (GV)", "#varepsilon", chNum[i]);
			formatAxis(h, chNum.size());
			formatTitle(h,chNum.size());
			formatZoomY(h,effName); //Contains also grid y option
			formatZoomX(h,effName);
			h->Draw();
			if (opt=="dat") {
				a->AddEntry(eff.first[j+nEff*i], "Data");
				setTitle(eff.first[j+nEff*i], effName, "R (GV)", "#varepsilon", chNum[i]);
				formatAxis(eff.first[j+nEff*i], chNum.size());
				formatTitle(eff.first[j+nEff*i],chNum.size());
				formatZoomY(eff.first[j+nEff*i], effName);
				getColor(eff.first[j+nEff*i], chNum[i], opt);
				eff.first[j+nEff*i]->Draw("P E3 SAME");
				if (i!=chNum.size()-1) eff.first[j+nEff*i]->GetXaxis()->SetTitle(""); //remove repeated label
			}
			if (opt == "mc") {
				a->AddEntry(eff.second[j+nEff*i], "Mc");
				setTitle(eff.second[j+nEff*i], effName, "R (GV)", "#varepsilon", chNum[i]);
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
