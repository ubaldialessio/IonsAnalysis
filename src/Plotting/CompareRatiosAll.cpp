#include "binning.h"
#include "utils.h"
#include "TLatex.h"

struct ChargeData {
    unsigned int charge;
	TH1D* data;
	TH1D* smc;
	TH1D* gmc;
	TH1D* data_smc_ratio;
	TH1D* data_gmc_ratio;
	TH1D* spline_smc_ratio;
	TH1D* spline_gmc_ratio;
};

int nEff = 7; //tof, track, trigger, l1, l1unb, track charge, daq

TString getEffName(int j);
TString getEffName2(int j);

template <typename T>
void print(std::vector<ChargeData> histograms, std::vector<unsigned int> charge_number,TString timePeriod);

std::vector<ChargeData> getHistos(std::vector<std::pair<TFile*,TFile*>> files, std::vector<unsigned int> charge_number);

//Function definitions
std::vector<ChargeData> getHistos(std::vector<std::pair<TFile*,TFile*>> files, std::vector<unsigned int> charge_number) {
	std::vector<ChargeData> histograms;
	auto empty = (TH1D*)hist_rig_highZ->Clone("empty");
	for (int j=0; j<files.size(); j++) {
		TFile* fileS = files[j].first; 
		TFile* fileG = files[j].second;
		for (int i=0; i<nEff; i++) {
			auto Da   			  = (TH1D*)fileS->Get(da[i].Data());
			auto SMc  			  = (TH1D*)fileS->Get(mc[i].Data());
			auto GMc 			  = (TH1D*)fileG->Get(mc[i].Data());
			auto Data_SMc_ratio   = (TH1D*)fileS->Get(damc[i].Data());
			auto Data_GMc_ratio   = (TH1D*)fileG->Get(damc[i].Data());
			auto Spline_Smc_ratio = (TH1D*)fileS->Get(fit[i].Data());
			auto Spline_GMc_ratio = (TH1D*)fileG->Get(fit[i].Data());
			if (Da && SMc && GMc && Data_SMc_ratio && Data_GMc_ratio && Spline_Smc_ratio && Spline_GMc_ratio)
				histograms.push_back( { charge_number[j], Da, SMc, GMc, Data_SMc_ratio, Data_GMc_ratio, Spline_Smc_ratio, Spline_GMc_ratio} );
			if (Da && !(SMc && GMc) ) //This is for DAQ
				histograms.push_back( { charge_number[j], Da, Da, Da, Da, Da, Spline_Smc_ratio, Spline_GMc_ratio } );
			if (!Da && !Spline_Smc_ratio && !SMc && !Data_SMc_ratio)
				histograms.push_back( { charge_number[j], empty,empty,empty,empty,empty,empty,empty } );
			else std::cerr << "Something wrong in retrieving histograms\n";
		}
	}
	return histograms;
}
void print(std::vector<ChargeData> histograms, std::vector<unsigned int> charge_number,TString timePeriod) {
	TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TString prefix = getPrefixFile(charge_number);
	for (int j=0; j<nEff; j++) { 
		TCanvas *b = new TCanvas("","",2048,1280);
		if(j==0) b->SaveAs( Form("../output/%s_ratiosAll_%s.pdf[", prefix.Data(),timePeriod.Data()) ,"RECREATE");		
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
			TString effName = getEffName(j);			
			TH2D *h1 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);			
			setTitle(h1, "", "R (GV)", "", charge_number[i]);	
			formatAxis(h1, charge_number.size());		
			adjustZoomY(h1,effName,charge_number[i]);			
			adjustZoomX(h1,effName);			
			formatTitle(h1,charge_number.size());			
			h1->Draw();			
//----------Data
			if (i!=charge_number.size()-1) h1->GetXaxis()->SetTitle(""); //remove repeated label			
			a1->AddEntry(histograms[j+nEff*i].data, "Data");
			//setTitle(histograms[j+nEff*i].data, effName, "R (GV)", "#varepsilon", charge_number[i]);					
			formatAxis(histograms[j+nEff*i].data, charge_number.size());						
			formatTitle(histograms[j+nEff*i].data,charge_number.size());				
			getColor(histograms[j+nEff*i].data, charge_number[i], "dat");						
			adjustZoomY(histograms[j+nEff*i].data, effName,charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].data,charge_number.size());		
			adjustZoomX(histograms[j+nEff*i].data, effName);
			formatMarkerSize(h1,charge_number.size());
			
			histograms[j+nEff*i].data->Draw("P E1 SAME");		
			if (i!=charge_number.size()-1) histograms[j+nEff*i].data->GetXaxis()->SetTitle(""); //remove repeated label		
//----------Single MC
			if (effName!="Daq") a1->AddEntry(histograms[j+nEff*i].smc, "Mc");	
			//setTitle(histograms[j+nEff*i].smc, effName, "R (GV)", "#varepsilon", charge_number[i]);	
			formatAxis(histograms[j+nEff*i].smc, charge_number.size());	
			formatTitle(histograms[j+nEff*i].smc,charge_number.size());	
			getColor(histograms[j+nEff*i].smc, charge_number[i], "mc");
			adjustZoomY(histograms[j+nEff*i].smc, effName, charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].smc,charge_number.size());
			adjustZoomX(histograms[j+nEff*i].smc, effName);
			formatMarkerSize(h1,charge_number.size());
			
			histograms[j+nEff*i].smc->Draw("P E1 SAME");		
			if (i!=charge_number.size()-1) histograms[j+nEff*i].smc->GetXaxis()->SetTitle(""); //remove repeated label

//----------Global MC
			if (effName!="Daq") a1->AddEntry(histograms[j+nEff*i].gmc, "G Mc");	
			//setTitle(histograms[j+nEff*i].gmc, effName, "R (GV)", "#varepsilon", charge_number[i]);	
			formatAxis(histograms[j+nEff*i].gmc, charge_number.size());	
			formatTitle(histograms[j+nEff*i].gmc,charge_number.size());	
			getColor(histograms[j+nEff*i].gmc, charge_number[i], "gmc");
			adjustZoomY(histograms[j+nEff*i].gmc, effName, charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].gmc,charge_number.size());
			adjustZoomX(histograms[j+nEff*i].gmc, effName);
			formatMarkerSize(h1,charge_number.size());
			
			histograms[j+nEff*i].gmc->Draw("P E1 SAME");		
			if (i!=charge_number.size()-1) histograms[j+nEff*i].gmc->GetXaxis()->SetTitle(""); //remove repeated label

			gPad->SetLogx();		
			a1->Draw();

			if (i==0) {
				TLegend *y = new TLegend(0.1,0.,0.3,0.13);
				y->SetBorderSize(0);
 				y->SetFillStyle(0);  // Transparent background
				y->AddEntry((TObject*)nullptr, effName, "");
				y->SetTextFont(62);
				y->SetTextSize(0.06);
				y->Draw();
			}	

		//----------bottom row
			
			auto a2 = formatLegendBottom(charge_number.size());		
			b->cd(i+1+charge_number.size());			
			b->Update();
			TString title2 = getIonName(i);		
			TString effName2 = getEffName2(j);			
			TH2D *h2 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);			
			if (i==charge_number.size() / 2) setTitle(h2, effName2, "R (GV)", "Data/MC", charge_number[i]);			
			h2->SetTitle("");
			formatAxis(h2, charge_number.size());			
			formatTitle(h2,charge_number.size());			
			adjustZoomY(h2,effName2,charge_number[i]-1);			
			adjustZoomX(h2,effName2);			
			h2->Draw("");	
//----------Data/Single MC ratio
			if (i!=charge_number.size()-1) h1->GetXaxis()->SetTitle(""); //remove repeated label			
			//a2->AddEntry(histograms[j+nEff*i].data_smc_ratio , "Data/Mc");	
			//setTitle(histograms[j+nEff*i].data_smc_ratio, effName2, "R (GV)", "Data/Mc", charge_number[i]);		
			formatAxis(histograms[j+nEff*i].data_smc_ratio, charge_number.size()-1);			
			formatTitle(histograms[j+nEff*i].data_smc_ratio,charge_number.size());		
			getColor(histograms[j+nEff*i].data_smc_ratio, charge_number[i], "mc");			
			adjustZoomY(histograms[j+nEff*i].data_smc_ratio, effName2,charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].data_smc_ratio,charge_number.size());
			adjustZoomX(histograms[j+nEff*i].data_smc_ratio, effName);

//----------Data/Global MC ratio
			if (i!=charge_number.size()-1) h1->GetXaxis()->SetTitle(""); //remove repeated label			
			//a2->AddEntry(histograms[j+nEff*i].data_gmc_ratio , "Data/Mc");	
			//setTitle(histograms[j+nEff*i].data_gmc_ratio, effName2, "R (GV)", "Data/Mc", charge_number[i]);		
			formatAxis(histograms[j+nEff*i].data_gmc_ratio, charge_number.size()-1);			
			formatTitle(histograms[j+nEff*i].data_gmc_ratio,charge_number.size());		
			getColor(histograms[j+nEff*i].data_gmc_ratio, charge_number[i], "gmc");			
			adjustZoomY(histograms[j+nEff*i].data_gmc_ratio, effName2,charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].data_gmc_ratio,charge_number.size());
			adjustZoomX(histograms[j+nEff*i].data_gmc_ratio, effName);

//----------Spline Data/Single MC 
			if (i!=charge_number.size()-1) histograms[j+nEff*i].spline_smc_ratio->GetXaxis()->SetTitle(""); //remove repeated label
			if (charge_number[i]!=15) {
				const char* input_cstr = histograms[j+nEff*i].spline_smc_ratio->GetTitle();
                TString input(input_cstr);
                TString result;
                int chi2Pos = input.Index("Chi2/ndf:");
                if (chi2Pos != kNPOS) {
                    TString chi2Part = input(chi2Pos + 9, input.Length() - chi2Pos - 9);
                    int commaPos = chi2Part.Index(",");
                    TString chi2ValueStr = TString(chi2Part(0, commaPos));
                    chi2ValueStr.Strip(TString::kBoth);
                    double chi2Value = chi2ValueStr.Atof();
                    result.Form("S, #chi^{2}/ndf = %.2f", chi2Value);
                }
                a2->AddEntry(histograms[j+nEff*i].data_smc_ratio, result);
			}
			if (charge_number[i]==15)	
				a2->AddEntry(histograms[j+nEff*i].spline_smc_ratio,"Interpolation", "LC");
			//setTitle(histograms[j+nEff*i].spline_smc_ratio, effName2, "R (GV)", "Data/Mc", charge_number[i]);
			formatAxis(histograms[j+nEff*i].spline_smc_ratio, charge_number.size()-1);
			formatTitle(histograms[j+nEff*i].spline_smc_ratio,charge_number.size());
			getColor(histograms[j+nEff*i].spline_smc_ratio, charge_number[i], "mc");
			adjustZoomY(histograms[j+nEff*i].spline_smc_ratio, effName2, charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].spline_smc_ratio,charge_number.size());
//----------Spline Data/Global MC
			if (i!=charge_number.size()-1) histograms[j+nEff*i].spline_gmc_ratio->GetXaxis()->SetTitle(""); //remove repeated label
			if (charge_number[i]!=15) {
				const char* input_cstr = histograms[j+nEff*i].spline_gmc_ratio->GetTitle();
                TString input(input_cstr);
                TString result;
                int chi2Pos = input.Index("Chi2/ndf:");
                if (chi2Pos != kNPOS) {
                    TString chi2Part = input(chi2Pos + 9, input.Length() - chi2Pos - 9);
                    int commaPos = chi2Part.Index(",");
                    TString chi2ValueStr = TString(chi2Part(0, commaPos));
                    chi2ValueStr.Strip(TString::kBoth);
                    double chi2Value = chi2ValueStr.Atof();
                    result.Form("G, #chi^{2}/ndf = %.2f", chi2Value);
                }
                a2->AddEntry(histograms[j+nEff*i].data_gmc_ratio, result);
			}
			if (charge_number[i]==15)	
				a2->AddEntry(histograms[j+nEff*i].spline_gmc_ratio,"Interpolation", "LC");
			//setTitle(histograms[j+nEff*i].spline_gmc_ratio, effName2, "R (GV)", "Data/Mc", charge_number[i]);
			formatAxis(histograms[j+nEff*i].spline_gmc_ratio, charge_number.size()-1);
			formatTitle(histograms[j+nEff*i].spline_gmc_ratio,charge_number.size());
			getColor(histograms[j+nEff*i].spline_gmc_ratio, charge_number[i], "mc");
			adjustZoomY(histograms[j+nEff*i].spline_gmc_ratio, effName2, charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].spline_gmc_ratio,charge_number.size());

//----------Drawing Spline single
			if (charge_number[i]!=0) histograms[j+nEff*i].data_smc_ratio->Draw("P E1 SAME");				
			histograms[j+nEff*i].spline_smc_ratio->SetFillStyle(0);
			histograms[j+nEff*i].spline_smc_ratio->DrawCopy("hist L SAME");
			histograms[j+nEff*i].spline_smc_ratio->SetFillStyle(1001);
			histograms[j+nEff*i].spline_smc_ratio->SetMarkerSize(0);
			histograms[j+nEff*i].spline_smc_ratio->Draw("E3 L SAME");
			//b->RedrawAxis();
			adjustZoomX(histograms[j+nEff*i].spline_smc_ratio, effName);
//----------Drawing Spline global
			if (charge_number[i]!=0) histograms[j+nEff*i].data_gmc_ratio->Draw("P E1 SAME");				
			histograms[j+nEff*i].spline_gmc_ratio->SetFillStyle(0);
			histograms[j+nEff*i].spline_gmc_ratio->DrawCopy("hist L SAME");
			histograms[j+nEff*i].spline_gmc_ratio->SetFillStyle(1001);
			histograms[j+nEff*i].spline_gmc_ratio->SetMarkerSize(0);
			histograms[j+nEff*i].spline_gmc_ratio->Draw("E3 L SAME");
			//b->RedrawAxis();
			adjustZoomX(histograms[j+nEff*i].spline_gmc_ratio, effName);

			
			if (i!=charge_number.size()-1) histograms[j+nEff*i].spline_gmc_ratio->GetXaxis()->SetTitle(""); //remove repeated label
			gPad->SetLogx();		
			a2->Draw();
			
			if (i==0) {
				TLegend *y = new TLegend(0.1,1.,0.3,0.87);
				y->SetBorderSize(0);
 				y->SetFillStyle(0);  // Transparent background
				y->AddEntry((TObject*)nullptr, "Data/Mc", "");
				y->SetTextFont(62);
				y->SetTextSize(0.06);
				y->Draw();
			}	
		}
		b->SaveAs( Form("../output/%s_ratiosAll_%s.pdf", prefix.Data(), timePeriod.Data()) );
		if(j==nEff-1 ) b->SaveAs( Form("../output/%s_ratiosAll_%s.pdf]", prefix.Data(), timePeriod.Data()) );
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
		A="L1 detection";
		break;
		case 5:
		A="Track Charge";
		break;
		case 6:
		A="Daq";
		break;
	}
	return A;
}
TString getEffName2(int j) {
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
		case 6:
		A="Daq";
		break;
	}
	return A;
}
//-------------------------MAIN----------------------------
int main(int argc, char **argv) {
	if (argc < 3) {
		printf("Usage: \n");
		printf("%s <time> <charge1> <charge2> ... <charge n> \n", argv[0]);
		return 1;
	}
	TString timePeriod = argv[1];
    std::vector<unsigned int> charge_number;

	//Get the files
	auto files = getFilesFromUnfoldingFactorAll(argc, argv, charge_number,timePeriod);
	auto ratios = getHistos(files, charge_number);
	print(ratios, charge_number,timePeriod);
    setLogon(); 
}