#include "binning.h"
#include "utils.h"
#include "TLatex.h"

struct HistogramPair {
    TH1D* da;
    TH1D* mc;
};

struct ChargeData {
    unsigned int charge;
    HistogramPair da_mc;
    HistogramPair ratio_fit;
};
int nEff = 8; //tof, track, trigger, l1, l1unb, track charge, daq, L1ChargeCut

void CanvasPartition(TCanvas *C, const Int_t Nx, const Int_t Ny, Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);
double XtoPad(double x);
double YtoPad(double x);
double legendSizeBasedOnDivisions(unsigned int n) {
	double l=0.1;
	switch(n) {
		case 3:
			l = 0.05;
			break;
		case 4:
			l = 0.05;
			break;
		case 5:
			l = 0.05;
			break;
		case 11:
		case 12:
		case 13:
			l=0.18;
			break;
	}
	return l;
}
void CanvasPartition(TCanvas *C, const Int_t Nx, const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin) {
   if (!C) return;

   // Setup Pad layout
   Float_t vSpacing = 0.0;
   Float_t vStep = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny;

   Float_t hSpacing = 0.0;
   Float_t hStep = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx;

   Float_t vposd, vposu, vmard, vmaru, vfactor;
   Float_t hposl, hposr, hmarl, hmarr, hfactor;

   for (Int_t i = 0; i < Nx; i++) {

      if (i == 0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr - hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx - 1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         if (hposr > 1.0) hposr = 1.0 - 1e-6;  // evita x=1.0
         hfactor = hposr - hposl;
         hmarl = 0.0;
         hmarr = rMargin / hfactor;
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr - hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }

      for (Int_t j = 0; j < Ny; j++) {

         if (j == 0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu - vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny - 1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            if (vposu > 1.0) vposu = 1.0 - 1e-6;  // evita y=1.0
            vfactor = vposu - vposd;
            vmard = 0.0;
            vmaru = tMargin / vfactor;
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu - vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }

         C->cd(0);

         TString name = Form("pad_%d_%d", i, j);
         if (auto oldPad = (TPad*)C->FindObject(name))
            delete oldPad;

         TPad *pad = new TPad(name, "", hposl, vposd, hposr, vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
         pad->Draw();
      }
   }
}

double XtoPad(double x)
{
   double xl, yl, xu, yu;
   gPad->GetPadPar(xl, yl, xu, yu);
   double pw = xu - xl;
   double lm = gPad->GetLeftMargin();
   double rm = gPad->GetRightMargin();
   double fw = pw - pw * lm - pw * rm;
   return (x * fw + pw * lm) / pw;
}
double YtoPad(double y)
{
   double xl, yl, xu, yu;
   gPad->GetPadPar(xl, yl, xu, yu);
   double ph = yu - yl;
   double tm = gPad->GetTopMargin();
   double bm = gPad->GetBottomMargin();
   double fh = ph - ph * bm - ph * tm;
   return (y * fh + bm * ph) / ph;
}

TString getEffName(int j);
TString getEffName2(int j);

template <typename T>
void print(std::vector<ChargeData> histograms, std::vector<unsigned int> charge_number,TString timePeriod);

std::vector<ChargeData> getHistos(std::vector<TFile*> files, std::vector<unsigned int> charge_number);

//Function definitions
std::vector<ChargeData> getHistos(std::vector<TFile *> files, std::vector<unsigned int> charge_number) {
	std::vector<ChargeData> histograms;
	auto empty = (TH1D*)hist_rig_highZ->Clone("empty");

	for (int j = 0; j < (int)files.size(); j++) {
		for (int i = 0; i < nEff; i++) {
			auto Da   = (TH1D*)files[j]->Get(da[i].Data());
			auto Mc   = (TH1D*)files[j]->Get(mc[i].Data());
			auto Damc = (TH1D*)files[j]->Get(damc[i].Data());
			auto Fit  = (TH1D*)files[j]->Get(fit[i].Data());

			if (Da && Mc && Damc && Fit) {
				histograms.push_back({ charge_number[j], {Da, Mc}, {Damc, Fit} });
			} 
			else if (Da && Fit && !(Mc && Damc)) { // DAQ && L1ChargeCut case
				histograms.push_back({ charge_number[j], {Da, Da}, {Da, Fit} });
			} 
			else if (!Da && !Fit && !Mc && !Damc) {
				histograms.push_back({ charge_number[j], {empty, empty}, {empty, empty} });
			} 
			else {
				std::cerr << "⚠️  [getHistos] Problem retrieving histograms for charge Z = " << charge_number[j]
				          << " (index j=" << j << ", i=" << i << ")\n"
				          << "    File: " << files[j]->GetName() << "\n"
				          << "    da[" << i << "]   = " << da[i]   << " -> " << (Da   ? "OK" : "MISSING") << "\n"
				          << "    mc[" << i << "]   = " << mc[i]   << " -> " << (Mc   ? "OK" : "MISSING") << "\n"
				          << "    damc[" << i << "] = " << damc[i] << " -> " << (Damc ? "OK" : "MISSING") << "\n"
				          << "    fit[" << i << "]  = " << fit[i]  << " -> " << (Fit  ? "OK" : "MISSING") << "\n"
				          << std::endl;
			}
		}
	}
	return histograms;
}

void print(std::vector<ChargeData> histograms, std::vector<unsigned int> charge_number,TString timePeriod, TString mcType) {
	TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TString prefix = getPrefixFile(charge_number);
	TString mc;
	mcType=="g" ? mc = "GlobalMC" : mc = "Mc";
	for (int j=0; j<nEff; j++) { 
		TCanvas *b = new TCanvas("","",2048,1280);
		if(j==0) b->SaveAs( Form("../output/%s_ratios_%s.pdf[", prefix.Data(),timePeriod.Data()) ,"RECREATE");		
		// Number of pads (columns × rows)
		int nx = charge_number.size();
		int ny = 2;
		
		// Margins
		Float_t lMargin = 0.065;
		Float_t rMargin = 0.001;
		Float_t bMargin = 0.075;
		Float_t tMargin = 0.03;

		double legTextSize = legendSizeBasedOnDivisions(charge_number.size());
		
		// Canvas setup
		CanvasPartition(b, nx, ny, lMargin, rMargin, bMargin, tMargin);
		TPad *pad[nx][ny];

		for (int i=0; i<charge_number.size(); i++) {	

			//----------top row		
			b->cd(0); // make sure the canvas is current	
			pad[i][1] = (TPad*)b->FindObject(Form("pad_%d_%d", i, 1));
			if (!pad[i][1]) {
				printf("❌ pad_%d_%d not found!\n", i, 1);
				continue;
			}
			pad[i][1]->Draw();
			pad[i][1]->cd();
			b->Update();
			// Calcola la larghezza relativa del pad rispetto a uno “standard”
			auto refPad = ((TPad*)b->FindObject("pad_0_1"));
			double padWidth = pad[i][1]->GetAbsWNDC(); // larghezza in coordinate canvas
			double refWidth = ((TPad*)b->FindObject("pad_0_1"))->GetAbsWNDC(); // per esempio, seconda colonna
			double totalWidth = pad[i][1]->GetAbsWNDC();
			double plotWidth = totalWidth * (1 - pad[i][1]->GetLeftMargin() - pad[i][1]->GetRightMargin());
			double refPlotWidth = refWidth * (1 - refPad->GetLeftMargin() - refPad->GetRightMargin());
			double scale = refPlotWidth / plotWidth;
			TString title = getIonName(i);			
			TString effName = getEffName(j);			
			TH2D *h1 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);			
			//setTitle(h1, "", "R (GV)", "", charge_number[i]);	
			//formatAxis(h1, charge_number.size());		
			adjustZoomY(h1,effName,charge_number[i]);			
			adjustZoomX(h1,effName);			
			 // Size factors
			Float_t xFactor = pad[0][1]->GetAbsWNDC() / pad[i][1]->GetAbsWNDC();
			Float_t yFactor = pad[0][1]->GetAbsHNDC() / pad[i][1]->GetAbsHNDC();
			//h1->SetTitleSize(0.21,"t");
			TH2D *hFrame = (TH2D *)h1->Clone(TString::Format("h_%d_%d", i, 1).Data());
			//hFrame->SetTitleSize(0.21,"t");
	
			// Format for y axis	
			hFrame->GetYaxis()->SetLabelSize(32);
			hFrame->GetYaxis()->SetLabelOffset(0.02);
			hFrame->GetYaxis()->SetTitleFont(63);
			hFrame->GetYaxis()->SetTitleSize(32);
			hFrame->GetYaxis()->SetTitleOffset(2);
	
			hFrame->GetYaxis()->CenterTitle();
	
			// TICKS Y Axis
			hFrame->GetYaxis()->SetTickLength(xFactor * 0.04 / yFactor);
	
			// Format for x axis
			hFrame->GetXaxis()->SetLabelSize(16);
			hFrame->GetYaxis()->SetLabelFont(63);
			hFrame->GetXaxis()->SetLabelFont(63);
			hFrame->GetXaxis()->SetLabelOffset(0.02);
			hFrame->GetXaxis()->SetTitleSize(16);
			//hFrame->GetXaxis()->SetTitleOffset(1);
			hFrame->GetXaxis()->CenterTitle();
	
			// TICKS X Axis
			hFrame->GetXaxis()->SetTickLength(yFactor * 0.03 / xFactor);
			hFrame->GetYaxis()->CenterTitle(true);
			hFrame->GetYaxis()->SetTitle( Form("%s efficiecy",effName.Data()) );
	
			// Draw cloned histogram with individual settings
			hFrame->Draw();		
			//h1->Draw();			
			
			formatTitle(h1,charge_number.size());	
			auto titleLegend = formatLegendToPad(XtoPad( i==0? 0.4-lMargin:0.4), YtoPad(0.99), XtoPad(0.6), YtoPad(1.09));
			titleLegend->AddEntry((TObject*)nullptr, Form("%s (Z=%d)",getIonPath(charge_number[i]).Data(),charge_number[i]),"" );
			titleLegend->SetTextSize(legTextSize*scale);
			titleLegend->Draw("SAME");
			auto a1 = formatLegendToPad(XtoPad(0.025), YtoPad(0.975), XtoPad(0.225), YtoPad(0.775));
			if (i!=charge_number.size()-1) h1->GetXaxis()->SetTitle(""); //remove repeated label			
			a1->AddEntry(histograms[j+nEff*i].da_mc.da, "Data","p");
			//setTitle(histograms[j+nEff*i].da_mc.da, effName, "R (GV)", "#varepsilon", charge_number[i]);					
			formatAxis(histograms[j+nEff*i].da_mc.da, charge_number.size());						
			formatTitle(histograms[j+nEff*i].da_mc.da,charge_number.size());				
			getColor(histograms[j+nEff*i].da_mc.da, charge_number[i], "dat");						
			adjustZoomY(histograms[j+nEff*i].da_mc.da, effName,charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].da_mc.da,charge_number.size());		
			adjustZoomX(histograms[j+nEff*i].da_mc.da, effName);
			formatMarkerSize(h1,charge_number.size());
			
			histograms[j+nEff*i].da_mc.da->Draw("P E1 SAME");		

			if (i!=charge_number.size()-1) histograms[j+nEff*i].da_mc.da->GetXaxis()->SetTitle(""); //remove repeated label		
			if (effName!="Daq" && effName !="L1 charge") a1->AddEntry(histograms[j+nEff*i].da_mc.mc, mc.Data(),"p");	
			//setTitle(histograms[j+nEff*i].da_mc.mc, effName, "R (GV)", "#varepsilon", charge_number[i]);	
			formatAxis(histograms[j+nEff*i].da_mc.mc, charge_number.size());	
			formatTitle(histograms[j+nEff*i].da_mc.mc,charge_number.size());	
			getColor(histograms[j+nEff*i].da_mc.mc, charge_number[i], "mc");
			adjustZoomY(histograms[j+nEff*i].da_mc.mc, effName, charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].da_mc.mc,charge_number.size());
			adjustZoomX(histograms[j+nEff*i].da_mc.mc, effName);
			formatMarkerSize(h1,charge_number.size());
			
			histograms[j+nEff*i].da_mc.mc->Draw("P E1 SAME");		

			if (i!=charge_number.size()-1) histograms[j+nEff*i].da_mc.mc->GetXaxis()->SetTitle(""); //remove repeated label
			gPad->SetLogx();	
			a1->SetTextSize(legTextSize*scale);	
			a1->Draw();

			if (i==0) {
				TLegend *y = formatLegendToPad(XtoPad(0.035), YtoPad(0.1), XtoPad(0.65), YtoPad(0.03));
				y->SetBorderSize(0);
 				y->SetFillStyle(0);  // Transparent background
				y->SetTextFont(62);
				y->SetTextSize(legTextSize*scale);
				y->AddEntry((TObject*)nullptr, effName, "");
				//y->Draw();
			}	

			//----------bottom row
		
			b->cd(0); // make sure the canvas is current
			pad[i][0] = (TPad*) b->FindObject(Form("pad_%d_%d", i, 0));
			if (!pad[i][0]) {
				printf("❌ pad_%d_%d not found!\n", 1, 0);
				continue;
			}
			pad[i][0]->SetTicks(1,0);
			pad[i][0]->Draw();
			pad[i][0]->cd();
			pad[i][0]->Update();
			TString title2 = getIonName(i);		
			TString effName2 = getEffName2(j);			
			TH2D *h2 = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);			
			if (i==charge_number.size() / 2) setTitle(h2, effName2, "R (GV)", "Data/MC", charge_number[i]);			
			h2->SetTitle("");
			//formatAxis(h2, charge_number.size());	
			formatTitle(h2,charge_number.size());			
			adjustZoomY(h2,effName2,charge_number[i]-1);			
			adjustZoomX(h2,effName2);	
			 // Size factors
			Float_t xFactor2 = pad[0][0]->GetAbsWNDC() / pad[i][0]->GetAbsWNDC();
			Float_t yFactor2 = pad[0][0]->GetAbsHNDC() / pad[i][0]->GetAbsHNDC();
	
			TH2D *hFrame2 = (TH2D *)h2->Clone(TString::Format("h2_%d_%d", i, 0).Data());
	
			// y axis range
			hFrame2->SetMinimum(0.0001); // do not show 0
			hFrame2->SetMaximum(1.2 * h1->GetMaximum());
	
			// Format for y axis
			hFrame2->GetYaxis()->SetLabelFont(43);
			hFrame2->GetYaxis()->SetLabelSize(32);
			hFrame2->GetYaxis()->SetLabelOffset(0.02);
			hFrame2->GetYaxis()->SetTitleFont(63);
			hFrame2->GetYaxis()->SetTitleSize(32);
			hFrame2->GetYaxis()->SetTitleOffset(2);
	
			hFrame2->GetYaxis()->CenterTitle();
	
			// TICKS Y Axis
			hFrame2->GetYaxis()->SetTickLength(xFactor2 * 0.04 / yFactor2);
	
			// Format for x axis
			hFrame2->GetXaxis()->SetLabelFont(43);
			hFrame2->GetXaxis()->SetLabelSize(32);
			hFrame2->GetXaxis()->SetLabelOffset(0.001);
			hFrame2->GetXaxis()->SetTitleFont(43);
			hFrame2->GetXaxis()->SetTitleSize(32);
			hFrame2->GetXaxis()->SetTitleOffset(1.5);
			//hFrame2->GetXaxis()->SetTitleOffset(1);
			hFrame2->GetXaxis()->CenterTitle();
	
			// TICKS X Axis
			hFrame2->GetXaxis()->SetTickLength(yFactor2 * 0.03 / xFactor2);
	
			// Draw cloned histogram with individual settings
			hFrame2->GetYaxis()->SetLabelFont(63);
			hFrame2->GetXaxis()->SetLabelFont(63);
			hFrame2->GetXaxis()->SetTitleFont(63);
			hFrame2->GetYaxis()->CenterTitle(true);
			if (effName!="Daq" && effName !="L1 charge")
				hFrame2->GetYaxis()->SetTitle("Ratio");
			else hFrame2->GetYaxis()->SetTitle("Spline");
			hFrame2->Draw();				
			//h2->Draw("");	

			auto a2 =  formatLegendToPad(XtoPad(0.025), YtoPad(0.975), XtoPad(0.225), YtoPad(0.775));
			if (i!=charge_number.size()-1) h1->GetXaxis()->SetTitle(""); //remove repeated label			
			//a2->AddEntry(histograms[j+nEff*i].ratio_fit.da , "Data/Mc");	
			//setTitle(histograms[j+nEff*i].ratio_fit.da, effName2, "R (GV)", "Data/Mc", charge_number[i]);		
			formatAxis(histograms[j+nEff*i].ratio_fit.da, charge_number.size()-1);			
			formatTitle(histograms[j+nEff*i].ratio_fit.da,charge_number.size());		
			getColor(histograms[j+nEff*i].ratio_fit.da, charge_number[i], "dat");			
			adjustZoomY(histograms[j+nEff*i].ratio_fit.da, effName2,charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].ratio_fit.da,charge_number.size());
			adjustZoomX(histograms[j+nEff*i].ratio_fit.da, effName);

			if (i!=charge_number.size()-1) histograms[j+nEff*i].ratio_fit.mc->GetXaxis()->SetTitle(""); //remove repeated label
			if (charge_number[i]!=15 || effName=="L1 charge") {
				const char* input_cstr = histograms[j+nEff*i].ratio_fit.mc->GetTitle();
                TString input(input_cstr);
                TString result;
                int chi2Pos = input.Index("Chi2/ndf:");
                if (chi2Pos != kNPOS) {
                    TString chi2Part = input(chi2Pos + 9, input.Length() - chi2Pos - 9);
                    int commaPos = chi2Part.Index(",");
                    TString chi2ValueStr = TString(chi2Part(0, commaPos));
                    chi2ValueStr.Strip(TString::kBoth);
                    double chi2Value = chi2ValueStr.Atof();
                    result.Form("= %.2f", chi2Value);
                }
				a2->AddEntry(histograms[j+nEff*i].ratio_fit.mc,Form("#chi^{2}/ndf%s",result.Data()),"L");
			}
			if (charge_number[i]==15 && effName!="L1 charge")
				a2->AddEntry(histograms[j+nEff*i].ratio_fit.mc,"Int", "LC");
			//setTitle(histograms[j+nEff*i].ratio_fit.mc, effName2, "R (GV)", "Data/Mc", charge_number[i]);
			formatAxis(histograms[j+nEff*i].ratio_fit.mc, charge_number.size()-1);
			formatTitle(histograms[j+nEff*i].ratio_fit.mc,charge_number.size());
			getColor(histograms[j+nEff*i].ratio_fit.mc, charge_number[i], "mc");
			adjustZoomY(histograms[j+nEff*i].ratio_fit.mc, effName2, charge_number[i]);
			formatMarkerSize(histograms[j+nEff*i].ratio_fit.mc,charge_number.size());

		 	TLegend *y = formatLegendToPad(XtoPad(0.025), YtoPad(0.975), XtoPad(0.225), YtoPad(0.775));
			y->SetBorderSize(0);
			y->SetFillStyle(0);   // Transparent background
			y->SetTextFont(62);
			y->SetTextSize(legTextSize*scale); // Large, consistent text
			y->AddEntry(histograms[j+nEff*i].ratio_fit.da, "Data/Mc","p"); 
			//y->Draw();

			if (charge_number[i]!=0) histograms[j+nEff*i].ratio_fit.da->Draw("P E1 SAME");		
			histograms[j+nEff*i].ratio_fit.mc->SetFillStyle(0);
			histograms[j+nEff*i].ratio_fit.mc->DrawCopy("hist L SAME");
			histograms[j+nEff*i].ratio_fit.mc->SetFillStyle(1001);
			histograms[j+nEff*i].ratio_fit.mc->SetMarkerSize(0);
			histograms[j+nEff*i].ratio_fit.mc->Draw("E3 L SAME");
			//b->RedrawAxis();
			adjustZoomX(histograms[j+nEff*i].ratio_fit.mc, effName);

			
			if (i!=charge_number.size()-1) histograms[j+nEff*i].ratio_fit.mc->GetXaxis()->SetTitle(""); //remove repeated label
			gPad->SetLogx();		
			a2->SetTextSize(legTextSize*scale);
			a2->Draw();
		}
		b->SaveAs( Form("../output/%s_ratios_%s.pdf", prefix.Data(), timePeriod.Data()) );
		if(j==nEff-1 ) b->SaveAs( Form("../output/%s_ratios_%s.pdf]", prefix.Data(), timePeriod.Data()) );
	}
}
TString getEffName(int j) {
	TString A;
	switch(j) {
		case 0:
		A="L1 pickup";
		break;
		case 1:
		A="Tof";
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
		A="Track charge";
		break;
		case 6:
		A="Daq";
		break;
		case 7:
		A = "L1 charge";
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
		A="Track charge Data/Mc";
		break;
		case 6:
		A="Daq";
		break;
		case 7:
		A = "L1 charge";
		break;
	}
	return A;
}
//-------------------------MAIN----------------------------
int main(int argc, char **argv) {
	if (argc < 4) {
		printf("Usage: \n");
		printf("%s <time> <mcType> <charge1> <charge2> ... <charge n> \n", argv[0]);
		return 1;
	}
	TString timePeriod = argv[1];
	TString mcType = argv[2];
    std::vector<unsigned int> charge_number;

	//Get the files
	auto files = getFilesFromUnfoldingFactor(argc, argv, charge_number,timePeriod,mcType);
	auto ratios = getHistos(files, charge_number);
	print(ratios, charge_number,timePeriod, mcType);
    setLogon(); 
}