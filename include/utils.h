#ifndef UTILS_H
#define UTILS_H
#include "includes.h"

template <typename T>
void setTransparentErrorBars(T hist, float alpha, int c);
std::vector<std::pair<int, int>> GetGraphCommonPoints(TGraph *gr1, TGraph *gr2);

TString getIonName(unsigned int charge) {
TString A  = "";
  switch (charge) {
  case 1:
    A = "Proton";
    break;
  case 2:
    A = "Helium";
    break;
  case 3:
    A = "Litium";
    break;
  case 4:
    A = "Beryllium";
    break;
  case 5:
    A = "Boron";
    break;
  case 6:
    A = "Carbon";
    break;
  case 7:
    A = "Nitrogen";
    break;
  case 8:
    A = "Oxygen";
    break;
  case 9:
    A = "Fluorine";
    break;
  case 10:
    A = "Neon";
    break;
  case 14:
    A = "Silicon";
    break;
  case 15:
    A = "Phosphorus";
    break;
  case 16:
   	A = "Sulfur";
    break;
  case 17:
    A = "Chlorine";
    break;
  }
  return A;
}
TString getIonPath(unsigned int charge) {
TString A  = "";
  switch (charge) {
  case 1:
    A = "Pr";
    break;
  case 2:
    A = "He";
    break;
  case 3:
    A = "Li";
    break;
  case 4:
    A = "Be";
    break;
  case 5:
    A = "B";
    break;
  case 6:
    A = "C";
    break;
  case 7:
    A = "N";
    break;
  case 8:
    A = "O";
    break;
  case 9:
    A = "F";
    break;
  case 10:
    A = "Ne";
    break;
  case 11:
    A = "Na";
    break;
  case 12:
    A = "Mg";
    break;
  case 13:
    A = "Al";
    break;
  case 14:
    A = "Si";
    break;
  case 15:
    A = "P";
    break;
  case 16:
    A = "S";
    break;
  case 17:
  	A = "Cl";
  	break;
  }
  return A;
}
TString getPrefixFile(std::vector<unsigned int> chargeNumber) {
  TString a;
  for (auto const &c : chargeNumber) {
    a+=getIonPath(c);
  }
  return a;
}
std::pair<int,int> getRange(unsigned int charge ) {
  std::pair<int,int> p {};
  switch (charge) {
    case 1:
    p.first = 0;
    p.second = 1;
    break;
  case 2:
    p.first = 0;
    p.second = 3;
    break;
  case 3:
    p.first = 1;
    p.second = 4;
    break;
  case 4:
    p.first = 2;
    p.second = 5;
    break;
  case 5:
    p.first = 3;
    p.second = 7;
    break;
  case 6:
    p.first = 4;
    p.second = 8;
    break;
  case 7:
    p.first = 5;
    p.second = 9;
    break;
  case 8:
    p.first = 6;
    p.second = 10;
    break;
  case 9:
    p.first = 7;
    p.second = 11;
    break;
  case 10:
    p.first = 8;
    p.second = 12;
    break;
  case 14:
    p.first = 12;
    p.second = 16;
    break;
  case 15:
    p.first = 13;
    p.second = 17;
    break;
  case 16:
   	p.first = 14;
    p.second = 18;
    break;
  }
  return p;
}
TStyle *effHistStyle (unsigned int charge) {
	auto MyStyle  = new TStyle("MyStyle","MyStyle");
	MyStyle->SetErrorX(0.00000001);
	MyStyle->SetFrameFillColor(0);
	MyStyle->SetCanvasColor(0);
  MyStyle->SetPadColor(0);
	MyStyle->SetStatBorderSize(1);
	MyStyle->SetCanvasBorderMode(0);
	MyStyle->SetTitleFillColor(0);
	MyStyle->SetOptStat(0);
	MyStyle->SetStatStyle(0);
	MyStyle->SetTitleStyle(0);
	MyStyle->SetCanvasBorderSize(0);
	MyStyle->SetFrameBorderMode(0);
	MyStyle->SetFrameBorderSize(0);
	MyStyle->SetLegendBorderSize(1);
	MyStyle->SetStatBorderSize(0);
	MyStyle->SetTitleBorderSize(0);
	//MyStyle->SetFillStyle(4000);
	MyStyle->SetNdivisions(14,"y");
	UInt_t NRGBs = 5;
	Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[5]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[5]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };	
	Int_t NCont=255;
	TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
	switch (charge) {
	  case 5: //Boron
	    MyStyle->SetLineColor(kRed+2);
	    MyStyle->SetMarkerColor(kRed+2);
	    MyStyle->SetHistLineColor(kRed+2);
	    MyStyle->SetTitleTextColor(kRed+2);
	    break;
	  case 6: //Carbon
	    MyStyle->SetLineColor(kBlue+2);
	    MyStyle->SetMarkerColor(kBlue+2);
	    MyStyle->SetHistLineColor(kBlue+2);
	    MyStyle->SetTitleTextColor(kBlue+2);
	    break;
	  case 8: //Oxygen
	    MyStyle->SetLineColor(kAzure-4);
	    MyStyle->SetMarkerColor(kAzure-4);
	    MyStyle->SetHistLineColor(kAzure-4);
	    MyStyle->SetTitleTextColor(kAzure-4);
	    break;
	  case 14://Silicon
	  	MyStyle->SetLineColor(kGreen-2);
	  	MyStyle->SetMarkerColor(kGreen-2);
	    MyStyle->SetHistLineColor(kGreen-2);
	 	  MyStyle->SetTitleTextColor(kGreen-2);
	  	break;
	  case 15://Silicon
	  	MyStyle->SetLineColor(kMagenta+3);
	  	MyStyle->SetMarkerColor(kMagenta+3);
	  	MyStyle->SetHistLineColor(kMagenta+3);
	  	MyStyle->SetTitleTextColor(kMagenta+3);
	  	break;
	  case 16://Sulfur
		MyStyle->SetLineColor(kOrange+2);
	  	MyStyle->SetMarkerColor(kOrange+2);
	  	MyStyle->SetHistLineColor(kOrange+2);
	  	MyStyle->SetTitleTextColor(kOrange+2);  
	    break;
	};
	MyStyle->SetTitleAlign(33);
	MyStyle->SetTitleX(.55);
	MyStyle->SetTitleY(0.97);
	MyStyle->SetTitleXSize(0.05);
  MyStyle->SetTitleYSize(0.05);
  MyStyle->SetTitleFont(62,"xyz");
	MyStyle->SetLabelSize(0.05,"XYZ");
	gROOT->ForceStyle();
	
	return MyStyle; 
}
void beautifyEffHisto(TH1D *h, TString t, unsigned int charge, int Ncharges) {//Specific for CompareEff.cpp
  switch (Ncharges) {
    case 1:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.01);
      h->GetXaxis()->SetLabelSize(0.07);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(0.4);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.047);
      h->GetYaxis()->SetTitleSize(0.08);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(0.4);
      h->GetYaxis()->SetTitle("#varepsilon");

      h->SetTitleFont(62,"XYZ");
      h->SetLabelFont(62,"XYZ");
      gStyle->SetTitleSize(0.065, "t");
      gStyle->SetTitleAlign(Ncharges*5.5);
      gStyle->SetTitleSize(0.07, "t");
      gStyle->SetTitleY(0.92);
      gStyle->SetLegendTextSize(0.06); 
      break;
    case 2:
      gStyle->SetTitleAlign(Ncharges*7.5);
      break;
    case 3:
      gStyle->SetTitleAlign(Ncharges*7.5);
      gStyle->SetTitleSize(0.07, "t");
      h->GetXaxis()->SetLabelOffset(-0.038);
      break;
    case 4:
      gStyle->SetTitleAlign(Ncharges*5.5);
      gStyle->SetTitleSize(0.07, "t");
      gStyle->SetTitleY(0.92);
      gStyle->SetLegendTextSize(0.06); 
      break;
    case 5:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.052);
      h->GetXaxis()->SetLabelSize(0.07);
      h->GetXaxis()->SetTitleSize(0.07);
      h->GetXaxis()->SetTitleOffset(0.35);

      //Y axis format
      h->GetYaxis()->SetLabelSize(0.047);
      h->GetYaxis()->SetTitleSize(0.1);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(0.25);
      h->GetYaxis()->SetTitle("#varepsilon");

      h->SetTitleFont(62,"XYZ");
      h->SetLabelFont(62,"XYZ");
      gStyle->SetTitleSize(0.065, "t");
      gStyle->SetTitleAlign(Ncharges*5);
      gStyle->SetLegendTextSize(0.065); 
      break;
  }
  //gStyle->SetTitleX(.99);
  
  if (charge == 5) {
	  //Boron
	    h->SetLineColor(kRed+2);
	    h->SetMarkerColor(kRed+2);
      gStyle->SetTitleTextColor(kRed+2);
      setTransparentErrorBars(h,0.15,kRed+2);
      if (t=="mc") {
        h->SetLineColor(kRed-2);
	      h->SetMarkerColor(kRed-2);
        setTransparentErrorBars(h,0.35,kRed-2);
      }
	}
	if (charge == 6) { //Carbon
	    h->SetLineColor(kBlue+2);
	    h->SetMarkerColor(kBlue+2);
      gStyle->SetTitleTextColor(kBlue+2);
      setTransparentErrorBars(h,0.15,kBlue+2);
      if (t=="mc") {
        h->SetLineColor(kBlue-2);
        h->SetMarkerColor(kBlue-2);
        setTransparentErrorBars(h,0.35,kBlue-2);
      }
	}
	if (charge == 8) { //Oxygen
	    h->SetLineColor(kAzure-4);
	    h->SetMarkerColor(kAzure-4);
      gStyle->SetTitleTextColor(kAzure-4);
      setTransparentErrorBars(h,0.15,kAzure-4);
      if (t=="mc") {
        h->SetLineColor(kAzure-2);
        h->SetMarkerColor(kAzure-2);
        setTransparentErrorBars(h,0.35,kAzure-2);
        }
	}
	if (charge == 14) {  //Silicon
	  	h->SetLineColor(kGreen-5);
	  	h->SetMarkerColor(kGreen-5);
      gStyle->SetTitleTextColor(kGreen-5);
      setTransparentErrorBars(h,0.15,kGreen-5);
      if (t=="mc") {
        h->SetLineColor(kGreen-3);
        h->SetMarkerColor(kGreen-3);
        setTransparentErrorBars(h,0.35,kGreen-3);
      }
	}
	if (charge == 15) {  //Phosphorus
	  	h->SetLineColor(kViolet);
	  	h->SetMarkerColor(kViolet);
      gStyle->SetTitleTextColor(kViolet);
      setTransparentErrorBars(h,0.15,kViolet);
      if (t=="mc") {
        h->SetLineColor(kViolet+6);
        h->SetMarkerColor(kViolet+6);
        setTransparentErrorBars(h,0.35,kViolet+6);
        }
	}
	if (charge == 16) { //Sulfur
		  h->SetLineColor(kOrange+2);
	  	h->SetMarkerColor(kOrange+2); 
      gStyle->SetTitleTextColor(kOrange+2);
      setTransparentErrorBars(h,0.15,kOrange+2);
      if (t=="mc") {
        h->SetLineColor(kOrange+4);
        h->SetMarkerColor(kOrange+4); 
        setTransparentErrorBars(h,0.35,kOrange+4);
        }
	}

  if (t=="mc") {
    h->SetMarkerStyle(32);
    h->SetMarkerSize(1.2);
  }
  if (t=="dat") {
    h->SetMarkerStyle(8);
    h->SetMarkerSize(1.2);
  }
}
template <typename T>
void setTransparentErrorBars(T hist, float alpha, int c) { 
    // Colore delle barre d'errore
    hist->SetFillColorAlpha(c,alpha);
    //hist->SetFillStyle(3002); // Stile di riempimento con trasparenza
}
void setLogon() {
   gStyle->SetErrorX(0.00000001);
   gStyle->SetCanvasColor(0);
   gStyle->SetPadColor(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetStatBorderSize(1);
   gStyle->SetFrameFillColor(0);
   gStyle->SetTitleFillColor(0);
   gStyle->SetMarkerStyle(20);

   UInt_t NRGBs = 5;
   Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[5]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[5]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };

   Int_t NCont=255;
   TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
   gROOT->ForceStyle();

   gStyle->SetTitleYOffset(0.98);
   gStyle->SetTitleFont(62,"XYZ");
   gStyle->SetTitleXSize(0.05);
   gStyle->SetTitleYSize(0.05);
   gStyle->SetLabelSize(0.05,"XYZ");
}
template <typename T>
void getColor(T *h, unsigned int charge, TString t) {
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.2);
  if (t=="mc") {
    h->SetMarkerStyle(32);
    h->SetMarkerSize(1.2);
  }
  if (t=="dat") {
    h->SetMarkerStyle(8);
    h->SetMarkerSize(1.2);
  }
  if (charge == 5) { //Boron
	    h->SetLineColor(kRed+2);
	    h->SetMarkerColor(kRed+2);
      gStyle->SetTitleTextColor(kRed+2);
      setTransparentErrorBars(h,0.15,kRed+2);
      if (t=="mc") {
        h->SetLineColor(kRed-2);
	      h->SetMarkerColor(kRed-2);
        gStyle->SetTitleTextColor(kRed-2);
        setTransparentErrorBars(h,0.35,kRed-2);
      }
	}
	if (charge == 6) { //Carbon
	    h->SetLineColor(kBlue+2);
	    h->SetMarkerColor(kBlue+2);
      gStyle->SetTitleTextColor(kBlue+2);
      setTransparentErrorBars(h,0.15,kBlue+2);
      if (t=="mc") {
        h->SetLineColor(kBlue-2);
        h->SetMarkerColor(kBlue-2);
        gStyle->SetTitleTextColor(kBlue-2);
        setTransparentErrorBars(h,0.35,kBlue-2);
      }
	}
	if (charge == 8) { //Oxygen
	    h->SetLineColor(kAzure-4);
	    h->SetMarkerColor(kAzure-4);
      gStyle->SetTitleTextColor(kAzure-4);
      setTransparentErrorBars(h,0.15,kAzure-4);
      if (t=="mc") {
        h->SetLineColor(kAzure-2);
        gStyle->SetTitleTextColor(kAzure-2);
        h->SetMarkerColor(kAzure-2);
        setTransparentErrorBars(h,0.35,kAzure-2);
        }
	}
	if (charge == 14) {  //Silicon
	  	h->SetLineColor(kGreen-5);
	  	h->SetMarkerColor(kGreen-5);
      gStyle->SetTitleTextColor(kGreen-5);
      setTransparentErrorBars(h,0.15,kGreen-5);
      if (t=="mc") {
        h->SetLineColor(kGreen-3);
        h->SetMarkerColor(kGreen-3);
        gStyle->SetTitleTextColor(kGreen-3);
        setTransparentErrorBars(h,0.35,kGreen-3);
      }
	}
	if (charge == 15) {  //Phosphorus
	  	h->SetLineColor(kViolet);
	  	h->SetMarkerColor(kViolet);
      gStyle->SetTitleTextColor(kViolet);
      setTransparentErrorBars(h,0.15,kViolet);
      if (t=="mc") {
        h->SetLineColor(kViolet+6);
        h->SetMarkerColor(kViolet+6);
        gStyle->SetTitleTextColor(kViolet+6);
        setTransparentErrorBars(h,0.35,kViolet+6);
        }
	}
	if (charge == 16) { //Sulfur
		  h->SetLineColor(kOrange+2);
	  	h->SetMarkerColor(kOrange+2); 
      gStyle->SetTitleTextColor(kOrange+2);
      setTransparentErrorBars(h,0.15,kOrange+2);
      if (t=="mc") {
        h->SetLineColor(kOrange+6);
        h->SetMarkerColor(kOrange+6); 
        gStyle->SetTitleTextColor(kOrange+6);
        setTransparentErrorBars(h,0.35,kOrange+6);
        }
	}
}
template <typename T>
void formatAxis(T *h, unsigned int nCanvasDivision) {
  gStyle->SetLabelFont(62,"XYZ");
  switch(nCanvasDivision) {
    case 1:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.01);
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(1);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.047);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(1);
      break;
    case 2:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.01);
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.043);
      h->GetXaxis()->SetTitleOffset(0.5);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.047);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(0.85);
      break;
    case 3:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.02);
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(0.6);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.047);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(0.9);
      break;
    case 4:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.05);
      h->GetXaxis()->SetLabelSize(0.07);
      h->GetXaxis()->SetTitleSize(0.06);
      h->GetXaxis()->SetTitleOffset(0.45);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.047);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(0.70);
      break;
    case 5:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.052);
      h->GetXaxis()->SetLabelSize(0.07);
      h->GetXaxis()->SetTitleSize(0.07);
      h->GetXaxis()->SetTitleOffset(0.35);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.047);
      h->GetYaxis()->SetTitleSize(0.1);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(0.25);
      break;
  }
}
template <typename T>
void formatTitle(T *h, unsigned int nCanvasDivision) {
  gStyle->SetTitleFont(62,"XYZ");
  switch(nCanvasDivision) {
    case 1:
      gStyle->SetTitleSize(0.05, "t");
      gStyle->SetTitleAlign(23);
      break;
    case 2:
      gStyle->SetTitleSize(0.052, "t");
      gStyle->SetTitleAlign(23);
      break;
    case 3:
      gStyle->SetTitleSize(0.052, "t");
      gStyle->SetTitleAlign(23);
      break;
    case 4:
      gStyle->SetTitleSize(0.052, "t");
      gStyle->SetTitleAlign(23);
      break;
    case 5:
      gStyle->SetTitleSize(0.055, "t");
      gStyle->SetTitleAlign(23);
      break;
  }
}
TLegend* formatLegend( unsigned int nCanvasDivision) {
  std::map<int, double> x1{{1, 0.80}, {2, 0.80}, {3, 0.80}, {4, 0.80}, {5, 0.80}};
  std::map<int, double> y1{{1, 0.90}, {2, 0.90}, {3, 0.90}, {4, 0.90}, {5, 0.90}};
  std::map<int, double> x2{{1, 1.00}, {2, 1.00}, {3, 1.00}, {4, 1.00}, {5, 1.00}};
  std::map<int, double> y2{{1, 1.00}, {2, 1.00}, {3, 1.00}, {4, 1.00}, {5, 1.00}};

  TLegend *a= new TLegend(x1[nCanvasDivision], y1[nCanvasDivision], x2[nCanvasDivision], y2[nCanvasDivision]);
  return a;

}
template <typename T>
void adjustLabelOffeset(T *h, unsigned int nCanvasDivision) {
  switch(nCanvasDivision) {
    case 0:
      h->GetXaxis()->SetTitleOffset(0.45);
      h->GetXaxis()->SetLabelOffset(0.005);
      break;
  }
}
template <typename T>
void setTitle(T *h, TString title, TString xLabel, TString yLabel, unsigned int charge) {
  TString nucleus = getIonPath(charge);
  TString t = Form("%s (Z=%d) %s", nucleus.Data() ,charge, title.Data() );
  h->SetTitle(t);
  h->GetYaxis()->SetTitle(yLabel);
  h->GetXaxis()->SetTitle(xLabel);
}
template <typename T>
void formatZoomY(T *h, TString nameEff) { 
  h->GetYaxis()->SetLabelFont(62);
  h->GetXaxis()->SetLabelFont(62);
  gPad->SetGridy();
	gStyle->SetGridColor(kGray+3);
  if (nameEff == "L1" )       h->GetYaxis()->SetRangeUser(0.56, 0.93);
  if (nameEff == "Tof" )      h->GetYaxis()->SetRangeUser(0.82, 1.03);
  if (nameEff == "Trigger" )  h->GetYaxis()->SetRangeUser(0.91,1.06);
  if (nameEff == "Track" )    h->GetYaxis()->SetRangeUser(0.5,1.01);
  if (nameEff == "Flux")      h->GetYaxis()->SetRangeUser(0.1,25);
}
template <typename T>
void formatZoomX(T *h, TString nameEff) { 
  if (nameEff == "Track" )    h->GetXaxis()->SetRangeUser(0,30);
}
template <typename T>
void formatGridY(T *h, TString nameEff) { 
	gPad->SetGridy();
	gStyle->SetGridColor(kBlack+3);
}
template <typename T>
void adjustZoomX(T *h, TString a) {
	if (a== "Track")    h->GetXaxis()->SetRangeUser(1.8,30);
	if (a== "Trigger")  h->GetYaxis()->SetRangeUser(0.85,1.099);
	if (a== "Tof")      h->GetYaxis()->SetRangeUser(0.75,1.099);
  if (a == "L1 Data/Mc")        h->GetXaxis()->SetRangeUser(1,10000);
}
template <typename T>
void adjustZoomY(T *h, TString title) {
	if (title == "L1 Data/Mc")        h->GetYaxis()->SetRangeUser(0.,2);
	if (title == "Tof Data/Mc")       h->GetYaxis()->SetRangeUser(0.7,1.2);
	if (title == "Trigger Data/Mc")   h->GetYaxis()->SetRangeUser(0.9,1.2);
  if (title == "Track Data/Mc")     h->GetXaxis()->SetRangeUser(0,30);
}


TGraphErrors *DivideGraph(TGraphErrors *gr1, TGraphErrors *gr2, bool check) {

  std::vector<std::pair<int, int>> common_points = GetGraphCommonPoints(gr1, gr2);
  TGraphErrors *graph = new TGraphErrors();
  unsigned int npoints = common_points.size();

  for (unsigned int i = 0; i < npoints; i++) {
    int i1 = common_points[i].first, i2 = common_points[i].second;
    double ratio = gr1->GetY()[i1] / gr2->GetY()[i2];
    double ey = (gr1->GetEY()[i1] / gr1->GetY()[i1] + gr2->GetEY()[i2] / gr2->GetY()[i2]) * ratio;

    graph->SetPoint(i, gr1->GetX()[i1], ratio);
    graph->SetPointError(i, gr1->GetEX()[i1], ey);
  }

  return graph;
}
TGraphAsymmErrors *DivideAsymmGraph(TGraphAsymmErrors *gr1, TGraphAsymmErrors *gr2, bool check) {

  std::vector<std::pair<int, int>> common_points = GetGraphCommonPoints(gr1, gr2);
  TGraphAsymmErrors *graph = new TGraphAsymmErrors();
  unsigned int npoints = common_points.size();

  for (unsigned int i = 0; i < npoints; i++) {
    int i1 = common_points[i].first, i2 = common_points[i].second;
    if (gr2->GetY()[i2] == 0) {
      common_points.erase(common_points.begin() + i);
      i--;
      continue;
    }
    double ratio = gr1->GetY()[i1] / gr2->GetY()[i2];
    double eyh = (gr1->GetY()[i1] + gr1->GetEYhigh()[i1]) / (gr2->GetY()[i2] - gr2->GetEYlow()[i2]) - ratio;
    double eyl = ratio - (gr1->GetY()[i1] - gr1->GetEYlow()[i1]) / (gr2->GetY()[i2] + gr2->GetEYhigh()[i2]);

    graph->SetPoint(i, gr1->GetX()[i1], ratio);
    graph->SetPointError(i, gr1->GetEXlow()[i1], gr1->GetEXhigh()[i1], eyl, eyh);
  }

  return graph;
}
TGraphAsymmErrors *MultiplyAsymmGraph(TGraphAsymmErrors *gr1, TGraphAsymmErrors *gr2, bool check) {

  std::vector<std::pair<int, int>> common_points = GetGraphCommonPoints(gr1, gr2);
  TGraphAsymmErrors *graph = new TGraphAsymmErrors();
  unsigned int npoints = common_points.size();

  for (unsigned int i = 0; i < npoints; i++) {
    int i1 = common_points[i].first, i2 = common_points[i].second;
    double y = gr1->GetY()[i1] * gr2->GetY()[i2];
    double eyh = (gr1->GetY()[i1] + gr1->GetEYhigh()[i1]) * (gr2->GetY()[i2] + gr2->GetEYhigh()[i2]) - y;
    double eyl = y - (gr1->GetY()[i1] - gr1->GetEYlow()[i1]) * (gr2->GetY()[i2] - gr2->GetEYlow()[i2]);

    graph->SetPoint(i, gr1->GetX()[i1], y);
    graph->SetPointError(i, gr1->GetEXlow()[i1], gr1->GetEXhigh()[i1], eyl, eyh);
  }

  return graph;
}
std::vector<std::pair<int, int>> GetGraphCommonPoints(TGraph *gr1, TGraph *gr2) {

  std::vector<std::pair<int, int>> result;

  for (int ip1 = 0; ip1 < gr1->GetN(); ip1++) {
    double x1 = gr1->GetX()[ip1];
    for (int ip2 = 0; ip2 < gr2->GetN(); ip2++) {
      double x2 = gr2->GetX()[ip2];
      if (TMath::Abs(x1 - x2) < 1e-10) {
        result.push_back(std::make_pair(ip1, ip2));
      }
    }
  }
  Info("GetGraphCommonPoints", "Found %lu points in common", result.size());
  return result;
}
#endif
