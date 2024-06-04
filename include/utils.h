#ifndef UTILS_H
#define UTILS_H
#include "includes.h"

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
  }
  return A;
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
#endif
