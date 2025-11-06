#ifndef PLOTUTILITY_H
#define PLOTUTILITY_H
#include "includes.h"

std::vector<TString> da   = {"l1_da","tf_da","tr_da","tk_da","l1u_da","tkch_da","daq_da","l1ch"};
std::vector<TString> mc   = {"l1_mc","tf_mc","tr_mc","tk_mc","l1u_mc","tkch_mc"};
std::vector<TString> damc = {"damc_l1","damc_tf","damc_tr","damc_tk","damc_l1u","damc_tkch"};
std::vector<TString> fit  = {"final_damc_l1","final_damc_tf","final_damc_tr",
                             "final_damc_tk","final_damc_l1u","final_damc_tkch","final_daq","final_l1ch"};

std::vector<TString> ratios = {"l1_da","tf_da","tr_da","tk_da","l1u_da","tkch_da",
								               "l1_mc","tf_mc","tr_mc","tk_mc","l1u_mc","tkch_mc",
  "damc_l1","damc_tf","damc_tr","damc_tk","damc_l1u","damc_tkch",
	"final_damc_l1","final_damc_tf","final_damc_tr","final_damc_tk","final_damc_l1u","final_damc_tkch","final_daq","final_l1ch"};	

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
  case 18:
    A = "Ar";
    break;
  case 19:
    A = "K";
    break;
  case 20:
    A = "Ca";
    break;
  case 21:
    A = "Sc";
    break;
  case 22:
    A = "Ti";
    break;
  case 23:
    A = "V";
    break;
  case 24:
    A = "Cr";
    break;
  case 25:
    A = "Mn";
    break;
  case 26:
    A = "Fe";
    break;
  case 27:
    A = "Co";
    break;
  case 28:
    A = "Ni";
    break;
  }
  return A;
}
template <typename T>
void getColor(T* h, unsigned int charge, TString t) {
  // Set marker styles based on type
  h->SetMarkerStyle((t == "mc") ? 24 : 20);
  if (t == "gmc") h->SetMarkerStyle(26);
  // Define default alpha values
  double alpha = 0.15;
  if (t == "mc" || t == "gmc") alpha = 0.35;

  // Define color schemes for each charge
  struct ColorScheme {
    Color_t color_data;
    Color_t color_mc;
  };

  static const std::map<unsigned int, ColorScheme> colorMap = {
      {9,   {kBlue,      kBlue+2} },     // Fluorine
      {14,  {kGreen-1,   kGreen+1}},     // Silicon
      {15,  {kViolet,    kViolet+6}},    // Phosphorus
      {16,  {kOrange+2,  kOrange+6}},    // Sulfur
      {17,  {kYellow+2,  kYellow-6}},    // Chlorine (golden yellow, bright and clear)
      {18,  {kRed+1,     kRed-9}},       // Argon
      {19,  {kCyan+2,    kCyan-6}},      // Potassium (cyan, color-blind safe)
      {20,  {kAzure+2,   kAzure-4}},     // Calcium
      {21,  {kOrange-3,  kOrange+1}},    // Scandium (orange-brown, color-blind safe)
      {22,  {kTeal+1,    kTeal-6}},      // Titanium (teal, color-blind safe)
      {23,  {kGray+2,    kGray+1}},      // Vanadium (gray, neutral, distinct)
      {24,  {kSpring+4,  kSpring-6}},    // Chromium (yellow-green, bright but safe)
      {25,  {kPink-6,    kPink+3}},      // Manganese (muted rose, color-blind safe)
      {26,  {kMagenta+2, kMagenta-6}},   // Iron
      {999, {kBlack+2,   kBlack}}        // Generic black
  };


  auto it = colorMap.find(charge);
  if (it != colorMap.end()) {
    Color_t c = it->second.color_data;
    if (t == "mc" || t == "gmc") c = it->second.color_mc;
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    gStyle->SetTitleTextColor(c);
    setTransparentErrorBars(h, alpha, c);
  }
}
template <typename T>
void formatAxis(T &h, unsigned int nCanvasDivision) {
  gStyle->SetLabelFont(62,"XYZ");
  //h->GetXaxis()->CenterTitle(true);
  switch(nCanvasDivision) {
    case 1:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.008);
      h->GetXaxis()->SetLabelSize(0.052);
      h->GetXaxis()->SetTitleSize(0.045);
      h->GetXaxis()->SetTitleOffset(1.015);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.052);
      h->GetYaxis()->SetTitleSize(0.045);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(1.05);
      break;
    case 2:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.01);
      h->GetXaxis()->SetLabelSize(0.037);
      h->GetXaxis()->SetTitleSize(0.043);
      h->GetXaxis()->SetTitleOffset(0.8);
      //Y axis format
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetLabelSize(0.035);
      h->GetYaxis()->SetTitleSize(0.035);
      h->GetYaxis()->SetTitleOffset(1.4);
      break;
    case 3:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.005);
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.045);
      h->GetXaxis()->SetTitleOffset(1.02);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.045);
      h->GetYaxis()->SetTitleSize(0.032);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(1.70);
      break;
    case 4:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.025);
      h->GetXaxis()->SetLabelSize(0.045);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(0.6);
      //Y axis format
      h->GetYaxis()->SetLabelSize(0.038);
      h->GetYaxis()->SetTitleSize(0.032);
      h->GetYaxis()->SetLabelOffset(0.005);
      h->GetYaxis()->SetTitleOffset(1.5);
      break;
    case 5:
      //X axis format
      h->GetXaxis()->SetLabelOffset(-0.039);
      h->GetXaxis()->SetLabelSize(0.075);
      h->GetXaxis()->SetTitleSize(0.07);
      h->GetXaxis()->SetTitleOffset(0.65);
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
  h->SetTitleFont(62, "XYZ");
  
  // Make titles large and consistent
  double titleSize = 0.12;  // stronger visual presence
  if (nCanvasDivision > 10) titleSize = 0.11;
  if (nCanvasDivision > 15) titleSize = 0.10;
  double scaleX = gPad->XtoPad(1.0) - gPad->XtoPad(0.0);
  double scaleY = gPad->YtoPad(1.0) - gPad->YtoPad(0.0);
  h->SetTitleSize(titleSize* TMath::Min(scaleX, scaleY), "t");
}
TLegend* formatLegend(unsigned int nCanvasDivision, bool bottom = false) {
  double textSize = 0.06;
  if (nCanvasDivision > 10) textSize = 0.045;
  
  // Legend width
  double x1 = 0.60;
  double x2 = 0.95;
  
  // Top and bottom positioning
  double y1, y2;
  if (bottom) {
    y1 = 0.15;  // low
    y2 = 0.28;
  } else {
    y1 = 0.72;  // below title but inside plot
    y2 = 0.88;
  }

  TLegend *a = new TLegend(x1, y1, x2, y2);
  a->SetTextSize(textSize);
  a->SetBorderSize(0);
  a->SetFillStyle(0);
  return a;
}
TLegend* formatLegendToPad(double x1, double y1, double x2, double y2) {
  double textSize = 0.09;
  TLegend *a = new TLegend(x1, y1, x2, y2);
  a->SetTextSize(textSize);
  a->SetBorderSize(0);
  a->SetFillStyle(0);
  return a;
}
TLegend* formatLegendBottom( unsigned int nCanvasDivision) {
  std::map<int, double> x1{{1, 0.58}, {2, 0.58}, {3, 0.58}, {4, 0.58}, {5, 0.48}, {6, 0.48}, {9, 0.48}, {13, 0.58}};
  std::map<int, double> y1{{1, 0.90}, {2, 0.90}, {3, 0.90}, {4, 0.90}, {5, 0.90}, {6, 0.90}, {9, 0.90}, {13, 0.90}};
  std::map<int, double> x2{{1, 1.00}, {2, 1.00}, {3, 1.00}, {4, 1.00}, {5, 1.00}, {6, 1.00}, {9, 1.00}, {13, 1.00}};
  std::map<int, double> y2{{1, 1.00}, {2, 1.00}, {3, 1.00}, {4, 0.80}, {5, 1.00}, {6, 1.00}, {9, 1.00}, {13, 1.00}};
  std::map<int, double> textSize{{1, 0.065}, {2, 0.04}, {3, 0.045}, {4, 0.04}, {5, 0.06}, {6, 0.06}, {9, 0.11}, {13, 0.12}};

  TLegend *a= new TLegend(x1[nCanvasDivision], y1[nCanvasDivision], x2[nCanvasDivision], y2[nCanvasDivision]);
  a->SetTextSize(textSize[nCanvasDivision]);
  a->SetBorderSize(0);
 	a->SetFillStyle(0);
  return a;

}
void formatLegendTextSize(TLegend *a,unsigned int nCanvasDivision) {
  switch(nCanvasDivision) {
    case 1:
     a->SetTextSize(0.045);
      break;
    case 2:
      a->SetTextSize(0.045);
      break;
    case 3:
      a->SetTextSize(0.045);
      break;
    case 4:
     a->SetTextSize(0.06);
      break;
    case 5:
     a->SetTextSize(0.06);
     break;
    case 6:
     a->SetTextSize(0.045);
     break;
    case 9:
     a->SetTextSize(0.045);
     break;
    default:
     a->SetTextSize(0.059);
     break;
  }
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
void formatMarkerSize(T *h, unsigned int nCanvasDivision) {
  switch(nCanvasDivision) {
    case 1:
      h->SetMarkerSize(2.2);
      break;
    case 2:
      h->SetMarkerSize(1.4);
      break;
    case 3:
      h->SetMarkerSize(1.75);
      break;
    case 4:
      h->SetMarkerSize(1.2);
      break;
    case 5:
      h->SetMarkerSize(1.);
      break;
    case 6:
      h->SetMarkerSize(1.);
      break;
    case 9:
      h->SetMarkerSize(1.3);
      break;
    default:
      h->SetMarkerSize(1.2);
      break;
  }
}
template <typename T>
void setTitle(T *h, TString title, TString xLabel, TString yLabel, unsigned int charge) {
  if (!h) {
    std::cerr << "Error: Histogram pointer is null in setTitle" << std::endl;
    return;
  }
  h->GetYaxis()->SetLabelFont(62);
  h->GetYaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelFont(62);
  h->GetXaxis()->SetTitleFont(62);
  TString nucleus = getIonPath(charge);
  TString t = Form("%s (Z=%d)", nucleus.Data() ,charge );
  h->SetTitle(t);
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->SetTitle(yLabel);
  h->GetXaxis()->SetTitle(xLabel);
  gStyle->SetTitleSize(1.2);
}
template <typename T>
void formatZoomY(T *h, TString nameEff) { 
  h->GetYaxis()->SetLabelFont(62);
  h->GetXaxis()->SetLabelFont(62);
  gPad->SetGridy();
	gStyle->SetGridColor(kGray+3);
  if (nameEff == "L1" )       h->GetYaxis()->SetRangeUser(0.56, 0.93);
  if (nameEff == "Tof" )      h->GetYaxis()->SetRangeUser(0.95, 1.05);
  if (nameEff == "Trigger" )  h->GetYaxis()->SetRangeUser(0.92,1.02);
  if (nameEff == "Track" )    h->GetYaxis()->SetRangeUser(0.,2.);
  if (nameEff == "Flux")      h->GetYaxis()->SetRangeUser(0.1,25);
  if (nameEff == "Tof spline")h->GetYaxis()->SetRangeUser(0.8,1);
  if (nameEff == "Fragm")     h->GetYaxis()->SetRangeUser(0.0000000,1);
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
  double xmin = 2.15;
  if (a == "L1 pickup Data/Mc") h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "L1 detection Data/Mc")h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "Tof Data/Mc")       h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "Track Data/Mc")     h->GetXaxis()->SetRangeUser(xmin,30);
  if (a == "Track charge Data/Mc")h->GetXaxis()->SetRangeUser(xmin,1000+499);
	if (a == "Trigger Data/Mc")   h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "L1 pickup")         h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "L1 detection")      h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "Tof")              h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "Inner tracker")     h->GetXaxis()->SetRangeUser(xmin,30);
  if (a == "Track charge")      h->GetXaxis()->SetRangeUser(xmin,1000+499);
	if (a == "Trigger")           h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "L1 charge")         h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "Unfolding Factor")  h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a == "Mc acceptance")     h->GetXaxis()->SetRangeUser(2.15,2101);
  if (a=="Spline mc acceptance")h->GetXaxis()->SetRangeUser(2.15,1000+499);
  if (a=="InnTrk")              h->GetXaxis()->SetRangeUser(xmin,16);
  if (a=="rate")                h->GetXaxis()->SetRangeUser(0.8,1000+499);
  if (a=="Flux")                h->GetXaxis()->SetRangeUser(0.8,1000+499);
  if (a=="unfFact")             h->GetXaxis()->SetRangeUser(xmin,1000+499);
  if (a=="Daq")                 h->GetXaxis()->SetRangeUser(xmin,1000+499);
}
std::pair<double,double> yAxisLimits(TString title, unsigned int charge) {
  double ymin,ymax;
  switch (charge) {
  case 14:
  case 15:
  case 16:
  default:
    if (title == "L1 pickup")      {ymin=0.95  ;ymax=1.02;}
    if (title == "Tof")           {ymin=0.91  ;ymax=1.041;}
    if (title == "Trigger")        {ymin=0.755 ;ymax=1.17;}    
    if (title == "Inner tracker")  {ymin=0.8  ;ymax=0.955;}   
    if (title == "L1 detection")   {ymin=0.711  ;ymax=0.831;}
    if (title == "Track charge")   {ymin=0.91  ;ymax=1.16;}  
    if (title == "Daq")            {ymin=0.975  ;ymax=1.005;}   
    if (title == "L1 charge")      {ymin=0.86  ;ymax=1.03;}   

    if (title == "L1 pickup Data/Mc")     {ymin=0.951;ymax=1.031;}     
    if (title == "Tof Data/Mc")           {ymin=0.931 ;ymax=1.041;}     
    if (title == "Track Data/Mc")         {ymin=0.935  ;ymax=1.061;}     
    if (title == "Trigger Data/Mc")       {ymin=0.9  ;ymax=1.224;} 
    if (title == "L1 detection Data/Mc")  {ymin=0.891;ymax=1.041;}  
    if (title == "Track charge Data/Mc")  {ymin=0.941;ymax=1.061;}  

    if (title == "Tof spline")     {ymin=0.95 ;ymax=1.05;}   
    if (title == "Data/Mc Tot")    {ymin=0.9  ;ymax=1.;}     
    if (title == "Total Data/Mc correction") {ymin=0.8 ;ymax=1.1;}   
    if (title == "Total Data/Mc correction overlap") {
      {ymin=0.87 ;ymax=1.145;}   
      if (charge==18) {ymin=0.9 ;ymax=1.1;}   
      if (charge==20) {ymin=0.9 ;ymax=1.1;}   
    }
    if (title == "Flux")            {ymin=0.0001 ;ymax=0.24;}    
    if (title == "Unfolding Factor"){ymin=0.5 ;ymax=1.3;}    
    if (title == "Mc acceptance")   {ymin=0.01 ;ymax=0.067;}   
    if (title == "Spline mc acceptance") {ymin=0.01 ;ymax=0.57;} 
    if (title == "rate")              {ymin=10e-5 ;ymax=50.5;} 
    if (title == "unfFact")           {ymin=0.45 ;ymax=1.65;} 
    break;
  }
  return std::make_pair(ymin, ymax);
}
template <typename T>
void adjustZoomY(T *h, TString title, unsigned int charge) {
  h->GetYaxis()->SetLabelFont(62);
  h->GetXaxis()->SetLabelFont(62);
  h->GetYaxis()->SetTitleFont(62);
  h->GetXaxis()->SetTitleFont(62);
  gPad->SetGridy();
	gStyle->SetGridColor(kGray);
  gStyle->SetGridStyle(1);
  gStyle->SetGridWidth(1);
  auto yLim = yAxisLimits(title,charge);
  h->GetYaxis()->SetRangeUser(yLim.first,yLim.second);
}
template <typename T>
void setFriendColor(T* h, int i) {
    if (i == 0) {h->SetLineColor(TColor::GetColor("#E69F00")); h->SetMarkerColor(TColor::GetColor("#E69F00"));} // orange
    if (i == 1) {h->SetLineColor(TColor::GetColor("#56B4E9")); h->SetMarkerColor(TColor::GetColor("#56B4E9"));} // sky blue
    if (i == 2) {h->SetLineColor(TColor::GetColor("#009E73")); h->SetMarkerColor(TColor::GetColor("#009E73"));} // bluish green
    if (i == 3) {h->SetLineColor(TColor::GetColor("#F0E442")); h->SetMarkerColor(TColor::GetColor("#F0E442"));} // yellow
    if (i == 4) {h->SetLineColor(TColor::GetColor("#0072B2")); h->SetMarkerColor(TColor::GetColor("#0072B2"));} // blue
    if (i == 5) {h->SetLineColor(TColor::GetColor("#D55E00")); h->SetMarkerColor(TColor::GetColor("#D55E00"));} // vermillion
    if (i == 6) {h->SetLineColor(TColor::GetColor("#CC79A7")); h->SetMarkerColor(TColor::GetColor("#CC79A7"));} // reddish purple
    if (i == 7) {h->SetLineColor(TColor::GetColor("#FC8D62")); h->SetMarkerColor(TColor::GetColor("#FC8D62"));} // salmon
    // Extended blind-friendly extras (for > 8 curves)
    if (i == 8)  {h->SetLineColor(TColor::GetColor("#999999")); h->SetMarkerColor(TColor::GetColor("#999999"));} // gray
    if (i == 9)  {h->SetLineColor(TColor::GetColor("#8DD3C7")); h->SetMarkerColor(TColor::GetColor("#8DD3C7"));} // light cyan-green
    if (i == 10) {h->SetLineColor(TColor::GetColor("#FFD92F")); h->SetMarkerColor(TColor::GetColor("#FFD92F"));} // bright yellow
    if (i == 11) {h->SetLineColor(TColor::GetColor("#A6D854")); h->SetMarkerColor(TColor::GetColor("#A6D854"));} // light green
    if (i == 12) {h->SetLineColor(TColor::GetColor("#E78AC3")); h->SetMarkerColor(TColor::GetColor("#E78AC3"));} // pink
    if (i == 13) {h->SetLineColor(TColor::GetColor("#000000")); h->SetMarkerColor(TColor::GetColor("#000000"));} // black

}
std::vector<TFile *> getFilesFromUnfoldingFactor(int argc, char **argv, std::vector<unsigned int> &chNum, TString timePeriod, TString mcType) {
	std::vector<TFile*> files;
  TString type;
	for (int i=3; i<argc; i++) { 
        TString ch = argv[i];
		chNum.push_back(atoi(argv[i]));
    if (mcType=="g") type="_global";
    else type ="";
        TString ionPath = getIonPath(atoi(argv[i]));
		TFile *file = new TFile( Form("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+"_UnfoldingFactor%s_%s.root",
                 type.Data(),timePeriod.Data()) );
		if (!file || file->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",argv[i]); }
		files.push_back(file);
	}
	return files;
}
std::vector<TFile*> getDistributionFiles(int argc, char **argv, 
                             std::vector<unsigned int> &charge_number, 
                             TString type, TString detector) {
    std::vector<TFile*> files;
    for (int i=2; i<argc; i++) { // skip mode argument
        charge_number.push_back(atoi(argv[i]));
        TString ionPath = getIonPath(atoi(argv[i]));
        TString filePath = Form("../IonsSelected/%s/Distributions/%s/%s.root",
                                ionPath.Data(), detector.Data(), type.Data());
        TFile *file = new TFile(filePath);
        if (!file || file->IsZombie()) 
            printf("Errore nell'apertura del file %s per carica %d \n",filePath.Data(),atoi(argv[i]));
        files.push_back(file);
        std::cout << "Get histogram from " << filePath << std::endl;
    }
    return files;
}
std::pair<int,int> getQRange(unsigned int charge) {
    int min = (int)charge - 10;
    int max = (int)charge + 10;
    return std::make_pair(min, max);
}
std::vector<std::pair<TFile*,TFile*>> getFilesFromUnfoldingFactorAll(int argc, char **argv, std::vector<unsigned int> &chNum, TString timePeriod) {
	std::vector<std::pair<TFile*,TFile*>>  files;
  TString type;
	for (int i=2; i<argc; i++) { 
    TString ch = argv[i];
		chNum.push_back(atoi(argv[i]));
    TString ionPath = getIonPath(atoi(argv[i]));
		TFile *file1 = new TFile( Form("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+"_UnfoldingFactor_%s.root",
                 timePeriod.Data()) );
    TFile *file2 = new TFile( Form("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+"_UnfoldingFactor_global_%s.root",
                 timePeriod.Data()) );
		if (!file1 || !file2 || file1->IsZombie() || file2->IsZombie()) { printf("Errore nell'apertura del file per carica %d \n",argv[i]); }
		files.push_back(std::make_pair(file1,file2));
	}
	return files;
}
std::pair<TPad*,TPad*> buildPadsForFitAndPull(double gap) {
    float totalHeight = 1.0;
    float bottomHeight = 0.30 * (1.0 - gap);
    float topHeight = 0.70 * (1.0 - gap);

    // Bottom pad: from y=0 to y=bottomHeight
    TPad* pad2 = new TPad("pad2", "Bottom Pad", 0, 0.0, 1, bottomHeight);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);  // For axis labels

    // Top pad: from y=bottomHeight + gap to y=1.0
    TPad* pad1 = new TPad("pad1", "Top Pad", 0, bottomHeight + gap, 1, 1.0);
    pad1->SetBottomMargin(0.1);  // Tight bottom margin

    return std::make_pair(pad1, pad2);

}
#endif