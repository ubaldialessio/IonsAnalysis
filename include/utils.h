#ifndef UTILS_H
#define UTILS_H
#include "includes.h"
#include "PlotUtility.h"
#include "FluxUtility.h"
#include "SplineUtility.h"

template <typename T>
void setTransparentErrorBars(T hist, float alpha, int c);
std::vector<std::pair<int, int>> GetGraphCommonPoints(TGraph *gr1, TGraph *gr2);
void Rebin(TH1 *h, int n);

TString getPrefixFile(std::vector<unsigned int> chargeNumber) {
  TString a;
  for (auto const &c : chargeNumber) {
    a+=getIonPath(c);
  }
  return a;
}
std::pair<int,int> getRange(unsigned int charge ) {
  std::pair<double,double> p {};
  switch (charge) {
  case 5:
    p.first = 3;
    p.second = 8;
    break;
  case 6:
    p.first = 4;
    p.second = 9;
    break;
  case 8:
    p.first = 6;
    p.second = 11;
    break;
  case 14:
    p.first = 12;
    p.second = 17;
    break;
  case 15:
    p.first = 14;
    p.second = 16;
  case 16:
    p.first = 15;
    p.second = 16;
  case 17:
    p.first = 16;
    p.second = 17;
  default:
    p.first = 14;
    p.second = 19;
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
    h->SetMarkerStyle(21);
    h->SetMarkerSize(1.2);
  }
  if (t=="dat") {
    h->SetMarkerStyle(8);
    h->SetMarkerSize(1.22);
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
TH1D* divide(TH1* hpass, TH1* htotal, TString option) {
  auto ratio = (TH1D*)htotal->Clone(); ratio->Reset();
  if (option.Contains("TEfficiency")) {
    if (!TEfficiency::CheckConsistency(*hpass, *htotal)) {
      std::cerr << "[divide] Histograms are not consistent for TEfficiency!" << std::endl;
      return ratio;
    }
    TEfficiency eff(*hpass, *htotal);
    for (int i = 1; i <= hpass->GetNbinsX(); ++i) {
      if (!eff.GetEfficiency(i)) {
        ratio->SetBinContent(i, 0);
        ratio->SetBinError(i, 0);
        continue;
      }
      double efficiency = eff.GetEfficiency(i);
      double err_low = eff.GetEfficiencyErrorLow(i);
      double err_up  = eff.GetEfficiencyErrorUp(i);
      double err_sym = std::max(err_low, err_up);  // conservative symmetric error
      ratio->SetBinContent(i, efficiency);
      ratio->SetBinError(i, err_sym);
    }
    return ratio;
  }
  else if (option.Contains("efficiency")) {
      for (int i = 1; i <= hpass->GetNbinsX(); ++i) {
          double pass = hpass->GetBinContent(i);
          double tot = htotal->GetBinContent(i);
          pass = std::min(pass, tot);  // Prevent pass > tot
          if (tot > 0 && pass/tot!=1) {
              double eff = pass / tot;
              ratio->SetBinContent(i, eff);
              double err = sqrt(eff * (1 - eff) / tot);
              ratio->SetBinError(i, err);
          } else {
              ratio->SetBinContent(i, 0);
              ratio->SetBinError(i, 0);
          }
      }
  }
  else if (option.Contains("mc")) {
    auto wNum = hpass->GetSumw2();
    auto wDen = htotal->GetSumw2();
    for (int i = 1; i <= hpass->GetNbinsX(); ++i) {
      double pass = hpass->GetBinContent(i);
      double tot  = htotal->GetBinContent(i);
      pass = std::min(pass, tot);  // safe check
      if (tot > 0) {
        double eff = pass / tot;
        double err = (1. / tot )*sqrt(wNum->At(i)+wDen->At(i)*eff*eff);  
        ratio->SetBinContent(i, eff);
        ratio->SetBinError(i, err);
      } else {
        ratio->SetBinContent(i, 0);
        ratio->SetBinError(i, 0);
      }
    }
  }
  else if(option.Contains("bayesian")){
    for (int i=1; i<=hpass->GetNbinsX(); ++i){
      double pass = hpass->GetBinContent(i);
      double tot = htotal->GetBinContent(i);
      double eff = tot>0 ? pass/tot : 0;
      ratio->SetBinContent(i, eff);
      double err = eff>0 ? sqrt((pass+1)*(pass+2)/((tot+2)*(tot+3))-(pass+1)*(pass+1)/((tot+2)*(tot+2))) : 0;
      ratio->SetBinError(i, err);
    }
  }
  else if (option.Contains("d'Agostini")) {
    for (int i = 1; i <= hpass->GetNbinsX(); ++i) {
        double pass = hpass->GetBinContent(i);
        double tot = htotal->GetBinContent(i);
        if (tot > 0) {
            double eff = (pass + 1.0) / (tot + 2.0);
            double var = ((pass + 1.0) * (tot - pass + 1.0)) / ((tot + 2.0) * (tot + 2.0) * (tot + 3.0));
            double err = sqrt(var);
            ratio->SetBinContent(i, eff);
            ratio->SetBinError(i, err);
        } else {
            ratio->SetBinContent(i, 0);
            ratio->SetBinError(i, 0);
        }
    }
  }
  else if(option.Contains("ratio")){
    for (int i=1; i<=hpass->GetNbinsX(); ++i){
      double num = hpass->GetBinContent(i);
      double den = htotal->GetBinContent(i);
      double e_num = hpass->GetBinError(i);
      double e_den = htotal->GetBinError(i);
      double value = den>0 ? num/den : 0;
      ratio->SetBinContent(i, value);
      double err = (num!=0 && den!=0) ? value*sqrt( (e_num/num)*(e_num/num)+(e_den/den)*(e_den/den) ): 0;
      ratio->SetBinError(i, err);
    }
  }
  else if(option.Contains("simple")){
    ratio = (TH1D*)hpass->Clone();
    ratio->Divide(htotal);
  }
  else if (option.Contains("alberto")) {
    for (int i = 1; i <= hpass->GetNbinsX(); ++i) {
      double A  = hpass->GetBinContent(i);
      double T  = htotal->GetBinContent(i);
      double B  = T-A;
      double eA = hpass->GetBinError(i);
      double eT = htotal->GetBinError(i);
      double eB = sqrt( pow(eT,2)-pow(eA,2) );

      double value = (T > 0) ? A / T : 0;
      double err = 0.0;

      if (A > 0 && T > 0) {
        err = sqrt(pow(eB, 2)*A + pow(eA, 2)*B) / (T * T);
      }

      ratio->SetBinContent(i, value);
      ratio->SetBinError(i, err);
    }
  }
  else if (option.Contains("trigger")) {
    for (int i = 1; i <= hpass->GetNbinsX(); ++i) {
      double Np  = hpass->GetBinContent(i);
      double sNp  = hpass->GetBinError(i);
      double Nunb = htotal->GetBinContent(i);
      double sNunb = htotal->GetBinError(i);
      double denom = Np + Nunb;
      if (denom <= 0) {
        ratio->SetBinContent(i, 0);
        ratio->SetBinError(i, 0);
        continue;
      }
      double eps = Np / denom;
      double dNp   =  Nunb / (denom * denom);
      double dNunb = -Np   / (denom * denom);
      double s_eps = std::sqrt(dNp*dNp * sNp*sNp + dNunb*dNunb * sNunb*sNunb);

      ratio->SetBinContent(i, eps);
      ratio->SetBinError(i, s_eps);
    }
  }
  return ratio;
}
void RebinHistAbove(TH1D* h, double x) {
    if (!h) return;

    int nbins = h->GetNbinsX();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    // Determine where x falls
    int cutBin = h->FindBin(x);
    if (x >= xmax || cutBin > nbins) return; // nothing to rebin

    // Build new bin edges
    std::vector<double> edges;
    edges.reserve(cutBin + 2);
    for (int i = 1; i <= cutBin; ++i)
        edges.push_back(h->GetBinLowEdge(i));
    edges.push_back(x);     // boundary bin
    edges.push_back(xmax);  // final merged bin

    // Rebin using variable binning
    std::unique_ptr<TH1D> hnew((TH1D*)h->Rebin(edges.size() - 1,
                                               (std::string(h->GetName()) + "_tmp").c_str(),
                                               edges.data()));

    // Copy back content and errors
    h->SetBins(edges.size() - 1, edges.data());
    for (int i = 1; i <= hnew->GetNbinsX(); ++i) {
        h->SetBinContent(i, hnew->GetBinContent(i));
        h->SetBinError(i,   hnew->GetBinError(i));
    }
}
void RebinByEfficiencyOverall(TH1D* h, int Eff, unsigned int Charge) {
    using namespace SplineUtility;

    int rebinFactor = 4; // Default for most cases

    switch (static_cast<SplineUtility::Efficiency>(Eff)) {
        case TofEff:
            switch (Charge) {
                case 14: rebinFactor = 3; break;
                case 15:
                case 16: rebinFactor = 5; break;
                case 18: rebinFactor = 6; break;
                case 20:
                default: rebinFactor = 6; break;
            }
            break;

        case TrackEff:
            switch (Charge) {
                case 14: rebinFactor = 2; break;
                case 15:
                case 16: rebinFactor = 2; break;
                case 18: rebinFactor = 2; break;
                case 20: rebinFactor = 4; break;
                default: rebinFactor = 4; break;
            }
            break;
        case TriggerEff:
            switch (Charge) {
                case 14: rebinFactor = 2; break;
                case 15: rebinFactor = 2; break;
                case 16: rebinFactor = 2; break;
                case 18:
                case 20:
                default: rebinFactor = 6; break;
            }
            break;
        case L1Eff:
            switch (Charge) {
                case 14: rebinFactor = 4; break;
                case 15: rebinFactor = 4; break;
                case 16: rebinFactor = 2; break;
                case 18: 
                case 26: rebinFactor = 2; break;
                case 20: rebinFactor = 4; break;
                default: rebinFactor = 6; break;
            }
            break;
        case L1UnbEff:
            switch (Charge) {
                case 14: rebinFactor = 4; break;
                case 15:
                case 16: rebinFactor = 4; break;
                case 18: 
                case 20: rebinFactor = 4; break;
                default: rebinFactor = 2; break;
            }
            break;
        case TrackChEff:
            switch (Charge) {
                case 14: rebinFactor = 4; break;
                case 15:
                case 16: rebinFactor = 3; break;
                case 18: rebinFactor = 6; break;
                case 20: rebinFactor = 4; break;
                default: rebinFactor = 2; break;
            }
            break;
        case L9Eff:
            // Do nothing for L9, as specified
            return;
        case DaqEff:
            switch (Charge) {
                case 14:
                case 15: rebinFactor = 4; break;
                case 16: rebinFactor = 2; break;
                case 18: rebinFactor = 4; break;
                case 20: rebinFactor = 6; break;
                default: rebinFactor = 1; break;
            }
            break;
        default:
            rebinFactor = 1;
            break;
    }
    Rebin(h, rebinFactor);
}
void RebinByEfficiencyAbove(TH1D* h, int Eff, unsigned int Charge) {
    using namespace SplineUtility;

    double cut = -1; // valore della cut in R
    int rebinFactor = 1;

    // Definizione cut in base all’efficienza
    switch (static_cast<SplineUtility::Efficiency>(Eff)) {
        case L1Eff:
            cut = 100; break;
        case TofEff:
            switch (Charge) {
                case 14: cut = 200; break;
                case 15: cut = 20; break;
                case 16: cut = 30; break;
                default: cut = 30; break;
            }
        case TrackEff:
            cut = 30;
            break;
        case TriggerEff:
            switch (Charge) {
                case 14: cut = 1000; break;
                case 15: cut = 90; break;
                case 16: cut = 400; break;
                default: cut = 1000; break;
            }
            break;
        case L1UnbEff:
            switch (Charge) {
                case 14: cut = 100; break;
                case 15: cut = 30; break;
                case 16: cut = 100; break;
                default: cut = 100; break;
            }
        case TrackChEff:
            switch (Charge) {
                case 14: cut = 200; break;
                case 15: cut = 10; break;
                case 16: cut = 30; break;
                default: cut = 1000; break;
            }
            cut = 300;
            break;
        case DaqEff:
            switch (Charge) {
                case 14: cut = 600; break;
                case 15: cut = 10; break;
                case 16: cut = 100; break;
                default: cut = 1000; break;
            }
        default:
            cut = -1;
            break;
    }

    // --- OPZIONE E: Adaptive rebinning sopra la cut ---
    if (cut > 0) {
        const int nbins = h->GetNbinsX();
        std::vector<double> newEdges;
        newEdges.push_back(h->GetXaxis()->GetXmin());

        double sumEntries = 0;
        const double minEntries = 50;  // minimo eventi per bin sopra la cut
        const double minRelErr = 0.05; // errore relativo massimo (≈5%)

        for (int i = 1; i <= nbins; ++i) {
            double xlow = h->GetBinLowEdge(i);
            double xup = h->GetBinLowEdge(i + 1);

            // Prima della cut → bin invariati
            if (xup <= cut) {
                newEdges.push_back(xup);
                continue;
            }

            sumEntries += h->GetBinContent(i);

            double err = h->GetBinError(i);
            double relErr = (sumEntries > 0) ? err / sumEntries : 1.0;

            // Unisci finché non hai abbastanza statistica
            if (sumEntries >= minEntries || relErr < minRelErr) {
                newEdges.push_back(xup);
                sumEntries = 0;
            }
        }

        // Se l’ultimo bin non è chiuso
        if (newEdges.back() < h->GetXaxis()->GetXmax())
            newEdges.push_back(h->GetXaxis()->GetXmax());

        // Crea istogramma rebinnato
        TH1D* hRebinned = (TH1D*)h->Rebin(newEdges.size() - 1, 
                                          Form("%s_rebinned", h->GetName()), 
                                          newEdges.data());
        *h = *hRebinned;
        delete hRebinned;
    } 
    else {
        // Caso standard: rebinning fisso
        switch (static_cast<SplineUtility::Efficiency>(Eff)) {
            case TofEff:
                switch (Charge) {
                    case 14: rebinFactor = 3; break;
                    case 15:
                    case 16: rebinFactor = 5; break;
                    case 18:
                    case 20:
                    default: rebinFactor = 6; break;
                }
                break;
            case TrackEff:
                switch (Charge) {
                    case 14: rebinFactor = 2; break;
                    case 15:
                    case 16: rebinFactor = 2; break;
                    case 18: rebinFactor = 2; break;
                    case 20:
                    default: rebinFactor = 4; break;
                }
                break;
            case TriggerEff:
                switch (Charge) {
                    case 14: rebinFactor = 4; break;
                    case 15:
                    case 16: rebinFactor = 5; break;
                    case 18:
                    case 20:
                    default: rebinFactor = 6; break;
                }
                break;
            case L1Eff:
                switch (Charge) {
                    case 14: rebinFactor = 2; break;
                    case 15: rebinFactor = 4; break;
                    case 16: rebinFactor = 2; break;
                    case 18: 
                    case 26: rebinFactor = 2; break;
                    case 20: rebinFactor = 4; break;
                    default: rebinFactor = 6; break;
                }
                break;
            case L1UnbEff:
                switch (Charge) {
                    case 14: rebinFactor = 4; break;
                    case 15:
                    case 16: rebinFactor = 4; break;
                    case 18:
                    case 20: rebinFactor = 4; break;
                    default: rebinFactor = 2; break;
                }
                break;
            case TrackChEff:
                switch (Charge) {
                    case 14: rebinFactor = 4; break;
                    case 15:
                    case 16: rebinFactor = 3; break;
                    case 18: rebinFactor = 6; break;
                    case 20: rebinFactor = 4; break;
                    default: rebinFactor = 2; break;
                }
                break;
            case DaqEff:
                switch (Charge) {
                    case 14:
                    case 15: rebinFactor = 4; break;
                    case 16: rebinFactor = 2; break;
                    case 18: rebinFactor = 4; break;
                    case 20: rebinFactor = 6; break;
                    default: rebinFactor = 1; break;
                }
                break;
            default:
                rebinFactor = 1;
                break;
        }

        Rebin(h, rebinFactor);
    }
}

void Rebin(TH1* h, int n) {
  h->Rebin(n);
  //h->Scale(1.0/n);
}
TH2D* GetHist2D(TFile* file, const TString& geom, const TString& name) {
    if (!file) return nullptr;
    TString path = geom + "/" + name;
    auto h = dynamic_cast<TH2D*>(file->Get(path));
    if (!h) {
        std::cerr << "[GetHist2D] Warning: histogram " << path << " not found in "
                  << file->GetName() << std::endl;
    }
    return h;
}
TH1D *emulateTrigger(unsigned int charge, int start_bin, int stop_bin, TString geom, TString rebin) {
  TString name = getIonName(charge);
  TString fileName="../IonsSelected/"+getIonPath(charge)+"/Efficiencies/Trigger/emu.root";
  TFile *file = TFile::Open(fileName.Data());
  if (!file || file->IsZombie()) 
    printf("Errore nell'aprire il file \n");
  auto nPhys1   = (TH1D*) GetHist2D(file, geom, "nPhys_1")
                      ->ProjectionY("p1", hist_time->FindBin(0),hist_time->FindBin(1456503197));

  auto nPhys2   = (TH1D*) GetHist2D(file, geom, "nPhys_2")
                      ->ProjectionY("p2", hist_time->FindBin(1456503197), stop_bin);

  auto acc5_pre = (TH1D*) GetHist2D(file, geom, "acc5_counters_pre")
                      ->ProjectionY("acc5_pre", hist_time->FindBin(0), hist_time->FindBin(1456503197));

  auto acc8     = (TH1D*) GetHist2D(file, geom, "acc8_counters")
                      ->ProjectionY("acc8", hist_time->FindBin(1456503197), stop_bin);

  auto acc5     = (TH1D*) GetHist2D(file, geom, "acc5_counters")
                      ->ProjectionY("acc5", hist_time->FindBin(1456503197), stop_bin);

  if (rebin=="rebin") {
    RebinByEfficiencyOverall(nPhys1,SplineUtility::Efficiency::TriggerEff,charge);
    RebinByEfficiencyOverall(nPhys2,SplineUtility::Efficiency::TriggerEff,charge);
    RebinByEfficiencyOverall(acc5_pre ,SplineUtility::Efficiency::TriggerEff,charge);
    RebinByEfficiencyOverall(acc8 ,SplineUtility::Efficiency::TriggerEff,charge);
    RebinByEfficiencyOverall(acc5 ,SplineUtility::Efficiency::TriggerEff,charge);
  }
  
  int nBins = acc5->GetNbinsX();
  TH1D *eff = (TH1D*)acc5->Clone("eff_trig");
  eff->Reset();

  for (int bin = 1; bin <= nBins; ++bin) {
      double Np1 = nPhys1->GetBinContent(bin);
      double Np2 = nPhys2->GetBinContent(bin);
      double Nacc = acc5->GetBinContent(bin);

      double eNp1 = nPhys1->GetBinError(bin);
      double eNp2 = nPhys2->GetBinError(bin);
      double eNacc = acc5->GetBinError(bin);

      // A and var(A)
      double A = Np1 + Np2;
      double varA = eNp1*eNp1 + eNp2*eNp2;

      // B and var(B)
      double B = 0.0;
      double varB = 0.0;
      double covA_B = 0.0;

      if (Nacc > 0.0) {
          B = Np1 * Np2 / Nacc + Np2;

          double dB_dNp1 = Np2 / Nacc;
          double dB_dNp2 = 1.0 + Np1 / Nacc;
          double dB_dNacc = - (Np1 * Np2) / (Nacc * Nacc);

          varB = dB_dNp1*dB_dNp1 * eNp1*eNp1
              + dB_dNp2*dB_dNp2 * eNp2*eNp2
              + dB_dNacc*dB_dNacc * eNacc*eNacc;

          // Cov(A,B) ≈ dB/dNp1 Var(Np1) + dB/dNp2 Var(Np2)
          covA_B = dB_dNp1 * (eNp1*eNp1) + dB_dNp2 * (eNp2*eNp2);
      }

      double eff_val = 0.0;
      double eff_err = 0.0;

      if (B > 0.0) {
          eff_val = A / B;

          // Var(R) = (1/B^2) Var(A) + (A^2/B^4) Var(B) - 2A/B^3 Cov(A,B)
          double varR = (1.0 / (B*B)) * varA
                      + (A*A) / (B*B*B*B) * varB
                      - 2.0 * A / (B*B*B) * covA_B;

          // numerical safety
          if (!std::isfinite(varR) || varR < 0.0) varR = 0.0;
          eff_err = std::sqrt(varR);

          // optional physical clamps
          if (eff_val < 0.0) eff_val = 0.0;
          if (eff_val > 1.0) eff_val = 1.0;
          if (eff_err > 1.0) eff_err = 1.0; // conservative cap
          if (eff_val + eff_err > 1.0) eff_err = 1.0 - eff_val;
          if (eff_val - eff_err < 0.0) eff_err = eff_val;
      }

      eff->SetBinContent(bin, eff_val);
      eff->SetBinError(bin, eff_err);
  }


  return eff;
}
TH1D* emulateDAQ(unsigned int charge, int start_bin, int stop_bin, TString geom, TString rebin) {
    TString name = getIonName(charge);
    TString fileName = "../IonsSelected/" + getIonPath(charge) + "/Efficiencies/Daq/daq.root";
    TFile *f = TFile::Open(fileName.Data());
    if (!f || f->IsZombie()) {
        printf("Errore nell'aprire il file \n");
        return nullptr;
    }

    // Get 2D histos
    auto j1_above = (TH1D*) GetHist2D(f, geom, "jinj1_aboveBuffer")
                        ->ProjectionY("j1_above", start_bin, stop_bin);

    auto j1_below = (TH1D*) GetHist2D(f, geom, "jinj1_belowBuffer")
                        ->ProjectionY("j1_below", start_bin, stop_bin);

    auto j2_above = (TH1D*) GetHist2D(f, geom, "jinj2_aboveBuffer")
                        ->ProjectionY("j2_above", start_bin, stop_bin);

    auto j2_below = (TH1D*) GetHist2D(f, geom, "jinj2_belowBuffer")
                        ->ProjectionY("j2_below", start_bin, stop_bin);

    // Numerator and denominator in raw counts
    auto num_counts = (TH1D*)j1_above->Clone("num_counts");
    auto den_counts = (TH1D*)j1_above->Clone("den_counts");

    for (int i = 1; i <= num_counts->GetNbinsX(); i++) {
        double Npass = j1_below->GetBinContent(i) +
                       j2_below->GetBinContent(i) +
                       j2_above->GetBinContent(i);

        double Ntotal = j1_below->GetBinContent(i) *
                        (1.0 + (j2_above->GetBinContent(i) / (j2_below->GetBinContent(i) > 0 ? j2_below->GetBinContent(i) : 1.0))) +
                        j2_below->GetBinContent(i) +
                        j2_above->GetBinContent(i);

        num_counts->SetBinContent(i, Npass);
        num_counts->SetBinError(i, std::sqrt(Npass)); // Poisson error

        den_counts->SetBinContent(i, Ntotal);
        den_counts->SetBinError(i, std::sqrt(Ntotal)); // Poisson error
    }

    // Rebin counts first
    if (rebin == "rebin") {
      RebinByEfficiencyOverall(num_counts,SplineUtility::Efficiency::DaqEff,charge);
      RebinByEfficiencyOverall(den_counts,SplineUtility::Efficiency::DaqEff,charge);
    }
    // Create efficiency histogram with proper binomial errors
    auto eff = (TH1D*)num_counts->Clone("daq_efficiency");
    eff->Divide(num_counts, den_counts, 1.0, 1.0, "B"); // "B" = binomial errors

    return eff;
}
TH1D* getL1chEff(unsigned int charge) {
    TString name = getIonName(charge);
    TString fileName = "../IonsSelected/" + getIonPath(charge) + "/Efficiencies/L1ChargeCut/l1ChargeCut.root";
    TFile *f = TFile::Open(fileName.Data());
    auto h = (TH1D*)f->Get("hL1Eff");
    if (!h) return nullptr;
    return h;
}
double findMaxRelDiff(TH1D *flux, TH1D *oldflux, int prod_min, int prod_max) {
  double max_rel_diff = 0;
  for(int i=prod_min; i<=flux->FindBin(prod_max); ++i) {
    double old_value = oldflux->GetBinContent(i);
    if (old_value != 0) {
      auto diff = round(10000 * fabs(flux->GetBinContent(i) - old_value) / old_value) / 100;
        if (diff > max_rel_diff) max_rel_diff = diff;
    }
 }
 return max_rel_diff;
}
TH1D *relDiff(TH1D *flux, TH1D *oldflux, int prod_min, int prod_max){
  TH1D* h = (TH1D*)hist_rig_highZ->Clone("h");
  for (int i=0; i<h->GetNbinsX(); i++) {
    double old_value = oldflux->GetBinContent(i);
    if (old_value != 0) {
      auto diff = round(10000 * fabs(flux->GetBinContent(i) - old_value) / old_value) / 100;
      h->SetBinContent(i,diff);
    }
  }
  return h;
}
std::pair<int, int> getTimePeriod(TString time) {
  if (time=="10y") {
    std::cout << "Choosen time period : 10y" << std::endl;
    return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1620259200) ); //19 May 2011 to 6 May 2021 (10y)
  }
  if (time=="11y") {
    std::cout << "Choosen time period : 11y" << std::endl;
    return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1654087477) ); //19 May 2011 to 1 June 2022 (11y)
  }
  if (time=="12.5y") {
    std::cout << "Choosen time period : 12.5y" << std::endl;
    return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1699574400) ); //19 May 2011 to 10 November 2023 (12.5y)
  }
  if (time=="13.5y") {
    std::cout << "Choosen time period : 13.5y" << std::endl;
    return std::make_pair(hist_time->FindBin(1305763200),hist_time->FindBin(1731798000) ); //19 May 2011 to 16 November 2024 (13.5y)
  }
  else return std::make_pair(0,0);
}
int getNucleiCount(unsigned int charge, TString timePeriod) {
  TString name = getIonName(charge);
  TString ionPath = getIonPath(charge);
  TString fileName="../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+Form("_UnfoldingFactor_%s.root",timePeriod.Data());
  TFile *file = TFile::Open(fileName.Data());
  if (!file || file->IsZombie()) 
    printf("Errore nell'aprire il file \n");
  auto counts = (TH1D *)file->Get("countsRaw");
  return counts->Integral();
}
double GetDataMass(unsigned int charge) {
  double A = 0;
  switch (charge) {
  case 1:
    A = 1;
    break;
  case 2:
    A = 4;
    break;
  case 3:
    A = 6.5;
    break;
  case 4:
    A = 7.6;
    break;
  case 5:
    A = 10.65;
    break;
  case 6:
    A = 12;
  case 7:
    A = 14;
    break;
  case 8:
    A = 16;
    break;
  case 9:
    A = 18;
    break;
  case 10:
    A = 20;
    break;
  case 11:
    A = 22;
    break;
  case 12:
    A = 24;
    break;
  case 13:
    A = 26;
    break;
  case 14:
    A = 28;
    break;
  case 15:
    A = 31;
	break;
  case 16:
    A = 32;
    break;
  case 17:
    A = 34;
    break;
  case 18:
    A = 36;
	  break;
  case 19:
    A = 38;
    break;
  case 20:
    A = 40;
    break;
  case 21:
    A = 42;
    break;
  case 22:
    A = 44;
    break;
  case 23:
    A = 46;
    break;
  case 24:
    A = 48;
    break;
  case 25:
    A = 50;
    break;
  case 26:
    A = 52;
    break;
  }
  return A;
}
void printDistributions(std::vector<TH1D*> dat,
           std::vector<TH1D*> mc,
           std::vector<unsigned int> charge_number,
           const TString &tag,
           const std::vector<TString> &histogram,
           const std::map<TString,TString> &name) {
    TStyle *style = effHistStyle(0);
    gROOT->SetStyle("MyStyle");
    gROOT->ForceStyle();
    setLogon();
    style->cd();

    TString prefix = getPrefixFile(charge_number);
    TCanvas *b = new TCanvas("","",2048,1280);
    b->SaveAs(Form("../output/%s_%s.pdf[", prefix.Data(), tag.Data()), "RECREATE");

    int nhist   = histogram.size();
    int ncharge = charge_number.size();
    auto a = new TLegend(0.67,0.67,0.9,0.9);
    gPad->SetLogy();

    for (int i=0; i<ncharge; i++) {
        auto qRange = getQRange(charge_number[i]);
        b->Update();
        TString title = getIonName(i);

        TH2D *h1 = new TH2D("", "", 550, 0, 28, 500, 0.1, 10e+6);
        setTitle(h1, "d", Form("Q_{%s}", tag.Data()), "Entries", charge_number[i]);
        formatAxis(h1, 1);
        formatTitle(h1,1);
        adjustZoomX(h1,tag);
        h1->GetYaxis()->SetLabelFont(62);
        h1->GetXaxis()->SetLabelFont(62);
        h1->GetXaxis()->SetRangeUser(qRange.first,qRange.second);

        TString titleHist = h1->GetTitle();
        TString t1 = Form("%s Data",titleHist.Data());
        h1->SetTitle(t1);
        getColor(h1, charge_number[i], "dat");
        h1->Draw();

        //-- Drawing data (check for nulls and bounds)
        for (int j=0; j<nhist; j++) {
            size_t idx = j + i * nhist;
            if (idx < dat.size() && dat[idx]) {
                a->AddEntry(dat[idx], name.at(histogram[j]));
                setFriendColor(dat[idx], j);
                dat[idx]->SetMarkerStyle(20);
                dat[idx]->SetMarkerSize(1.6);
                dat[idx]->Draw("SAME");
            } else {
                printf("Warning: Data histogram missing for charge idx %d, hist '%s' (vec idx %zu)\n",
                       i, histogram[j].Data(), idx);
            }
        }

        a->Draw("SAME");
        a->SetTextSize(0.028);
        b->SaveAs(Form("../output/%s_%s.pdf", prefix.Data(), tag.Data()), "RECREATE");

        //-- Setup for MC page title & background
        TString t2 = Form("%s Mc",titleHist.Data());
        h1->SetTitle(t2);
        h1->Draw();
        getColor(h1, charge_number[i], "mc");
        a->Clear();

        //-- Drawing mc (only if available)
        for (int j=0; j<nhist; j++) {
            size_t idx = j + i * nhist;
            if (idx < mc.size() && mc[idx]) {
                a->AddEntry(mc[idx], name.at(histogram[j]));
                setFriendColor(mc[idx], j);
                mc[idx]->SetMarkerStyle(20);
                mc[idx]->SetMarkerSize(1.6);
                mc[idx]->Draw("SAME");
            } else {
                if (mc.empty()) {
                    if (i==0 && j==0)
                        printf("Info: mc vector empty; skipping mc plots.\n");
                } else {
                    printf("Warning: MC histogram missing for charge idx %d, hist '%s' (vec idx %zu)\n",
                           i, histogram[j].Data(), idx);
                }
            }
        }

        a->Draw("SAME");
        a->SetTextSize(0.028);
        b->SaveAs(Form("../output/%s_%s.pdf", prefix.Data(), tag.Data()), "RECREATE");
        a->Clear();
    }

    b->SaveAs(Form("../output/%s_%s.pdf]", prefix.Data(), tag.Data()), "RECREATE");
}
void printEff(std::vector<TH2D*> dat,
              std::vector<TH2D*> mc,
              std::vector<unsigned int> charge_number,
              const TString &tag,
              const std::vector<TString> &histogram,
              const std::map<TString,TString> &name) {
    TStyle *style = effHistStyle(0);
    gROOT->SetStyle("MyStyle");
    gROOT->ForceStyle();
    setLogon();
    style->cd();

    TString prefix = getPrefixFile(charge_number);
    int nhist   = histogram.size();
    int ncharge = charge_number.size();
    auto a = formatLegend(nhist);

    for (int i=0; i<ncharge; i++) {
        TCanvas *b = new TCanvas("","",2048,1280);
        if(i==0) b->SaveAs(Form("../output/%s_%s_eff.pdf[", prefix.Data(), tag.Data()),"RECREATE");

        b->Update();
        b->Divide(nhist,1,0,0);
        TString effName = Form("%s efficiency", tag.Data());

        TH2D *h = new TH2D("", "", nRbins_HighZ - 1, Rbins_HighZ, 200, 0.2, 1.15);
        setTitle(h, effName, "R (GV)", effName+" #bf{#varepsilon}", charge_number[i]);
        formatAxis(h, nhist);
        formatTitle(h,nhist);
        adjustZoomY(h,effName,nhist);
        adjustZoomX(h,effName);
        TString firstTitle = h->GetTitle();

        for (int j=0; j<nhist; j++) {
            b->Update();
            b->cd(j+1);
            gPad->SetLogx();
            a->Clear();

            TString t1 = Form("%s %s",firstTitle.Data(),name.at(histogram[j]).Data());
            h->SetTitle(t1);

            int idx = 2*j + i*nhist*2;
            TH1D *d_samp = nullptr, *d_pass = nullptr, *m_samp = nullptr, *m_pass = nullptr;
            if ((size_t)idx < dat.size()) d_samp = dat[idx] ? dat[idx]->ProjectionY(Form("a%u", idx)) : nullptr;
            if ((size_t)(idx+1) < dat.size()) d_pass = dat[idx+1] ? dat[idx+1]->ProjectionY(Form("b%u", idx)) : nullptr;
            if ((size_t)idx < mc.size()) m_samp = mc[idx] ? mc[idx]->ProjectionY(Form("c%u", idx)) : nullptr;
            if ((size_t)(idx+1) < mc.size()) m_pass = mc[idx+1] ? mc[idx+1]->ProjectionY(Form("d%u", idx)) : nullptr;

            if (!d_samp || !d_pass) {
                printf("Warning: missing data hist for eff (charge idx %d, hist %s). Skipping.\n",
                       i, histogram[j].Data());
                continue;
            }

            auto d = divide(d_pass,d_samp,"efficiency");
            d->Draw("P E1 SAME");
            a->AddEntry(d,"Data");
            getColor(d, charge_number[i], "dat");

            if (m_samp && m_pass) {
                auto m = divide(m_pass,m_samp,"efficiency");
                m->Draw("P E1 SAME");
                a->AddEntry(m,"Mc");
                getColor(m, charge_number[i], "mc");
            }

            if (j==nhist-1) a->Draw();
        }

        b->SaveAs(Form("../output/%s_%s_eff.pdf", prefix.Data(), tag.Data()),"RECREATE");
        if(i==ncharge-1) b->SaveAs(Form("../output/%s_%s_eff.pdf]", prefix.Data(), tag.Data()),"RECREATE");
    }
}

#endif
