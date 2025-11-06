#ifndef FLUXUTILITY_H
#define FLUXUTILITY_H
#include "includes.h"
#include "SplineUtility.h"

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
    A = "Lithium";
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
  case 11:
    A = "Sodium";
    break;
  case 12:
    A = "Magnesium";
    break;
  case 13:
    A = "Aluminum";
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
  case 18:
    A = "Argon";
    break;
  case 19:
    A = "Potassium";
    break;
  case 20:
    A = "Calcium";
    break;
  case 21:
    A = "Scandium";
    break;
  case 22:
    A = "Titanium";
    break;
  case 23:
    A = "Vanadium";
    break;
  case 24:
    A = "Chromium";
    break;
  case 25:
    A = "Manganese";
    break;
  case 26:
    A = "Iron";
    break;
  case 27:
    A = "Cobalt";
    break;
  case 28:
    A = "Nickel";
    break;
  }
  return A;
}
void ApplyPowerTransformation(std::vector<double>& x, std::vector<double>& y, 
							  std::vector<double>& exl, std::vector<double>& exh, double power) {
    for (size_t i = 0; i < y.size(); ++i) {
        y[i] *= std::pow(x[i], power);
		exl[i] *= std::pow(x[i], power);
		exh[i] *= std::pow(x[i], power);
    }
}
void MultiplyByXPower(TH1 *hh, double power) {
  for (unsigned int ibin = 0; ibin < hh->GetNbinsX(); ibin++) {
    hh->SetBinContent(ibin + 1, hh->GetBinContent(ibin + 1) * pow(hh->GetBinCenter(ibin + 1), power));
    hh->SetBinError(ibin + 1, hh->GetBinError(ibin + 1) * pow(hh->GetBinCenter(ibin + 1), power));
  }
}
TF1* MultiplyTF1ByXPower(TF1* f, double power, const char* newName = "f_modified") {
    int npar = f->GetNpar();

    auto wrapped = new TF1(newName, [=](double* x, double* p) {
        // Set parameters of original function to current p
        std::vector<double> parVec(p, p + npar);
        TF1 fcopy(*f);
        for (int i = 0; i < npar; ++i) fcopy.SetParameter(i, parVec[i]);

        return pow(x[0], power) * fcopy.Eval(x[0]);
    }, f->GetXmin(), f->GetXmax(), npar);

    // Set parameter names and initial values
    for (int i = 0; i < npar; ++i) {
        wrapped->SetParameter(i, f->GetParameter(i));
        wrapped->SetParName(i, f->GetParName(i));
    }

    return wrapped;
}
TF1* BuildRecoveredFlux(TF1* fluxModel, double Rmin, double Rmax, const std::string& name = "recoveredFlux") { //This one multiplies by R^(-2.7)
    if (!fluxModel) return nullptr;
    TF1* recoveredFlux = new TF1(name.c_str(),
        [fluxModel](double* x, double*) -> double {
            double R = x[0];
            if (R <= 0) return 0.0;
            return fluxModel->Eval(R) / std::pow(R, 2.7);
        },
        Rmin, Rmax, 0);
    return recoveredFlux;
}
std::string getCSVfilename(unsigned int charge) {
  std::string s ="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Published/";
  switch (charge) { 
    case 8:
      s+="table-sm-ii-O.csv";
      break;
    case 9:
      s+="table-siii.csv";
      break;
    case 14:
      s+="table-sm-v-Si.csv";
      break;
    case 15:
      s+="P_Sdiat.csv";
      break;
    case 16:
      s+="table-sm-i-S.csv";
      break;
  }
  return s;
}
void readFluxTable(std::vector<double> &lower_bounds, std::vector<double> &upper_bounds, std::vector<double> &values, std::vector<double> &errors,
                  unsigned int charge, std::vector<double> &syst) {
    // Get the CSV table file
    std::string filename = getCSVfilename(charge);
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
    }
    // Read the CSV file
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token = "";
        double lower_bound, upper_bound, value;
        double err1, err2, err3, err4, err5;
        // Read the lower bound
        std::getline(ss, token, ',');
        lower_bound = std::stod(token);
        // Read the upper bound
        std::getline(ss, token, ',');
        upper_bound = std::stod(token);
        // Read the value
        std::getline(ss, token, ',');
        value = std::stod(token);
        // Read the 5 errors
        std::getline(ss, token, ',');
        err1 = std::stod(token);
        std::getline(ss, token, ',');
        err2 = std::stod(token);
        std::getline(ss, token, ',');
        err3 = std::stod(token);
        std::getline(ss, token, ',');
        err4 = std::stod(token);
        std::getline(ss, token, ',');
        err5 = std::stod(token);
        // Compute the combined error (assuming they are statistical errors that add in quadrature)
        double combined_error = std::sqrt(err1*err1+err5*err5);
        // Store the data
        lower_bounds.push_back(lower_bound);
        upper_bounds.push_back(upper_bound);
        values.push_back(value);
        syst.push_back(err5);
        errors.push_back(combined_error);
    }
    file.close();

}
TH1D* buildPublished(std::vector<double> lower_bounds, std::vector<double> upper_bounds, std::vector<double> values, std::vector<double> errors) {
    int n = values.size();
    // Creiamo un vettore con i bordi dei bin
    std::vector<double> bin_edges(n + 1);
    for (int i = 0; i < n; ++i) {
        bin_edges[i] = lower_bounds[i];
    }
    bin_edges[n] = upper_bounds[n - 1];  // L'ultimo bordo superiore è fuori dal loop

    // Creiamo l'istogramma con i bin definiti
    TH1D* h = new TH1D("histogram", "Histogram from buildPublished", n, bin_edges.data());

    // Riempimento dei bin dell'istogramma
    for (int i = 0; i < n; ++i) {
        double bin_center = (lower_bounds[i] + upper_bounds[i]) / 2.0;
        h->SetBinContent(i + 1, values[i]);
        h->SetBinError(i + 1, errors[i]);
    }
    MultiplyByXPower(h,2.7);
    return h;
}
TH1D* GetPublishedFlux(unsigned int charge){
    if (charge > 16) {
      auto empty = (TH1D*)hist_rig_highZ->Clone("empty");
      return empty;
    }
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> values;
    std::vector<double> syst;
    std::vector<double> errors;
    readFluxTable(lower_bounds,upper_bounds,values,errors,charge,syst);
    return buildPublished(lower_bounds,upper_bounds,values,errors);
}
TH1D* GetQYanFlux(unsigned int charge){
  TString s ="../Published/QYan/hzbin/";
  TString h = Form("Z%ufluxh_totalA",charge);
  if (charge == 9) {
    s+="fluorine69_20250304P8_B1236602RAMCKY13COMBUNFOLDTOINNB_totalQYAN.root";
  }
  if (charge == 14) {
    s+="silicon69_20250304P8_B1236602RAMCKY13COMBUNFOLDNNB_totalQYAN.root";
  }
  if (charge == 15) {
    s+="phosphorus69_20250304P8_B1236602RAMCKY13COMBUNFOLDTOINNB_totalQYAN.root";
  }
  if (charge == 16) {
    s+="sulfur69_20250304P8_B1236602RAMCKY13COMBUNFOLDTOINNB_totalQYAN.root";
  }
  if (charge == 17) {
    s+="chlorine69_20250304P8_B1236602RAMCKY13COMBUNFOLDTOINNB_totalQYAN.root";
  }
  if (charge == 18) {
    s+="argon69_20250304P8_B1236602RAMCKY13COMBUNFOLDTOINNB_totalQYAN.root";
  }
  if (charge == 19) {
    s+="potassium69_20250304P8_B1236602RAMCKY13COMBUNFOLDTOINNB_totalQYAN.root";
  }
  if (charge == 20) {
    s+="calcium69_20250304P8_B1236602RAMCKY13COMBUNFOLDTOINNB_totalQYAN.root";
  }
  if (charge == 21) {
    s+="ferrum68_20230304P8_B1236602RAMCKY11COMBUNFOLDNNB_totalQYAN.root";
    h = Form("Z%ufluxh_totalA",26);
  }
  if (charge == 22) {
    s+="ferrum68_20230304P8_B1236602RAMCKY11COMBUNFOLDNNB_totalQYAN.root";
    h = Form("Z%ufluxh_totalA",26);
  }
  if (charge == 23) {
    s+="ferrum68_20230304P8_B1236602RAMCKY11COMBUNFOLDNNB_totalQYAN.root";
    h = Form("Z%ufluxh_totalA",26);
  }
  if (charge == 24) {
    s+="ferrum68_20230304P8_B1236602RAMCKY11COMBUNFOLDNNB_totalQYAN.root";
    h = Form("Z%ufluxh_totalA",26);
  }
  if (charge == 25) {
    s+="ferrum68_20230304P8_B1236602RAMCKY11COMBUNFOLDNNB_totalQYAN.root";
    h = Form("Z%ufluxh_totalA",26);
  }
  if (charge == 26) {
    s+="ferrum68_20230304P8_B1236602RAMCKY11COMBUNFOLDNNB_totalQYAN.root";
  }
  std::cout << s << std::endl;
  auto file =  new TFile(s.Data());
  TH1D *flux = (TH1D*)file->Get(h.Data());
  MultiplyByXPower(flux,2.7);
  return flux;
}
TH1D* GetJoseFlux(unsigned int charge) {
  TString s ="../Published/Jose/toJack.root";
  TString h = "Z15_Flux_11yr_L1Inn";
  auto file =  new TFile(s.Data());
  TH1D *flux = (TH1D*)file->Get(h.Data());
  MultiplyByXPower(flux,2.7);
  return flux;
}
TH1D* GetFlux(unsigned int charge, TString timePeriod) {  
  TString name = getIonName(charge);
  TString fileName="../IonsSelected/"+getIonPath(charge)+"/Flux/"+getIonName(charge)+Form("_Flux_%s.root",timePeriod.Data());
  std::cout << fileName << std::endl;
  TFile *file = TFile::Open(fileName.Data());
  if (!file || file->IsZombie()) 
    printf("Errore nell'aprire il file \n");
  else 
    std::cout << "Flux from : " << fileName << std::endl;
  TH1D *flux= (TH1D*)file->Get("flux"); 
  MultiplyByXPower(flux,2.7);
  return flux;
}
TH1D* getJackFlux(unsigned int charge) {
  TString s ="../Published/Jack/mit_flux.root";
  TString h = "phosphorus_flux";
  auto file =  new TFile(s.Data());
  TH1D *flux = (TH1D*)file->Get(h.Data());
  MultiplyByXPower(flux,2.7);
  return flux;
}
TH1D* GetZhenFlux(unsigned int charge) {
  TString s ="/storage/gpfs_ams/ams/users/zliu/public/share/flux_zhen/flux";
  TString f ="/";
  switch(charge) {
    case 14:
      f+="Silicon_unfolded_BinS.root";
      break;
    case 15:
      f+="Phosphorus_unfolded_BinS.root";
      break;
    case 16:
      f+="Sulfur_unfolded_BinS.root";
      break;
    case 18:
      f+="Argon_unfolded.root";
      break;
    case 20:
      f+="Calcium_unfolded.root";
      break;
  }
  TFile* file = TFile::Open(s+f);
    TTree* tree = (TTree*)file->Get("unfolded_flux");

    double x, y, y_err;

    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("y_err", &y_err);

    Long64_t nentries = tree->GetEntries();

    std::vector<double> vx(nentries), vy(nentries), vey(nentries);
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        vx[i] = x;
        vy[i] = y;
        vey[i] = y_err;
    }

    // Create histogram with custom binning
    TH1D* h = (TH1D*)hist_rig_highZ->Clone("h");

    for (int i = 0; i < nentries; ++i) {
        int bin = h->FindBin(vx[i]);
        h->SetBinContent(bin, vy[i]);
        h->SetBinError(bin, vey[i]);
    }

    return h;

}
TH1D* DivideFluxes(TH1D *published, TH1D *my) {
  TH1D *ratio;
  for (int ii=0; ii<published->GetNbinsX(); ii++) {
    ratio->SetBinContent(ii, published->GetBinContent(ii)/my->GetBinContent(ii) );
  }
  return ratio;
}
TH1D *FluxOverAverage(const TH1D *flux, std::vector<double> average, TString option ) {
  auto h = (TH1D*)flux->Clone();
  for (int ii=0; ii<hist_rig_highZ->GetNbinsX(); ii++) {
    if (average[ii]!=0 ) {
      h->SetBinContent(ii, flux->GetBinContent(ii)/average[ii]);
      if (option == "mine")
        h->SetBinError(ii, flux->GetBinError(ii)/average[ii]);
      else 
        h->SetBinError(ii, flux->GetBinError(ii)/average[ii]);
    }
  }
  return h;
}
TH1D *FluxOverErrors(const TH1D *mine, const TH1D *pub) {
  auto h = (TH1D*)hist_rig_highZ->Clone();
  for (int ii=0; ii<h->GetNbinsX(); ii++) {
    h->SetBinContent(ii, (pub->GetBinContent(ii)-mine->GetBinContent(ii))/mine->GetBinError(ii)  );
    h->SetBinError(ii, 0.00001);
  }
  return h;
}
TH1D *ErrorContour(const TH1D *pub, std::vector<double> average) {
  auto h = (TH1D*)hist_rig_highZ->Clone();
  for (int ii=0; ii<h->GetNbinsX(); ii++) {
    h->SetBinContent(ii,1);
    h->SetBinError(ii, pub->GetBinError(ii)/average[ii] );
  }
  return h;
}
TH1D* RebinHistogramWeighted(TH1D* h_orig) {
    // Retrieve the bin edges of the new histogram
    int nNewBins = hist_rig_highZ->GetNbinsX();
    std::vector<double> newBins;
    for (int i = 1; i <= nNewBins + 1; ++i) {
        newBins.push_back(hist_rig_highZ->GetBinLowEdge(i));
    }

    // Create the new histogram with the desired binning
    TH1D* h_rebinned = new TH1D("", h_orig->GetTitle(), nNewBins, newBins.data());

    // Loop over the new bins
    for (int i = 1; i <= nNewBins; ++i) {
        double binLow = newBins[i - 1];
        double binHigh = newBins[i];
        
        double weightedSum = 0;
        double weightSum = 0;
        double errorSum = 0;

        // Loop over old bins and accumulate weighted sum
        for (int j = 1; j <= h_orig->GetNbinsX(); ++j) {
            double origBinCenter = h_orig->GetBinCenter(j);
            double origBinContent = h_orig->GetBinContent(j);
            double origBinError = h_orig->GetBinError(j);
            double origBinWidth = h_orig->GetBinWidth(j);

            if (origBinCenter >= binLow && origBinCenter < binHigh) {
                double weight = origBinWidth;  // Use bin width as weight
                weightedSum += origBinContent * weight;
                weightSum += weight;

                // Propagate error (weighted average error formula)
                if (origBinError > 0) {
                    errorSum += (origBinError * origBinError) * (weight * weight);
                }
            }
        }

        // Assign the weighted average to the new bin
        if (weightSum > 0) {
            h_rebinned->SetBinContent(i, weightedSum / weightSum);
            h_rebinned->SetBinError(i, sqrt(errorSum) / weightSum);
        } else {
            h_rebinned->SetBinContent(i, 0);
            h_rebinned->SetBinError(i, 0);
        }
    }

    return h_rebinned;
}
TH1D* RebinHistogramWeightedPub(TH1D* h_orig) {
    // Retrieve the bin edges of the new histogram
    int nNewBins = hist_rig_pub->GetNbinsX();
    std::vector<double> newBins;
    for (int i = 1; i <= nNewBins + 1; ++i) {
        newBins.push_back(hist_rig_pub->GetBinLowEdge(i));
    }

    // Create the new histogram with the desired binning
    TH1D* h_rebinned = new TH1D("", h_orig->GetTitle(), nNewBins, newBins.data());

    // Loop over the new bins
    for (int i = 1; i <= nNewBins; ++i) {
        double binLow = newBins[i - 1];
        double binHigh = newBins[i];
        
        double weightedSum = 0;
        double weightSum = 0;
        double errorSum = 0;

        // Loop over old bins and accumulate weighted sum
        for (int j = 1; j <= h_orig->GetNbinsX(); ++j) {
            double origBinCenter = h_orig->GetBinCenter(j);
            double origBinContent = h_orig->GetBinContent(j);
            double origBinError = h_orig->GetBinError(j);
            double origBinWidth = h_orig->GetBinWidth(j);

            if (origBinCenter >= binLow && origBinCenter < binHigh) {
                double weight = origBinWidth;  // Use bin width as weight
                weightedSum += origBinContent * weight;
                weightSum += weight;

                // Propagate error (weighted average error formula)
                if (origBinError > 0) {
                    errorSum += (origBinError * origBinError) * (weight * weight);
                }
            }
        }

        // Assign the weighted average to the new bin
        if (weightSum > 0) {
            h_rebinned->SetBinContent(i, weightedSum / weightSum);
            h_rebinned->SetBinError(i, sqrt(errorSum) / weightSum);
        } else {
            h_rebinned->SetBinContent(i, 0);
            h_rebinned->SetBinError(i, 0);
        }
    }

    return h_rebinned;
}
std::vector<double> getAverage(const TH1D *published, const TH1D *my) {
  std::vector<double> v;
  for (int ii=0; ii<hist_rig_highZ->GetNbinsX(); ii++) {
    v.push_back( (published->GetBinContent(ii)+my->GetBinContent(ii))/2 );
  }
  return v;
}
std::vector<double> getAverageErrors(TH1D *published, TH1D *my) {
  std::vector<double> v;
  for (int ii=0; ii<published->GetNbinsX(); ii++) {
    v.push_back( (published->GetBinError(ii)+my->GetBinError(ii) )/2 );
  }
  return v;
}
std::vector<double> getWeightedAverage(const TH1D *published, const TH1D *my) {
    std::vector<double> v;
    for (int ii = 1; ii <= published->GetNbinsX(); ii++) { // Bins start at 1 in ROOT
        double y1 = published->GetBinContent(ii);
        double y2 = my->GetBinContent(ii);
        double e1 = published->GetBinError(ii);
        double e2 = my->GetBinError(ii);

        // Avoid division by zero
        if (e1 == 0 || e2 == 0) {
            v.push_back(0);
            continue;
        }

        double w1 = 1.0 / (e1 * e1);
        double w2 = 1.0 / (e2 * e2);
        double weightedMean = (w1 * y1 + w2 * y2) / (w1 + w2);

        v.push_back(weightedMean);
    }
    return v;
}
std::vector<double> getWeightedAverageErrors(const TH1D *published, const TH1D *my) {
    std::vector<double> v;
    for (int ii = 1; ii <= published->GetNbinsX(); ii++) { // Bins start at 1 in ROOT
        double e1 = published->GetBinError(ii);
        double e2 = my->GetBinError(ii);

        // Avoid division by zero
        if (e1 == 0 || e2 == 0) {
            v.push_back(0);
            continue;
        }

        double w1 = 1.0 / (e1 * e1);
        double w2 = 1.0 / (e2 * e2);
        double weightedError = sqrt(1.0 / (w1 + w2));

        v.push_back(weightedError);
    }
    return v;
}
double get_R_from_flux(unsigned int charge, double low_edge, double high_edge) {
  auto qi = GetQYanFlux(charge);
  //-------Creating model-----------------
  //auto flux_model = autospline(qi,0.8,3000.,3,4);
  auto flux_model = spfit(qi, 4, 2., 3000.);
  flux_model->SetRange(2.,3000);

  auto recoveredFlux = BuildRecoveredFlux(flux_model, 2., 3000.);

  // Define weighted φ(R) * R
  TF1* weightedFlux = new TF1("weightedFlux",
      [recoveredFlux](double *x, double *) {
          return x[0] * recoveredFlux->Eval(x[0]);  // R * φ(R)
      }, 
      low_edge, high_edge, 0);

  //Compute the integral over [low_edge, high_edge]
  double integral = weightedFlux->Integral(low_edge, high_edge);
  double width = recoveredFlux->Integral(low_edge, high_edge);
  return integral/width;
}


#endif