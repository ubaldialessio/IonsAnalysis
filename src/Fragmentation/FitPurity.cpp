#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"
#include "SplineUtility.h"

int Z;
double A,mass;
const double prod_min = 0.9957, prod_max = 2001;
double RigToEkn(double Rig, double Z, double A, double M) {
  return (sqrt((Rig * Z) * (Rig * Z) + M * M) - M) / A;
}
double EknToRig(double Ekn, double Z, double A, double M) {
  // return sqrt((Ekn * A + M) * (Ekn * A + M) - (M * M)) / Z;
  return sqrt(A * Ekn * (A * Ekn + 2. * M)) / Z;
}
//original
double force_field(double* x, double* p) {
  double Rig      = x[0];
  double Norm     = p[0];
  double phi      = p[1];
  double exponent = p[2];

  double Etot = RigToEkn(Rig, Z, A, mass) + mass;

  double SM_factor = (Etot * Etot - mass * mass) / ((Etot + phi) * (Etot + phi) - mass * mass);

  double flux = pow(EknToRig(Etot - mass + phi, Z, A, mass), exponent);

  return Norm * flux * SM_factor;
}
struct PurityGraph {
    TGraphErrors *purity;
    TString label;
};
struct ContributionGraph {
    std::vector<TGraphErrors*> components;  // each bkg fraction
    std::vector<TString> labels;           // names for legend
};
TGraphErrors *getPurity(unsigned int charge, TString type, ContributionGraph* contrib=nullptr);
void fitPurity(std::vector<PurityGraph> pur, unsigned int charge, ContributionGraph* contrib=nullptr);
double get_x_from_flux(unsigned int charge, double low_edge, double high_edge);
double weightedIntegral(TH1D* hist, double xmin, double xmax);
struct PurityValues {
    double value;
    double value_err;
    double bin_low;
    double bin_high;
};
std::pair<double,double> rangeFromCharge(unsigned int charge);

int main(int argc, char *argv[]) {
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    TH1::SetDefaultSumw2(true);
	if (argc < 2) {
		printf("Usage: \n");
	    printf("%s <charge> \n", argv[0]);
	    return 1;
	}
    unsigned int charge=atoi(argv[1]);

    ContributionGraph contrib;
    std::vector<PurityGraph> pur = {
        { getPurity(charge, "bin_center", nullptr), "Bin Center" },
        { getPurity(charge, "int",        &contrib), "Integral"  }
    };
    
    fitPurity(pur, charge, &contrib);

    //fitPurity(purity, charge);
}
TGraphErrors* getPurity(unsigned int charge, TString type, ContributionGraph* contrib) {
    const TString nucleus = getIonName(charge);
    const TString basePath = "../Fragmentation/BelowL1/Fractions/";
    const TString inpTxt = basePath + nucleus + "_purity.txt";
    std::ifstream infile(inpTxt);

    std::vector<double> x, y, e_x, e_y;
    // Etichette attese per i contributi: (Z-1, Z+1, Z+2)
    std::vector<TString> contribLabels = {
        Form("%s (Z-1)", getIonPath(charge-1).Data()),
        Form("%s (Z+1)", getIonPath(charge+1).Data()),
        Form("%s (Z+2)", getIonPath(charge+2).Data())
    };

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        double val, err, bin_low, bin_high;
        iss >> val >> err >> bin_low >> bin_high;

        double x_from_flux = (type == "bin_center")
                                ? 0.5 * (bin_high + bin_low)
                                : get_x_from_flux(charge, bin_low, bin_high);

        // background totale
        x.push_back(x_from_flux);
        y.push_back(1-val);
        e_x.push_back(0.);
        e_y.push_back(err);

        // contributi singoli
        if (contrib) {
            double frac, frac_err;
            int idx = 0;
            while (iss >> frac >> frac_err) {
                if (contrib->components.size() <= (size_t)idx) {
                    contrib->components.push_back(new TGraphErrors());

                    // assegna label specifica in base all'ordine
                    if (idx < (int)contribLabels.size()) {
                        contrib->labels.push_back(contribLabels[idx]);
                    } else {
                        contrib->labels.push_back(Form("Bkg%d", idx+1));
                    }
                }

                int pIndex = contrib->components[idx]->GetN();
                contrib->components[idx]->SetPoint(pIndex, x_from_flux, 1-frac);
                contrib->components[idx]->SetPointError(pIndex, 0.0, frac_err); // <-- now plot error
                idx++;
            }
        }
    }
    infile.close();

    TGraphErrors* graph = new TGraphErrors(x.size(),
                                           x.data(), y.data(),
                                           e_x.data(), e_y.data());
    return graph;
}
void fitPurity(std::vector<PurityGraph> pur, unsigned int charge, ContributionGraph* contrib) {
    //Check the integral one
    int nPoints = pur[1].purity->GetN();
    for (int i = 0; i < nPoints; ++i) {
        double x, y;
        pur[1].purity->GetPoint(i, x, y);
        double ex = pur[1].purity->GetErrorX(i);
        double ey = pur[1].purity->GetErrorY(i);
        
        std::cout << "Point " << i 
                << ": x = " << x << " ± " << ex
                << ", y = " << y << " ± " << ey << std::endl;
    }

    //double rMax = 1000.;
    //double rMin = 0.872587; P
    //double rMax = 744.694; P
    double rMin = 1.00135;
    double rMax = 486.939;


    TF1 *F_to_Si_function_delta2 = new TF1("F_to_Si_function_delta2", "1 - [0]*pow(x/185., +0.19)",10.,rMax);
    F_to_Si_function_delta2->SetLineColor(kGreen);
    F_to_Si_function_delta2->SetParameters(20.);

    //Output ROOT
    const TString nucleus = getIonName(charge);
    TString fileOut = "../Fragmentation/BelowL1/Fractions/contamination_"+nucleus+".root";
    TFile *output = new TFile(fileOut.Data(), "RECREATE");

    //Output PDF
    TCanvas *c1 = new TCanvas("c1", "", 1024, 640);
    c1->SetLogx(); 
    c1->SetGrid();
    c1->SaveAs(Form("../output/fitPurity_%d.pdf[",charge));
    //fraction->Fit(purityFit,"LR");
    //fraction->Draw("AP");
    auto leg = new TLegend(0.1,0.1,0.4,0.3);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);   
    leg->SetTextSize(0.045);
    leg->SetTextFont(62);
    bool first = true;
    for (int i=0; i<pur.size(); i++) {
        pur[i].purity->SetMarkerSize(1.5);
        pur[i].purity->SetMarkerColor(kRed+3*i);
        pur[i].purity->SetLineColor(kRed+3*i);
        pur[i].purity->SetTitle(Form("Z = %u;R (GV);Purity (%%)", charge));
        if (first) {
            pur[i].purity->SetMarkerStyle(24);
            pur[i].purity->Draw("AP");
            pur[i].purity->GetXaxis()->SetRangeUser(0.001,300000);
            first = false;
        } else {
            pur[i].purity->SetMarkerStyle(33);
            pur[i].purity->Draw("P SAME");
            pur[i].purity->GetXaxis()->SetRangeUser(0.001,300000);
        }
        leg->AddEntry(pur[i].purity,pur[i].label.Data());
    }
    leg->Draw("SAME");
    c1->SaveAs(Form("../output/fitPurity_%d.pdf",charge));


    //---Drawing only the integral one and the various fit
    TH1D* base = (TH1D*)hist_rig_highZ->Clone("base");
    formatAxis(base,1);
    base->GetXaxis()->SetRangeUser(0.001,1000.);
    base->SetTitle(Form("Z = %u;R (GV);Purity (%%)", charge));
    base->SetLabelFont(62,"XYZ");
    base->SetTitleFont(62,"XYZ");
    //pur[1].purity->Fit(F_to_Si_function_delta2,"LR");

    auto leg1 = new TLegend(0.1,0.15,0.4,0.40);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);   
    leg1->SetTextSize(0.045);
    leg1->SetTextFont(62);
    leg1->SetLineColor(0);      // Removes the frame around the legend box (around the whole legend, not the entry box)
    leg1->SetFillStyle(0); 

    std::vector<double> xvals, yerr68;
    int Npoints = pur[1].purity->GetN();
    for (int i = 0; i < nRbins_HighZ - 1; ++i) {
        double x_center = 0.5 * (Rbins_HighZ[i] + Rbins_HighZ[i+1]);
        xvals.push_back(x_center);
    }
    //custom knots
    //std::vector<double> knots = {rMin, 15, 1000}; P
    std::vector<double> knots = {rMin,17};
    //TF1* purityFit = spfit(pur[1].purity, 2, rMin, rMax, knots, "", &xvals, &yerr68);
    //TF1* purityFit = spfit_autoBIC(pur[1].purity, rMin, rMax, 2,5,"lin", -1, &xvals, &yerr68);
    TF1* purityFit = new TF1("constant", "pol0", rMin, 1000);
    purityFit->FixParameter(0, 0.99);
    pur[1].purity->Fit(purityFit);
    purityFit->SetNpx(10000);
    double xsplit = 12; // metti il valore desiderato
    double xmin = rMin, xmax = 1000; // range

    // definisco la funzione: parabola centrata in xsplit e poi costante C
    /*TF1* purityFit = new TF1("purityFit",
    Form("((x<%f) ? ([0]*(x-%f)*(x-%f) + [1]) : ([1]))", xsplit, xsplit, xsplit),
    xmin, xmax);

    // inizializza a e C (valori di partenza ragionevoli)
    purityFit->SetParameter(0, 1e-2); // a (se la curva scende, a<0)
    purityFit->SetParameter(1, 0.95);  // C plateau*/

    // limiti se vuoi (opzionale)
    /*purityFit->SetParLimits(0, -1e-2, 1e-2);
    purityFit->SetParLimits(1, 0.5, 1.0);*/
    //pur[1].purity->Fit(purityFit);
    //auto purityFit = spfit_autoBIC(pur[1].purity, rMin, 1000, 1, 3, "lin", 20, &xvals, &yerr68);
    purityFit->SetMarkerColor(kRed+3);

    // Optionally plot confidence band using TraphAsymmErrors
    /*TGraphAsymmErrors* band = new TGraphAsymmErrors(xvals.size());
    for (size_t i = 0; i < xvals.size(); ++i) {
        double x = xvals[i];
        //int bin = purityFit->FindBin(x);
        double y = purityFit->Eval(x);
        band->SetPoint(i, x, y);
        band->SetPointError(i, 0, 0, yerr68[i], yerr68[i]); // symmetric error
    }
    band->SetFillColorAlpha(kRed, 0.3);

    for (int i = 0; i < Npoints; ++i) {
        double x,y;
        pur[1].purity->GetPoint(i, x, y);
        std::cout << "x = " << x << ", y = " << y << std::endl;
    }*/

    purityFit->SetRange(rMin,3000);
    // Compute scaling factor to ensure continuity at R = rMax
    /*double R_match = 10.5;
    double scale = purityFit->Eval(R_match) / F_to_Si_function_delta2->Eval(R_match);

    TF1 *combinedPurity = new TF1("combinedPurity",
        [=](double *x, double *) {
            double R = x[0];
            if (R <= R_match)
                return purityFit->Eval(R);
            else
                return scale * F_to_Si_function_delta2->Eval(R);
        },
        1., 3000., 0);
        
    combinedPurity->SetLineColor(kBlue);
    F_to_Si_function_delta2->SetRange(1.,3000);*/

    leg1->AddEntry(pur[1].purity,"Purity","LP");
    leg1->AddEntry(purityFit,Form("%s",purityFit->GetTitle()),"L");
    
    base->SetStats(0);
    auto range = rangeFromCharge(charge);
    base->GetYaxis()->SetRangeUser(range.first,range.second);
    base->Draw("HIST");
    leg1->Draw();
    pur[1].purity->Draw("SAME P");
    pur[1].purity->GetXaxis()->SetRangeUser(0.001,1000.);
    purityFit->Draw("SAME L");
    //band->Draw("3 same");
    //F_to_Si_function_delta2->Draw("SAME");
    //combinedPurity->Draw("SAME");
    purityFit->Draw("SAME");
    c1->SaveAs(Form("../output/fitPurity_%d.pdf",charge));

    // dopo aver disegnato i grafici principali
    if (contrib) {
        int colorIdx = 2;
        for (size_t j=0; j<contrib->components.size(); ++j) {
            auto g = contrib->components[j];
            g->SetMarkerStyle(20);
            setFriendColor(g,j);
            g->Draw("P SAME");
            leg1->AddEntry(g, contrib->labels[j], "P");
            colorIdx++;
        }
    }
    
    c1->SaveAs(Form("../output/fitPurity_%d.pdf",charge));
    c1->SaveAs(Form("../output/fitPurity_%d.pdf]",charge));

    auto Purity = (TH1D*)hist_rig_highZ->Clone("Purity");
    for (int i=1; i<=Purity->GetNbinsX(); i++) {
        double xlow  = Purity->GetXaxis()->GetBinLowEdge(i);
        double xhigh = Purity->GetXaxis()->GetBinUpEdge(i);
        double val   = purityFit->Integral(xlow, xhigh) / (xhigh - xlow);
        Purity->SetBinContent(i, val);
    }

    output->WriteTObject(purityFit,"purityFit");
    output->WriteTObject(Purity,"Purity");


    output->Close();
}
double get_x_from_flux(unsigned int charge, double low_edge, double high_edge) {
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

    MultiplyByXPower(qi,-2.7);

    //Output PDF
    TString out = Form("../output/ff_fit_%d.pdf",charge);
    TCanvas *c1 = new TCanvas("c1", "", 1024, 640);
    c1->SetLogx(); 
    c1->SetLogy();
    c1->SetGrid();
    qi->GetXaxis()->SetTitle("R (GV)");
    qi->GetYaxis()->SetTitle("#phi R^{2.7} (m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    flux_model->SetMarkerStyle(20);
    flux_model->SetTitle(Form("Z = %u;R (GV);Purity (%%)", charge));
    //qi->GetYaxis()->SetRangeUser(-0.2,0.8);
    qi->SetMarkerColor(kRed);
    qi->SetLineColor(kRed);
    qi->Draw("P");
    //flux_model->Draw("SAME");
    recoveredFlux->Draw("SAME");
    c1->SaveAs(out.Data());

    // Compute the integral over [low_edge, high_edge]
    double integral = weightedFlux->Integral(low_edge, high_edge);
    double width = recoveredFlux->Integral(low_edge, high_edge);
    std::cout << " ----- Integral = " << integral << std::endl;
    return integral/width;

    //double integral = weightedIntegral(flux_model, low_edge, high_edge);
}
double weightedIntegral(TH1D* hist, double xmin, double xmax) {
    double result = 0.0;
    double content= 0.0;
    double bin_center = (xmax+xmin)*0.5;
    for (int i = hist->FindBin(xmin); i<= hist->FindBin(xmin); i++ ) {
        content+=hist->GetBinContent(i);
    }

    /*for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double binCenter = hist->GetBinCenter(i);
        double binContent = hist->GetBinContent(i);

        // Check if bin center is inside the desired x-range
        if (binCenter >= xmin && binCenter <= xmax) {
            result += binCenter * binContent;
        }
    }*/

    return bin_center*content;
}
std::pair<double,double> rangeFromCharge(unsigned int charge) {
    double ymin,ymax;
    switch(charge) {
        case 14:
            ymin = 0.99;
            ymax = 1.01;
            break;
        case 15:
            ymin = 0.4;
            ymax = 1.1;
            break;
        case 16:
            ymin = 0.85;
            ymax = 1;
            break;
        case 18:
            ymin = 0.88;
            ymax = 1.0;
            break;
        case 20:
            ymin = 0.65;
            ymax = 1.05;
            break;
    }
    return std::make_pair(ymin,ymax);
}