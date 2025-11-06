#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "SplineUtility.h"
//Roofit headers
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooHistConstraint.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooCmdArg.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooCrystalBall.h"
#include "RooHist.h"
#include "RooProduct.h"
//additional header
#include <string_view>
#include <TLatex.h>
#include "TRandom3.h"
#include <TFitResult.h>
static constexpr auto EXENAME = R"(fitTemplates)";
using namespace RooFit;

// data structure: [charge][layer/distribution][histogram]
std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> templatesByCharge;
void fitHistograms(std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>>& templatesByCharge,
                   unsigned int charge, unsigned int n, TString option);
std::map<std::string, std::vector<TH1D*>> getTemplates(unsigned int charge, TString binType); 
std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> loadTemplates(unsigned int charge, unsigned int n, TString binType);
double computeFraction(const RooRealVar &coeffSig, const RooArgList &bgCoeffList,
                       std::vector<double> &den, std::vector<double> &bkg,
                       std::vector<double> &h, std::vector<double> &h_err,
                       double sigIntegral, std::vector<double> bkgIntegral,
                       double err_sig, std::vector<double> err_bkg,
                       RooFitResult* fitRes,
                       std::vector<std::vector<double>> &singleContribution,
                       std::vector<std::vector<double>> &singleContributionErr,
                       int uncertaintyOption=1, // 1 = fit error only, 2 = template stats, 3 = quadrature
                       int nToys=0);
double get_x_from_flux(unsigned int charge, double low_edge, double high_edge);
std::pair<double,double> getBkgRange(unsigned int charge);
std::pair<double,double> getSigRange(unsigned int charge);
TH1D* shiftHistogram(const TH1D* inputHist, double shiftAmount);
void fitPurity(std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>>& templatesByCharge,
               unsigned int charge, unsigned int n, TString binType,
               int firstNCombine, int groupSize, double lowerThresholdGV, TString rebinTemplates);
double initialFractionValue(unsigned int bkgCharge, unsigned int charge, double rig,std::vector<double> Rbins_purity);
TH1D* removeZeroBins(const TH1D* originalHist);
int findRigidityBin(double rigidity, const std::vector<double>& bins);
TH1D* shiftHistogramResample(const TH1D* inputHist, double shiftAmount);
TH1D* shiftHistogramInterpolation(const TH1D* inputHist, double shiftAmount);
TH1D* smoothHistogram(const TH1D* inputHist, int nSamples = 1e5);
std::tuple<double, int, double> computeChi2FromFrame(RooPlot* frame, const char* modelName, const char* dataHistName, int nFitParams);
std::tuple<double, int, double> computeChi2WithSigmaTemplate(
    RooRealVar x,
    RooAddPdf model,
    TH1D *data,
    TH1D *sigTemplate,
    std::vector<TH1D*> bkgTemplate,
    int nFitParams,
    RooPlot* frame, const char* modelName, const char* dataHistName);
std::tuple<double, double, double, double, int> computeChi2WithBootstrap(
    RooRealVar x,
    TH1D* data,
    TH1D* histSignal,
    std::vector<TH1D*> histBkg,
    unsigned int charge,
    int binIndex,
    const std::vector<double>& Rbins_purity,
    int nFitParams,
    int nBootstrap = 2);
std::tuple<RooAddPdf*,RooPlot*> buildAndFitModel(
    RooRealVar x,
    TH1D* data,
    TH1D* histSignal,
    std::vector<TH1D*> histBkg,
    unsigned int charge,
    int i,
    const std::vector<double>& Rbins_purity,
    RooFitResult*& fitResult,
    int bin,
    int nBoot);
double stepByCharge(unsigned int charge);
std::pair<double,double> rangeShiftByCharge(unsigned int charge);   
TH1D* RooHistToTH1D(const RooHist* rh, const TString& name = "pullHist");
struct frameValues {
    double chi2;
    int ndf;
    double chi2ndf;
    size_t bin;
    double shift;
    double shift_err;
};
void formatFrameSinglePage(std::unique_ptr<RooPlot> &frame, frameValues &,std::vector<double> Rbins_purity, TCanvas &c,
                           unsigned int charge);
std::pair<double,double> cutLimits(unsigned int charge);
static std::pair<double,double> computeL1ChargeEfficiencyFromPdf(
        RooHistPdf* signalPdf,
        RooRealVar& x,
        RooDataHist dataRoo,   // histData convertito in RooDataHist
        double cutLow,
        double cutHigh);
std::vector<TH1D*> rebinAndCombine(
    const std::vector<TH1D*>& histVec,
    const std::vector<double>& Rbins,
    int firstNCombine,
    int groupSize,
    double upperThresholdGV,
    int charge,
    const std::string& type,
    std::vector<double>& newRbins);
std::vector<TH1D*> autoRebinAndCombine(
    const std::vector<TH1D*>& histVec,
    const std::vector<double>& Rbins,
    int charge,
    const std::string& type,
    std::vector<double>& newRbins,
    int groupSize = 3,
    double lowerThresholdGV = 85.,
    double upperThresholdGV = 3000.);
int smoothAmount(unsigned int charge);
std::map<unsigned int, std::vector<double>> charge14_initialFractions = {
    {13, {0.001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001}},
    {14, {3, 3, 3, 3, 3, 3, 2.8, 2.5, 2.2, 2.0, 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.}},
    {15, {0.001, 0.001, 0.001, 0.005, 0.005, 0.005, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.005, 0.005, 0.005, 0.01, 0.01, 0.01, 0.01, 0.01}},
    {16, {0.003, 0.003, 0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001}}
};
std::map<unsigned int, std::vector<double>> charge15_initialFractions = {
    {14, {50, 50, 50, 50, 50, 0.01, 0.005, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,40,40}},
    {15, {20000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ,20 ,20 ,10,10,10 }},
    {16, {700, 1, 1, 1, 1, 1, 1, 1, 1, 1, 10, 1, 10 ,0.01 ,0.001 ,0.001,0.001,0.001 }},
    {17, {70, 10, 10, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.055, 0.055, 0.055, 0.055, 0., 0.,0.,.0}}
};
/*std::map<unsigned int, std::vector<double>> charge15_initialFractions = {
    {14, std::vector<double>(73, 40)},
    {15, std::vector<double>(73, 500)},
    {16, std::vector<double>(73, 300)},
    {17, std::vector<double>(73, 20)}
};*/
/*std::map<unsigned int, std::vector<double>> charge15_initialFractions = {
    {14, {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}},
    {15, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
    {16, {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8}},
    {17, {10e-4, 1, 1, 1, 0.5, 0,5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}}
};*/
std::map<unsigned int, std::vector<double>> charge16_initialFractions = {
    {15, {0.002, 0.002, 0.002, 0.05, 0.01, 0.01, 0.01, 0.005, 0.005, 0.002, 0.001, 0.001, 0.001, 0.0001, 0.0001, 0.0001, 0.0001}},
    {16, {30, 30, 3, 3, 3, 3, 3, 3, 3, 3, 3., 3., 3., 3., 3., 3., 3., 3.,3., 3., 2.}},
    {17, {0.001, 0.001, 0.001, 0.01, 0.02, 0.01, 0.01, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.001, 0.002, 0.001, 0.01, 0.0, 0.01, 0.01}},
    {18, {0.01, 0.001, 0.01, 0.004, 0.004, 0.002, 0.002, 0.001, 0.0005, 0.0005, 0.0005, 0.001, 0.001, 0.0, 0.0, 0., 0.00, 0.001, 0.001, 0.001, 0.001}}
};
std::map<unsigned int, std::vector<double>> charge18_initialFractions = {
    {17, {0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.002, 0.00, 0.00, 0.00, 0.000, 0.000, 0.000, 0.0001}},
    {18, {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10., 10., 10., 10., 10., 10., 10., 10.,10., 10., 2.}},
    {19, {0.02, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.001, 0.002, 0.001, 0.01, 0.01, 0.01, 0.01}},
    {20, {0.01, 0.005, 0.01, 0.01, 0.01, 0.002, 0.002, 0.001, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0, 0.0, 0., 0.001, 0.001, 0.001, 0.001, 0.001}}
};
std::map<unsigned int, std::vector<double>> charge20_initialFractions = {
    {19, {0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 2, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,0.6}},
    {20, {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8}},
    {21, {2.5, 2.5, 2.5, 2.5, 2.5, 1, 0.5, 0.5, 1.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 1}},
    {22, {0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.3, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.5}}
};

struct FitStatus {
    int status;
    int covQual;
    double chiSquare;
    FitStatus(int s, int cq, double chiS) : status(s), covQual(cq), chiSquare(chiS) {}
};

int main(int argc, char *argv[]) {
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    TH1::SetDefaultSumw2(true);

    if (argc < 7) {
        printf("Usage:\n");
        printf("%s <charge> <binType> <rebinTemplates> <firstNCombine> <groupSize> <upperThresholdGV>\n", argv[0]);
        printf(" - binType: pub / mine\n");
        printf(" - rebinTemplates: y / n\n");
        printf(" - firstNCombine: number of lowest bins to merge (e.g. 3)\n");
        printf(" - groupSize: group size for middle bins (e.g. 2)\n");
        printf(" - lowerThresholdGV: merge all bins above this rigidity (e.g. 3000)\n");
        return 1;
    }
    unsigned int charge     = atoi(argv[1]);
    TString binType         = argv[2];
    TString rebinTemplates  = argv[3];
    int firstNCombine       = atoi(argv[4]);
    int groupSize           = atoi(argv[5]);
    double lowerThresholdGV = atof(argv[6]);
    if (rebinTemplates == "y") {
        if (groupSize < 1) {
            std::cerr << "Error: groupSize must be >= 1\n";
            return 1;
        }
    }
    auto templates = loadTemplates(charge, 2, binType);
     fitPurity(templates, charge, 2, binType, firstNCombine, groupSize, lowerThresholdGV,rebinTemplates);

}
std::map<std::string, std::vector<TH1D*>> getTemplates(unsigned int charge, TString binType) {
    std::map<std::string, std::vector<TH1D*>> templates;

    TString path1 = "../Fragmentation/BelowL1/";
    TString nucleus = getIonName(charge);
    TString fileName;
    if (binType == "pub")
        fileName = path1 + nucleus + "/" + nucleus + "_publishedBinning.root";
    else if (binType == "mine")
        fileName = path1 + nucleus + "/" + nucleus + ".root";

    std::cout << "Getting templates: " << fileName << std::endl;

    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        printf("❌ Errore nell'aprire il file\n");
        return templates;
    }

    // Directories inside IL1
    TString baseDir = "IL1/";
    TString dirL1 = baseDir + "L1TemplateList_" + TString::Format("%d", charge);
    TString dirL2 = baseDir + "L2TemplateList_" + TString::Format("%d", charge);
    TString dirD  = baseDir + "L1Distribution_" + TString::Format("%d", charge);

    std::vector<std::pair<TString, std::string>> dirs = {
        {dirL1, "L1"},
        {dirL2, "L2"},
        {dirD,  "D"}
    };

    for (auto &entry : dirs) {
        TDirectory *dir = dynamic_cast<TDirectory*>(file->Get(entry.first));
        if (!dir) {
            printf("⚠️  Errore nel recuperare la directory %s\n", entry.first.Data());
            continue;
        }

        TIter nextkey(dir->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)nextkey())) {
            TObject *obj = key->ReadObj();
            TH1D *hist = dynamic_cast<TH1D*>(obj);
            if (hist) {
                templates[entry.second].push_back((TH1D*)hist->Clone());
            }
        }
    }

    return templates;
}

std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> loadTemplates(unsigned int charge, unsigned int n, TString binType) {
    std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> templatesByCharge;
    //Original with charge, charge+1, charge+2
    /*for (unsigned int i = 0; i <= n; i++) {
        unsigned int currentCharge = charge + i;
        templatesByCharge[currentCharge] = getTemplates(currentCharge); 
    }*/
    //Now including also charge -1 
    std::vector<unsigned int> allowedZ = {charge-1,charge,charge+1,charge+2};
    for (auto k : allowedZ)
        templatesByCharge[k] = getTemplates(k,binType); 
    return templatesByCharge;
}
void fitPurity(std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>>& templatesByCharge,
               unsigned int charge, unsigned int n, TString binType,
               int firstNCombine, int groupSize, double lowerThresholdGV, TString rebinTemplates) {
    auto ion = getIonName(charge).Data();
    TString outPath = Form("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/%s/%s_chi.root",
                          ion,ion);
    auto outfile = new TFile(outPath.Data(), "recreate");
    
    //binning based on charge
    std::vector<double> Rbins_purity;
    if (binType=="mine") {
        std::cout << "Using bin for P" << std::endl;
        Rbins_purity = {0.8, 3.64, 5.37, 7.76, 11.00, 15.3, 21.1, 38.9, 70.00, 150.00, 300.00, 500.00, 1000.00};
    } else {
        std::cout << "Using bin for Si, S" << std::endl;
        Rbins_purity = {
            0.8, 1.00, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.40, 2.67,
            2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.90, 6.47, 7.09,
            7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0, 14.1, 15.3, 16.6,
            18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1,
            38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9,
            80.5, 86.5, 93.0, 100., 108., 116., 125., 135., 147., 160.,
            175., 192., 211., 233., 259., 291., 330., 379., 441., 525.,
            660., 860., 1200., 3000.
        };
    }
    int nRbins_purity = 74;
    std::vector<double> Rbins_purity_rebinned;
    if (rebinTemplates == "y") {
        for (const auto& [charge, templateMap] : templatesByCharge) {
            for (const auto& [type, histVec] : templateMap) {
                std::vector<TH1D*> rebinnedVec;

                for (size_t i = 0; i < histVec.size(); i += groupSize) {
                    TH1D* hsum = (TH1D*) histVec[i]->Clone(
                        Form("%s_rebinned_charge%d_bin%zu", type.c_str(), charge, i / groupSize));
                    hsum->SetDirectory(0);

                    // Add the remaining histograms in the group
                    for (size_t j = 1; j < groupSize; ++j) {
                        if (i + j < histVec.size()) {
                            hsum->Add(histVec[i + j]);
                        }
                    }

                    rebinnedVec.push_back(hsum);
                }

                templatesByCharge[charge][type] = rebinnedVec;
            }
        }

        // Build new binning
        for (size_t i = 0; i < Rbins_purity.size(); i += groupSize) {
            Rbins_purity_rebinned.push_back(Rbins_purity[i]);
        }

        // Ensure final upper edge is pushed only if not already included
        if (Rbins_purity_rebinned.back() != Rbins_purity.back()) {
            Rbins_purity_rebinned.push_back(Rbins_purity.back());
        }
        Rbins_purity = Rbins_purity_rebinned;
    }

    if (rebinTemplates == "y") {
        double upperThresholdGV = 3000.0;

        // Find index where R > lowerThresholdGV
        size_t startIdx = 0;
        while (startIdx < Rbins_purity.size() && Rbins_purity[startIdx] <= lowerThresholdGV) {
            ++startIdx;
        }

        // Find index where R > upperThresholdGV
        size_t endIdx = startIdx;
        while (endIdx < Rbins_purity.size() && Rbins_purity[endIdx] <= upperThresholdGV) {
            ++endIdx;
        }

        for (auto& [charge, templateMap] : templatesByCharge) {
            for (auto& [type, histVec] : templateMap) {
                std::vector<TH1D*> combinedVec;

                // --- 1) Combine the first N bins
                if (histVec.size() >= static_cast<size_t>(firstNCombine)) {
                    TH1D* hfirst = (TH1D*) histVec[0]->Clone(
                        Form("%s_combined_charge%d_bin0to%d", type.c_str(), charge, firstNCombine - 1));
                    hfirst->SetDirectory(0);
                    for (int i = 1; i < firstNCombine; ++i) {
                        hfirst->Add(histVec[i]);
                    }
                    combinedVec.push_back(hfirst);
                } else {
                    // Fallback: not enough bins to combine
                    for (auto* h : histVec)
                        combinedVec.push_back(h);
                    templateMap[type] = combinedVec;
                    continue;
                }

                // --- 2) Keep bins from index firstNCombine up to startIdx - 1
                for (size_t i = firstNCombine; i < startIdx - 1 && i < histVec.size(); ++i) {
                    combinedVec.push_back(histVec[i]);
                }

                // --- 3) Combine bins between startIdx and endIdx - 1
                if (startIdx - 1 < histVec.size()) {
                    TH1D* hsum = (TH1D*) histVec[startIdx - 1]->Clone(
                        Form("%s_combined_charge%d_%.1fto%.0fGV", 
                            type.c_str(), charge, lowerThresholdGV, upperThresholdGV));
                    hsum->SetDirectory(0);
                    for (size_t i = startIdx; i < endIdx - 1 && i < histVec.size(); ++i) {
                        hsum->Add(histVec[i]);
                    }
                    combinedVec.push_back(hsum);
                }

                // --- 4) Add remaining histograms above upperThresholdGV
                for (size_t i = endIdx - 1; i < histVec.size(); ++i) {
                    combinedVec.push_back(histVec[i]);
                }

                // --- 5) Replace old vector with combined version
                templateMap[type] = combinedVec;
            }
        }


        // --- Build new Rbins_purity
        std::vector<double> newRbins;

        // --- 1) Combine first N bins
        if (Rbins_purity.size() >= static_cast<size_t>(firstNCombine + 1)) {
            // Add lower edge of first bin
            newRbins.push_back(Rbins_purity[0]);
            // Add upper edge of the last combined bin
            newRbins.push_back(Rbins_purity[firstNCombine]);
        } else {
            // Fallback: not enough bins to combine
            newRbins.push_back(Rbins_purity.front());
            newRbins.push_back(Rbins_purity.back());
            Rbins_purity = newRbins;
            return;
        }

        // --- 2) Add edges from firstNCombine to startIdx - 1
        for (size_t i = firstNCombine + 1; i < startIdx && i < Rbins_purity.size(); ++i) {
            newRbins.push_back(Rbins_purity[i]);
        }

        // --- 3) Add upper edge of combined bin from startIdx to endIdx
        if (endIdx - 1 < Rbins_purity.size()) {
            newRbins.push_back(Rbins_purity[endIdx - 1]);
        } else {
            newRbins.push_back(Rbins_purity.back());
        }

        // --- 4) Add remaining edges after upperThresholdGV
        for (size_t i = endIdx; i < Rbins_purity.size(); ++i) {
            newRbins.push_back(Rbins_purity[i]);
        }

        // --- 5) Replace
        Rbins_purity = newRbins;

    }

    /*const int nP_rigBins = 13; 
	const double P_rigBins[nP_rigBins] = {0.8,3.64,5.37,7.76,11.00,15.3,21.1,38.9,70.00,150.00,300.00,500.00,1000.00};*/
    // --- 1. Check that the main charge exists and has the necessary templates.
    if (templatesByCharge.find(charge) == templatesByCharge.end()) {
        std::cerr << "Missing templates for charge " << charge << "\n";
        return;
    }
    const auto &mainTemplates = templatesByCharge.at(charge);
    if (mainTemplates.find("L1") == mainTemplates.end() || mainTemplates.at("L1").empty()) {
        std::cerr << "No L1 templates for charge " << charge << "\n";
        return;
    }
    if (mainTemplates.find("L2") == mainTemplates.end() || mainTemplates.at("L2").empty()) {
        std::cerr << "No L2 templates for charge " << charge << "\n";
        return;
    }
    
    //double chMin = 2, chMax = 2;   
    auto cutLim = cutLimits(charge);
    double chMin = cutLim.first; double chMax = cutLim.second;


    // --- 2. Set up the templates for data and signal.
    auto dataTemplates   = templatesByCharge[charge]["D"];
    auto signalTemplates = templatesByCharge[charge]["L2"];
    const int nBins = dataTemplates.size();

    // --- 3. Define the shift scan
    double step = stepByCharge(charge);
    double shift_min = rangeShiftByCharge(charge).first, shift_max = rangeShiftByCharge(charge).second;
    int nShiftSteps = static_cast<int>((shift_max - shift_min) / step) + 1;

    // --- Set up the output
    const TString nucleus = getIonName(charge);
    const TString basePath = "../Fragmentation/BelowL1/Fractions/";
    const TString outPDF  = basePath + nucleus + "_purityFit.pdf";
    const TString outTXT = basePath + nucleus + "_purity.txt";
    std::ofstream outTxt(outTXT);
    TCanvas canvas("canvasBin", "Bin", 1024,640);
    canvas.SetLogy();
    canvas.Print((outPDF+"[").Data());
    auto shift_vs_R  = new TGraphErrors();
    shift_vs_R->SetTitle(";R (GV); Best shift (c.u.)");

    TCanvas* c1 = new TCanvas("c1", "Best shift vs R", 1024,640);
    const TString shiftPDF = basePath+nucleus+ "_chi.pdf";
    c1->SaveAs((shiftPDF+"[").Data());

    //debug
    std::vector<TH1D*> bkgVec;
    std::vector<std::unique_ptr<RooFitResult>> bestFits;
    std::unique_ptr<RooFitResult> bestSingleFit;
    std::vector<double> backgroundIntegral;
    std::vector<double> bestBackgroundIntegrals;
    std::vector<double> BkgErr;
    std::vector<double> bestBkgErr;
    double sigIntegral = 0;
    double bestSigIntegral = 0;
    double bestSigErr = 0;
    std::vector<double> sig;
    std::vector<double> h,h_err, den, bkg; 
    std::vector<std::vector<double>> singleContribution,singleContribution_err;
    std::unique_ptr<RooAddPdf> bestModel;

    //smooth amount
    int smt = smoothAmount(charge);


    // --- Choose rigidity threshold (GV) above which to reuse the same signal template
    double Rthr = 120.0; // <--- set your threshold here (GV)

    // Find the first bin index whose bin center is > Rthr.
    // Assumes Rbins_purity has at least nBins+1 entries (bin edges).
    size_t thresholdBinIndex = nBins; // default: no threshold
    for (size_t ib = 0; ib < static_cast<size_t>(nBins); ++ib) {
        if (ib + 1 >= Rbins_purity.size()) break; // safety
        double Rlow = Rbins_purity[ib];
        double Rhigh = Rbins_purity[ib + 1];
        double Rcenter = 0.5 * (Rlow + Rhigh);
        if (Rcenter > Rthr) { thresholdBinIndex = ib; break; }
    }
    int useTemplateIndex =thresholdBinIndex;


    // store L1-efficiency per rigidity bin (filled inside per-bin loop)
    std::vector<double> L1eff;
    std::vector<double> L1effErr;
    L1eff.reserve(nBins);
    L1effErr.reserve(nBins);

    // Loop over bins (each bin corresponds to one pair of data/signal templates)
    for (size_t i = 0; i < static_cast<size_t>(nBins); ++i) {
        bkgVec.clear();

        RooRealVar x("x", "Q_{L1}", charge-2.3, charge+1.7);
        
        // --- Data and Signal RooFit objects.
        TH1D* histData   = dataTemplates.at(i);
        //TH1D* histSignal = (i < signalTemplates.size()) ? signalTemplates.at(i) : signalTemplates.back();

        TH1D* histSignal = nullptr;
        TH1D* histSignalOriginal = nullptr;

        bool aboveThreshold = (i >= thresholdBinIndex);

        if (aboveThreshold) {
            // Choose a fixed template above threshold
            if (useTemplateIndex >= 0 && static_cast<size_t>(useTemplateIndex) < signalTemplates.size()) {
                histSignalOriginal = signalTemplates.at(useTemplateIndex);
            } else if (!signalTemplates.empty()) {
                histSignalOriginal = signalTemplates.back(); // default: reuse last template
            } else {
                histSignalOriginal = nullptr;
            }
        } else {
            // bin-dependent template (or fallback to last if not available)
            if (i < signalTemplates.size()) histSignalOriginal = signalTemplates.at(i);
            else if (!signalTemplates.empty()) histSignalOriginal = signalTemplates.back();
            else histSignalOriginal = nullptr;
        }

        // --- Printing
        TCanvas* cAll = new TCanvas("cAll", "All fits", 1024,640);
        if (i==0) cAll->Print("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"+getIonName(charge)+
            Form("/%s_shift.pdf[",ion));
        cAll->Divide(5, 4, 0.,0.);
        //aspect ratio
		int nx = 5;
		int ny = 4;
		float padW = 1.0 / nx;
		float padH = 1.0 / ny;
		// Loop over pads and set their size and margins
		for (int ix = 0; ix < nx; ++ix) {
			for (int iy = 0; iy < ny; ++iy) {
				int padNumber = ix + 1 + iy * nx;
				cAll->cd(padNumber);
				TPad *pad = (TPad*)gPad;
				
				float xlow = ix * padW;
				float ylow = 1.0 - (iy + 1) * padH;
				float xup;
				xup  = xlow + padW;
				if (ix==nx-1) xup  = xlow + 0.989*padW;
				float yup  = ylow + 0.989*padH;

				pad->SetPad(xlow, ylow, xup, yup);
				//pad->SetMargin(0.1, 0., 0., 0.); // custom margins
			}
		}
        for (int p = 1; p <= nShiftSteps; ++p) {
            cAll->cd(p);
            gPad->SetLogy(1);  // <- Set logY for each sub-pad
        }

        //TH1D* histSignalOriginal = (i < signalTemplates.size()) ? signalTemplates.at(i) : signalTemplates.back();
        // Now make histSignal the working copy (you use histSignal later after shifting)
        if (histSignalOriginal) {
            histSignal = static_cast<TH1D*>(histSignalOriginal->Clone(Form("histSignal_bin%zu", i)));
        } else {
            std::cerr << "No signal template available for bin " << i << " (Rcenter ~ "
                    << 0.5*(Rbins_purity[i] + Rbins_purity[i+1]) << " GV)\n";
            continue;
        }

        double minChi2 = 999999;
        int bestShiftIndex = -1;
        double bestShift = 0;

        // Oggetti da salvare per il best shift
        std::unique_ptr<RooRealVar> bestCoeffSig;
        std::unique_ptr<RooArgList> bestBgCoeffList;
        std::unique_ptr<RooPlot> bestFrame;
        std::unique_ptr<TLegend> bestLegend;
        std::unique_ptr<frameValues> bestFrameValues; 

        auto chi_vs_shift = new TGraph();

        //Loop over shift values
        for (int s = 0; s < nShiftSteps; ++s) {
            double currentShift = shift_min + s * step;

            // Remove zero bins
            histSignal->Smooth(smt);
            histSignal = removeZeroBins(histSignal);

            TH1D* shiftedTemplate = shiftHistogramResample(histSignalOriginal, currentShift);
            TH1D* originalSignal = (TH1D*)histSignal->Clone( Form("OriginalSignal") );
            if (shiftedTemplate) histSignal = shiftedTemplate;

            if (histSignal->GetEntries() == 0) {
                std::cerr << "Signal histogram is empty!" << std::endl;
                continue;
            }
            TString title = histData->GetTitle();
            RooDataHist dataRoo("dataRoo", "Data", RooArgList(x), Import(*histData));
            RooDataHist signalRoo("signalRoo", "Signal", RooArgList(x), Import(*histSignal));
            RooHistPdf *signalPdf = new RooHistPdf("signalPdf", "Signal PDF", RooArgSet(x), signalRoo);
            
            // --- 3. Build the background RooFit objects.
            // We will keep the individual background PDFs and their coefficients separated.
            RooArgList bgPdfList;    // List of individual background PDFs
            RooArgList bgCoeffList;  // Corresponding list of background coefficients
            
            // Loop over background charges (from charge+1 to charge+n).
            std::vector<unsigned int> allowedZ = {charge-1,charge+1,charge+2};
            for (auto j : allowedZ) {
                // Check that background templates for this charge exist.
                if (templatesByCharge.find(j) == templatesByCharge.end())
                    continue;
                const auto &temp = templatesByCharge.at(j);
                if (temp.find("L2") == temp.end() || temp.at("L2").empty())
                    continue;
                
                // Select the appropriate histogram (i-th template if available; otherwise the last one).
                //TH1D* histBkg = (i < temp.at("L1").size()) ? temp.at("L1").at(i) : temp.at("L1").back();
                TH1D* histBkg = nullptr;
                bool aboveThreshold_bkg = (i >= useTemplateIndex);

                std::string key;
                if (charge == 14 || charge == 15 || charge == 16) { // example: Si, P, S
                    key = "L1";
                } else if ( charge == 18 || charge == 20) { //Ar, Ca
                    key = "L2";
                } else {
                    key = "L1"; // default/fallback
                }

                if (aboveThreshold_bkg) {
                    if (useTemplateIndex >= 0 && static_cast<size_t>(useTemplateIndex) < temp.at(key).size()) {
                        histBkg = temp.at(key).at(useTemplateIndex);
                    } else {
                        histBkg = temp.at(key).back();
                    }
                } else {
                    if (i < temp.at(key).size()) histBkg = temp.at(key).at(i);
                    else histBkg = temp.at(key).back();
                }
                bkgVec.push_back(histBkg);

                // Normalizza manualmente il template del background
                //histBkg->Scale(1.0 / histBkg->Integral());

                // Remove zero bins
                histBkg->Smooth(smt);
                histBkg = removeZeroBins(histBkg);
                if (histBkg->GetEntries() == 0 || histBkg->Integral() == 0) {
                    std::cerr << "Skipping empty background histogram for Z=" << j << ", bin=" << i << "\n";
                    continue;
                }
                
                // Build the RooDataHist and then the RooHistPdf for this background component.
                RooDataHist* bgData = new RooDataHist(Form("bgData_%u_%zu", j, i), "Background Data", RooArgList(x), Import(*histBkg));
                RooHistPdf* bgPdf  = new RooHistPdf(Form("bgPdf_%u_%zu", j, i), "Background PDF", RooArgSet(x), *bgData);
                
                // Add this background PDF to the list.
                bgPdfList.add(*bgPdf);
                
                // Initialize the coefficient for this background component from the histogram integral.
                double initVal = histData->Integral()/histBkg->Integral();
                RooRealVar* coeff = new RooRealVar(Form("coeffBg_%u_%zu", j, i),
                                                    "Background Coefficient", initialFractionValue(j,charge,Rbins_purity[i],Rbins_purity), 
                                            0., 10e+4*initialFractionValue(j,charge,Rbins_purity[i],Rbins_purity));
                bgCoeffList.add(*coeff);
            }
            
            // Check that at least one background component was found.
            if (bgPdfList.getSize() == 0) {
                std::cerr << "No background components found for bin " << i << std::endl;
                continue;
            }
            
            // --- 4. Set up coefficient RooFit variables for the overall signal and for each background component.
            double sigIntegral = histSignal->Integral();
            if (sigIntegral == 0 || !std::isfinite(sigIntegral)) {
                std::cerr << "Signal histogram has zero or non-finite integral\n";
                continue;
            }
            //double initSig = histData->Integral() / sigIntegral;
            double initSig = histData->Integral()/histSignal->Integral();
            //initSig = std::min(initSig, 100.0); // Prevent very large initial guesses

            RooRealVar *coeffSig = new RooRealVar("coeffSig", "Signal Coefficient",initialFractionValue(charge,charge,Rbins_purity[i],Rbins_purity), 
                                            0., 10e+4*initialFractionValue(charge,charge,Rbins_purity[i],Rbins_purity));

            // --- 5. Build the overall composite model as the sum of the signal PDF and each background PDF,
            // each multiplied by its own coefficient.
            RooArgList modelList;
            RooArgList modelCoeffs;
            
            // Add the signal PDF and its coefficient.
            modelList.add(*signalPdf);
            modelCoeffs.add(*coeffSig);
            
            // Loop over each background component and add it with its corresponding coefficient.
            for (int j = 0; j < bgPdfList.getSize(); j++) {
                // The j-th background PDF and its coefficient are already stored in the lists.
                modelList.add(*(dynamic_cast<RooAbsPdf*>(bgPdfList.at(j))));
                modelCoeffs.add(*(dynamic_cast<RooRealVar*>(bgCoeffList.at(j))));
            }
            
            // Build the composite model (sum of signal + individual backgrounds).
            RooAddPdf model("model", "Composite Model", modelList, modelCoeffs);
            
            // --- 6. Fit the model to the data.
            x.setRange("fitRange", chMin, chMax);
            RooFitResult* fitResult = model.fitTo(dataRoo, SumW2Error(false) ,Save(true), Timer(true),
                                                        Range("fitRange"), PrintLevel(-1), Verbose(false) );
            if (!fitResult) {
                std::cerr << "Fit failed for bin " << i << std::endl;
                continue;
            }
            // (Optional) Retrieve and print fit parameters.
            std::cout << "Fit results for bin " << i << ":\n";
            coeffSig->Print("v");
            for (int j = 0; j < bgCoeffList.getSize(); j++) {
                dynamic_cast<RooRealVar*>(bgCoeffList.at(j))->Print("v");
            }

            // --- 7. Compute the integral (yield) for each component over the full x range.
            // (Assuming the PDFs are normalized over x's full range.)

            // (Optional) Plotting.
            x.setRange("defaultRange",charge-2.3, charge+1.7);
            RooPlot* frame = x.frame(Range("fitRange"));
            frame->GetYaxis()->SetTitle("Counts");
            frame->GetXaxis()->CenterTitle(true);
            gStyle->SetLabelFont(62,"XYZ");
            gStyle->SetTitleFont(62,"XYZ");
            dataRoo.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.4),Name("dataHist"));
            model.plotOn(frame,
                        LineColor(kGreen+1),
                        LineWidth(1.2),
                        FillColor(kGreen - 9),
                        FillStyle(1001),  // Solid fill
                        //DrawOption("F"),  // Fill area
                        Precision(1e-2),
                        Range("fitRange"), NormRange("fitRange"),
                        Name("model_filled"));
            model.plotOn(frame,
                        Components(*signalPdf),
                        LineColor(kAzure),
                        LineWidth(1.2),
                        FillColor(kAzure - 9),
                        FillStyle(1001),
                        //DrawOption("F"),
                        Precision(1e-2),
                        Range("fitRange"), NormRange("fitRange"),
                        Name("signal_filled"));
            // Loop over background components to plot them individually.
            for (int j = 0; j < bgPdfList.getSize(); j++) {
                RooAbsPdf* bg = dynamic_cast<RooAbsPdf*>(bgPdfList.at(j));
                model.plotOn(frame,
                            Components(*bg),
                            LineColor(kRed + j),
                            LineWidth(1.2),
                            FillColor(kRed + j - 9),  // Slightly lighter fill
                            FillStyle(1001),
                            //DrawOption("F"),
                            Precision(1e-2),
                            Range("fitRange"), NormRange("fitRange"),
                            Name(Form("bg%d_filled", j)));
            }
            model.plotOn(frame,LineColor(kGreen+1),
            LineWidth(1.2),
            Range("fitRange"), NormRange("fitRange"),Name("model"));
            // Compute fit quality
            int nFitParams = 1 + 3; //1 for signal and 3 for background
            if (fitResult) {
                nFitParams = fitResult->floatParsFinal().getSize();
            } else {
                nFitParams = 1 + bgPdfList.getSize(); // fallback
            }
            auto [chi2, ndf, chi2ndf] = computeChi2FromFrame(frame, "model", "dataHist", nFitParams); //ORIGINAL
            chi_vs_shift->SetPoint(s, currentShift, chi2);
            //std::tie(chi2, ndf, chi2ndf, chi2OverNdfErr, nValid) = computeChi2WithBootstrap(
                //x, histData, originalSignal, bkgVec, charge, i, Rbins_purity, nFitParams, 100);
            
            
            frame->SetTitle( Form("%s @ %s - #chi^{2}/ndf = %.2f/%d = %.2f - Fit Range [%.1f - %.1f]",nucleus.Data(),title.Data(),
                                    chi2,ndf,chi2ndf,charge-2.,charge+2. ) );
            frame->GetYaxis()->SetRangeUser(0.1, frame->GetMaximum()*10.);
        
            // if (currentShift==0) double fr = computeFraction(coeffSig, bgCoeffList, fitResult, den, bkg, h, h_err);
            cAll->Update();
            frame->SetName(Form("fit_bin%zu_shift%.2f", i, currentShift));
            frame->GetYaxis()->SetTitle("Counts");
            frame->GetXaxis()->CenterTitle(true);
            if (s!=nShiftSteps-1) frame->GetXaxis()->SetTitle("");
            frame->GetYaxis()->SetLabelSize(0.058);
            frame->GetYaxis()->SetLabelOffset(0.001);
            frame->GetXaxis()->SetLabelSize(0.08);
            frame->GetYaxis()->SetTitleSize(0.062);
            frame->GetYaxis()->SetTitleOffset(0.8);
            frame->GetXaxis()->SetTitleSize(0.085);
            frame->GetXaxis()->SetTitleOffset(0.4);
            frame->GetXaxis()->SetNdivisions(306);
            frame->GetYaxis()->SetRangeUser(0.1, frame->GetMaximum()*10.);
            gStyle->SetLabelFont(62,"XYZ");
            gStyle->SetTitleFont(62,"XYZ");
            frame->SetTitle("");
            //allFrames.push_back(frame);
            cAll->cd(s+1);
            frame->Draw();
            TLegend* leg = new TLegend(0.15, 0.71, 0.70, 0.96);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.074);
            leg->SetTextFont(62);
            leg->AddEntry((TObject*)nullptr, Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2,ndf,chi2ndf), "");
            leg->AddEntry((TObject*)nullptr, Form("Shift = %.3f c.u.", currentShift), "");
            leg->Draw();
            if (s==0) leg->AddEntry((TObject*)nullptr, Form("[%.1f,%.1f] GV",Rbins_purity[i],Rbins_purity[i+1]),"");


            if (chi2 < minChi2) {
                minChi2 = chi2;
                bestShiftIndex = s;
                bestShift = currentShift;
                // Salvo una copia dei coefficienti migliori
                bestCoeffSig = std::make_unique<RooRealVar>(*coeffSig);
                
                bestBgCoeffList = std::make_unique<RooArgList>(bgCoeffList);
                
                bestFrame = std::unique_ptr<RooPlot>(static_cast<RooPlot*>(frame->Clone("bestFrame")));
                
                bestLegend = std::unique_ptr<TLegend>(static_cast<TLegend*>(leg->Clone("bestLegend")));
                
                bestFrameValues = std::make_unique<frameValues>(frameValues{chi2, ndf, chi2ndf, i, currentShift, 0.});
                
                bestSingleFit = std::unique_ptr<RooFitResult>(static_cast<RooFitResult*>(fitResult->Clone( Form("bestFit_bin_%d",i)  )));
                
                bestModel = std::unique_ptr<RooAddPdf>(static_cast<RooAddPdf*>(model.Clone(Form("bestModel_bin_%d", i))));
                
                RooAbsReal *integral_sig = signalPdf->createIntegral(x, NormSet(x), Range("fitRange"));
                
                RooProduct sig_yield("sig", "sig_yield", RooArgList(*integral_sig, *bestCoeffSig));
                double value_sigIntegral = sig_yield.getVal(RooArgSet(x));
                sigIntegral = integral_sig->getVal();
                double sigErr = integral_sig->getVal()*bestCoeffSig->getError();

                backgroundIntegral.clear();
                BkgErr.clear();
                
                for (int j = 0; j < bgPdfList.getSize(); j++) {
                    auto* bgPdf = dynamic_cast<RooAbsPdf*>(bgPdfList.at(j));
                    auto* coeffPtr = dynamic_cast<RooRealVar*>(bestBgCoeffList->at(j));
                    if (!bgPdf || !coeffPtr) {
                        std::cerr << "Null pointer when casting background PDF or coefficient!" << std::endl;
                        continue;
                    }
                    RooAbsReal* integral_bkg = bgPdf->createIntegral(x, NormSet(x), Range("fitRange"));
                    if (!integral_bkg) {
                        std::cerr << "Failed to create integral for background PDF " << j << std::endl;
                        continue;
                    }
                    RooRealVar& coeff = *coeffPtr;
                    RooProduct background_yield("background_yield", "background_yield", RooArgList(*integral_bkg, coeff));
                    double value_bkgIntegral = background_yield.getVal(RooArgSet(x));
                    double errorBkg = background_yield.getPropagatedError(*fitResult, x);
                    backgroundIntegral.push_back(integral_bkg->getVal());
                    BkgErr.push_back(integral_bkg->getVal()*coeff.getError());
                }
                
                
                bestBackgroundIntegrals = backgroundIntegral; // Aggiungi questo (serve dichiarazione std::vector<double> bestBackgroundIntegrals;)
                
                bestSigIntegral = sigIntegral; // Salva anche questa
                
                bestBkgErr = BkgErr;
                
                bestSigErr = sigErr;
                
            }
        } //shift loop
        bestFits.emplace_back(static_cast<RooFitResult*>(bestSingleFit->Clone(Form("bestFit_bin%d", i))));


        c1->cd();
        chi_vs_shift->SetMarkerStyle(20);
        chi_vs_shift->Draw("AP");
        chi_vs_shift->SetTitle( Form("[%.1f,%.1f] GV;shift (c.u.); chi2",Rbins_purity[i],Rbins_purity[i+1])  );
        TF1 *parabola = new TF1("parabola", "pol2", shift_min, shift_max);
        auto fitParabola = chi_vs_shift->Fit(parabola, "SLR"); // Fit and save result

        // Get parameters of the parabola
        double a = parabola->GetParameter(0);
        double b = parabola->GetParameter(1);
        double c = parabola->GetParameter(2);

        // Compute the minimum of the parabola
        double xmin = -b / (2.0 * c);
        double ymin = parabola->Eval(xmin);

        // Compute the x values where chi2 = min + 1
        double deltaChi2 = 1.0;
        double y_target = ymin + deltaChi2;
        double xmin_err;

        // Solve c*(x - xmin)^2 = 1 => x = xmin ± sqrt(1/c)
        if (c <= 0) {
            std::cerr << "Fit curvature is non-positive — no minimum!" << std::endl;
            xmin_err=0;
        } else xmin_err = std::sqrt(deltaChi2 / c);
        c1->SaveAs(shiftPDF.Data());


        double Rcenter = 0.5 * (Rbins_purity[i] + Rbins_purity[i+1]);
        shift_vs_R->SetPoint(i, Rcenter, xmin);
        shift_vs_R->SetPointError(i, 0, xmin_err);

        delete parabola;

   //------------------------------Refit with the optimal shift given by the parabola fit-----------------------
        double currentShift = xmin;

            // Remove zero bins
            histSignal->Smooth(smt);
            histSignal = removeZeroBins(histSignal);

            TH1D* shiftedTemplate = shiftHistogramResample(histSignalOriginal, currentShift);
            TH1D* originalSignal = (TH1D*)histSignal->Clone( Form("OriginalSignal") );
            if (shiftedTemplate) histSignal = shiftedTemplate;

            if (histSignal->GetEntries() == 0) {
                std::cerr << "Signal histogram is empty!" << std::endl;
                continue;
            }
            TString title = histData->GetTitle();
            RooDataHist dataRoo("dataRoo", "Data", RooArgList(x), Import(*histData));
            RooDataHist signalRoo("signalRoo", "Signal", RooArgList(x), Import(*histSignal));
            RooHistPdf *signalPdf = new RooHistPdf("signalPdf", "Signal PDF", RooArgSet(x), signalRoo);
            
            // --- 3. Build the background RooFit objects.
            // We will keep the individual background PDFs and their coefficients separated.
            RooArgList bgPdfList;    // List of individual background PDFs
            RooArgList bgCoeffList;  // Corresponding list of background coefficients
            
            // Loop over background charges (from charge+1 to charge+n).
            std::vector<unsigned int> allowedZ = {charge-1,charge+1,charge+2};
            for (auto j : allowedZ) {
                // Check that background templates for this charge exist.
                if (templatesByCharge.find(j) == templatesByCharge.end())
                    continue;
                const auto &temp = templatesByCharge.at(j);
                if (temp.find("L1") == temp.end() || temp.at("L1").empty())
                    continue;
                
                // Select the appropriate histogram (i-th template if available; otherwise the last one).
                TH1D* histBkg = (i < temp.at("L1").size() ) ? temp.at("L1").at(i) : temp.at("L1").back();
                bkgVec.push_back(histBkg);

                // Normalizza manualmente il template del background
                //histBkg->Scale(1.0 / histBkg->Integral());

                // Remove zero bins
                histBkg->Smooth(smt);
                histBkg = removeZeroBins(histBkg);
                if (histBkg->GetEntries() == 0 || histBkg->Integral() == 0) {
                    std::cerr << "Skipping empty background histogram for Z=" << j << ", bin=" << i << "\n";
                    continue;
                }
                
                // Build the RooDataHist and then the RooHistPdf for this background component.
                RooDataHist* bgData = new RooDataHist(Form("bgData_%u_%zu", j, i), "Background Data", RooArgList(x), Import(*histBkg));
                RooHistPdf* bgPdf  = new RooHistPdf(Form("bgPdf_%u_%zu", j, i), "Background PDF", RooArgSet(x), *bgData);
                
                // Add this background PDF to the list.
                bgPdfList.add(*bgPdf);
                
                // Initialize the coefficient for this background component from the histogram integral.
                double initVal = histData->Integral()/histBkg->Integral();
                RooRealVar* coeff = new RooRealVar(Form("coeffBg_%u_%zu", j, i),
                                                    "Background Coefficient", initialFractionValue(j,charge,Rbins_purity[i],Rbins_purity), 
                                            0., 10e+4*initialFractionValue(j,charge,Rbins_purity[i],Rbins_purity));
                bgCoeffList.add(*coeff);
            }
            
            // Check that at least one background component was found.
            if (bgPdfList.getSize() == 0) {
                std::cerr << "No background components found for bin " << i << std::endl;
                continue;
            }
            
            // --- 4. Set up coefficient RooFit variables for the overall signal and for each background component.
            double sigIntegral = histSignal->Integral();
            if (sigIntegral == 0 || !std::isfinite(sigIntegral)) {
                std::cerr << "Signal histogram has zero or non-finite integral\n";
                continue;
            }
            //double initSig = histData->Integral() / sigIntegral;
            double initSig = histData->Integral()/histSignal->Integral();
            //initSig = std::min(initSig, 100.0); // Prevent very large initial guesses

            RooRealVar *coeffSig = new RooRealVar("coeffSig", "Signal Coefficient", histSignal->Integral(), 
            0., 10e+4);

            // --- 5. Build the overall composite model as the sum of the signal PDF and each background PDF,
            // each multiplied by its own coefficient.
            RooArgList modelList;
            RooArgList modelCoeffs;
            
            // Add the signal PDF and its coefficient.
            modelList.add(*signalPdf);
            modelCoeffs.add(*coeffSig);
            
            // Loop over each background component and add it with its corresponding coefficient.
            for (int j = 0; j < bgPdfList.getSize(); j++) {
                // The j-th background PDF and its coefficient are already stored in the lists.
                modelList.add(*(dynamic_cast<RooAbsPdf*>(bgPdfList.at(j))));
                modelCoeffs.add(*(dynamic_cast<RooRealVar*>(bgCoeffList.at(j))));
            }
            
            // Build the composite model (sum of signal + individual backgrounds).
            RooAddPdf model("model", "Composite Model", modelList, modelCoeffs);
            
            // --- 6. Fit the model to the data.
            x.setRange("fitRange", charge-1.3, charge+0.9);
            RooFitResult* fitResult = model.fitTo(dataRoo, SumW2Error(false) ,Save(true), Timer(true),
                                                        Range("fitRange"), PrintLevel(-1), Verbose(false) );
            if (!fitResult) {
                std::cerr << "Fit failed for bin " << i << std::endl;
                continue;
            }
            // (Optional) Retrieve and print fit parameters.
            std::cout << "Fit results for bin " << i << ":\n";
            coeffSig->Print("v");
            for (int j = 0; j < bgCoeffList.getSize(); j++) {
                dynamic_cast<RooRealVar*>(bgCoeffList.at(j))->Print("v");
            }

            // --- 7. Compute the integral (yield) for each component over the full x range.
            // (Assuming the PDFs are normalized over x's full range.)

            // (Optional) Plotting.
            x.setRange("defaultRange",charge-2.3, charge+1.7);
            RooPlot* frame = x.frame();
            frame->GetYaxis()->SetTitle("Counts");
            frame->GetXaxis()->CenterTitle(true);
            gStyle->SetLabelFont(62,"XYZ");
            gStyle->SetTitleFont(62,"XYZ");
            dataRoo.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.4),Name("dataHist"));
            model.plotOn(frame,
                        LineColor(kGreen+1),
                        LineWidth(1.2),
                        FillColor(kGreen - 9),
                        FillStyle(1001),  // Solid fill
                        //DrawOption("F"),  // Fill area
                        Precision(1e-2),
                        Range("defaultRange"), NormRange("fitRange"),
                        Name("model_filled"));
            model.plotOn(frame,
                        Components(*signalPdf),
                        LineColor(kAzure),
                        //LineWidth(1.2),
                        FillColor(kAzure - 9),
                        FillStyle(1001),
                        //DrawOption("F"),
                        Precision(1e-2),
                        Range("defaultRange"), NormRange("fitRange"),
                        Name("signal_filled"));
            // Loop over background components to plot them individually.
            for (int j = 0; j < bgPdfList.getSize(); j++) {
                RooAbsPdf* bg = dynamic_cast<RooAbsPdf*>(bgPdfList.at(j));
                model.plotOn(frame,
                            Components(*bg),
                            LineColor(kRed + 1.2*j),
                            //LineWidth(1.2),
                            FillColor(kRed + j - 9),  // Slightly lighter fill
                            FillStyle(1001),
                            //DrawOption("F"),
                            Precision(1e-2),
                            Range("defaultRange"), NormRange("fitRange"),
                            Name(Form("bg%d_filled", j)));
            }
            model.plotOn(frame,LineColor(kGreen+1),
            LineWidth(1.2),
            Range("defaultRange"), NormRange("fitRange"),Name("model"));
            // Compute fit quality
            int nFitParams = 1 + 3; //1 for signal and 3 for background
            //double chi2, chi2ndf, chi2OverNdfErr;
            //int ndf, nValid;
            auto [chi2, ndf, chi2ndf] = computeChi2FromFrame(frame, "model", "dataHist", nFitParams); //ORIGINAL
            //std::tie(chi2, ndf, chi2ndf, chi2OverNdfErr, nValid) = computeChi2WithBootstrap(
                //x, histData, originalSignal, bkgVec, charge, i, Rbins_purity, nFitParams, 100);
            
            
            frame->SetTitle( Form("%s @ %s - #chi^{2}/ndf = %.2f/%d = %.2f - Fit Range [%.1f - %.1f]",nucleus.Data(),title.Data(),
                                    chi2,ndf,chi2ndf,charge-2.,charge+2. ) );
            frame->GetYaxis()->SetRangeUser(0.1, frame->GetMaximum()*10.);
        
            // if (currentShift==0) double fr = computeFraction(coeffSig, bgCoeffList, fitResult, den, bkg, h, h_err);
            cAll->Update();
            frame->SetName(Form("fit_bin%zu_shift%.2f", i, currentShift));
            frame->GetYaxis()->SetTitle("Counts");
            frame->GetXaxis()->CenterTitle(true);
            //if (s!=nShiftSteps-1) frame->GetXaxis()->SetTitle("");
            frame->GetYaxis()->SetLabelSize(0.058);
            frame->GetYaxis()->SetLabelOffset(0.001);
            frame->GetXaxis()->SetLabelSize(0.08);
            frame->GetYaxis()->SetTitleSize(0.062);
            frame->GetYaxis()->SetTitleOffset(0.8);
            frame->GetXaxis()->SetTitleSize(0.085);
            frame->GetXaxis()->SetTitleOffset(0.4);
            frame->GetXaxis()->SetNdivisions(306);
            frame->GetYaxis()->SetRangeUser(0.1, frame->GetMaximum()*10.);
            gStyle->SetLabelFont(62,"XYZ");
            gStyle->SetTitleFont(62,"XYZ");
            frame->SetTitle("");
            //allFrames.push_back(frame);
            //cAll->cd(s+1);
            //frame->Draw();
            TLegend* leg = new TLegend(0.15, 0.71, 0.70, 0.96);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.074);
            leg->SetTextFont(62);
            leg->AddEntry((TObject*)nullptr, Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2,ndf,chi2ndf), "");
            leg->AddEntry((TObject*)nullptr, Form("Shift = %.3f ± %.3f c.u.", currentShift, xmin_err), "");
            //leg->Draw();
            //if (s==0) leg->AddEntry((TObject*)nullptr, Form("[%.1f,%.1f] GV",Rbins_purity[i],Rbins_purity[i+1]),"");
            bestShift = currentShift;
                // Salvo una copia dei coefficienti migliori
                bestCoeffSig = std::make_unique<RooRealVar>(*coeffSig);
                
                bestBgCoeffList = std::make_unique<RooArgList>(bgCoeffList);
                
                bestFrame = std::unique_ptr<RooPlot>(static_cast<RooPlot*>(frame->Clone("bestFrame")));
                
                bestLegend = std::unique_ptr<TLegend>(static_cast<TLegend*>(leg->Clone("bestLegend")));
                
                bestFrameValues = std::make_unique<frameValues>(frameValues{chi2, ndf, chi2ndf, i, currentShift, xmin_err});
                
                bestSingleFit = std::unique_ptr<RooFitResult>(static_cast<RooFitResult*>(fitResult->Clone( Form("bestFit_bin_%d",i)  )));
                
                bestModel = std::unique_ptr<RooAddPdf>(static_cast<RooAddPdf*>(model.Clone(Form("bestModel_bin_%d", i))));
                
                RooAbsReal *integral_sig = signalPdf->createIntegral(x, NormSet(x), Range("fitRange"));
                
                RooProduct sig_yield("sig", "sig_yield", RooArgList(*integral_sig, *bestCoeffSig));
                double value_sigIntegral = sig_yield.getVal(RooArgSet(x));
                sigIntegral = integral_sig->getVal();
                double sigErr = bestCoeffSig->getError();

                backgroundIntegral.clear();
                BkgErr.clear();
                
                for (int j = 0; j < bgPdfList.getSize(); j++) {
                    auto* bgPdf = dynamic_cast<RooAbsPdf*>(bgPdfList.at(j));
                    auto* coeffPtr = dynamic_cast<RooRealVar*>(bestBgCoeffList->at(j));
                    if (!bgPdf || !coeffPtr) {
                        std::cerr << "Null pointer when casting background PDF or coefficient!" << std::endl;
                        continue;
                    }
                    RooAbsReal* integral_bkg = bgPdf->createIntegral(x, NormSet(x), Range("fitRange"));
                    if (!integral_bkg) {
                        std::cerr << "Failed to create integral for background PDF " << j << std::endl;
                        continue;
                    }
                    RooRealVar& coeff = *coeffPtr;
                    RooProduct background_yield("background_yield", "background_yield", RooArgList(*integral_bkg, coeff));
                    double value_bkgIntegral = background_yield.getVal(RooArgSet(x));
                    double errorBkg = background_yield.getPropagatedError(*fitResult, x);
                    backgroundIntegral.push_back(integral_bkg->getVal());
                    BkgErr.push_back(coeff.getError());
                }
                
                
                bestBackgroundIntegrals = backgroundIntegral; // Aggiungi questo (serve dichiarazione std::vector<double> bestBackgroundIntegrals;)
                
                bestSigIntegral = sigIntegral; // Salva anche questa
                
                bestBkgErr = BkgErr;
                
                bestSigErr = sigErr;
   //------------------------------End of refit--------------------------------

        computeFraction(*bestCoeffSig, *bestBgCoeffList,den, bkg, h, h_err,
                         bestSigIntegral, bestBackgroundIntegrals,bestSigErr,bestBkgErr,fitResult,singleContribution,singleContribution_err);


    // ---------------- compute L1 charge-cut efficiency for this bin ----------------
        // 'originalSignal' is the pre-shift copy you created before applying the best shift.
        // use the cut limits you computed earlier: cutLim = cutLimits(charge);
        double L1_low = cutLim.first;
        double L1_high = cutLim.second;

        // choose which template to use for the efficiency:
        // prefer originalSignal if available (pre-shift copy), otherwise histSignalOriginal
        TH1D* effTemplate = nullptr;
        if (histSignal) effTemplate = histSignal;
        /*else if (histSignalOriginal) effTemplate = histSignalOriginal;
        else effTemplate = histSignal; // fallback*/

        auto [effL1, effErrL1] = computeL1ChargeEfficiencyFromPdf(signalPdf, x, dataRoo, L1_low, L1_high);

        // store results
        L1eff.push_back(effL1);
        L1effErr.push_back(effErrL1);

        // optional debug print
        std::cout << Form("Bin %zu : L1 eff = %.5f ± %.5f (using template %s)", i, effL1, effErrL1,
                        (effTemplate ? effTemplate->GetName() : TString("null").Data())) << std::endl;



        cAll->Print("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"+getIonName(charge)+
                    Form("/%s_shift.pdf",ion));
        outfile->WriteTObject(chi_vs_shift,Form("chi_vs_shift_bin%zu",i));
        if (i== nBins-1) cAll->Print("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"+getIonName(charge)+
                                    Form("/%s_shift.pdf]",ion));
        canvas.cd();
        bestModel->Print("v");  // or "t" for tree structure
        //auto pull = buildPullHist(histData,bestModel.get(),x);
        formatFrameSinglePage(bestFrame,*bestFrameValues,Rbins_purity, canvas, charge);
        canvas.Print(outPDF.Data());

        outTxt << h.at(i) << " " << h_err.at(i) << " "
            << Rbins_purity[i] << " " << Rbins_purity[i+1];

        for (size_t j = 0; j < singleContribution[i].size(); ++j) {
            outTxt << " " << singleContribution[i][j]
                << " " << singleContribution_err[i][j];
        }

        outTxt << std::endl;
    }//bin
    outTxt.close();
    c1->cd();
    c1->SetLogx();
    shift_vs_R->SetMarkerStyle(20);
    shift_vs_R->GetXaxis()->SetRangeUser(2.15,100);
    shift_vs_R->GetYaxis()->SetRangeUser(-0.1,0.1);
    shift_vs_R->Draw("AP");
    TF1 *bestShift_fit = new TF1("bestShift_fit", "(x <= 55.45) ? [0]*log(x) + [1] : [0]*log(55.45) + [1]", 2.15, 1000);
    shift_vs_R->Fit(bestShift_fit,"LR");
    bestShift_fit->Draw("SAME");
    c1->SaveAs(shiftPDF.Data());
    c1->SaveAs((shiftPDF+"]").Data());
    outfile->WriteTObject(shift_vs_R,"shift_vs_R");

    canvas.cd();
    canvas.SetLogx();
    canvas.SetLogy(0);
    canvas.Print(outPDF.Data());
    //TF1 *purityFit = new TF1("purityFit", "1.0 / (1.0 + [0]*pow(x, [1]))", 0.8, 300.);
    TF1 *purityFit = new TF1("purityFit", "[0] + [1]*log10(x) + [2]*pow(log10(x), 2)", 0.8, 80.);
    //purityFit->SetParameters(1e-2, 1.0);
    purityFit->SetRange(2.2,1000);
    purityFit->Draw();

    canvas.Print((outPDF+"]").Data());
    TString fileOut = "../Fragmentation/BelowL1/Fractions/contamination_"+nucleus+".root";
    TFile *output = new TFile(fileOut.Data(), "RECREATE");
    //print fit
    for (const auto& fit : bestFits) {
        std::cout << " ------------------------------------ " << std::endl;
        fit->Print("vm"); 
        fit->Print("a");
        const RooArgList params = fit->floatParsFinal();
        for (int i = 0; i < params.getSize(); ++i) {
            RooRealVar* var = dynamic_cast<RooRealVar*>(params.at(i));
            if (var) {
                std::cout << "Parameter: " << var->GetName()
                        << "\n  Value   = " << var->getVal()
                        << "\n  Error   = " << var->getError()
                        << "\n  Range   = [" << var->getMin() << ", " << var->getMax() << "]\n"
                        << std::endl;
            }
        }
        std::cout << "Covariance matrix : \n";
        fit->covarianceMatrix().Print();
        std::cout << "Correlation matrix:\n";
        fit->correlationMatrix().Print();

    }
    //print fraction
    for (int i=0; i<den.size(); i++) 
        std::cout << "fraction : " << bkg.at(i)/den.at(i) << std::endl;
    //print fraction err
    for (int i=0; i<den.size(); i++) 
        std::cout << "fraction err: " << h_err.at(i) << std::endl;
    //print numerator
    std::cout << " ----------- " << std::endl;
    for (int i=0; i<bkg.size(); i++)
        std::cout << "numerator : " << bkg.at(i) << std::endl;
    //print den
    std::cout << " ----------- " << std::endl;
    for (int i=0; i<den.size(); i++)
        std::cout << "denominator : " << den.at(i) << std::endl;
    output->WriteTObject(purityFit,"purityFit");
    output->Close();


    TString fileOutL1 = "../IonsSelected/"+getIonPath(charge)+"/Efficiencies/L1ChargeCut/l1ChargeCut.root";
    TFile *outputL1 = new TFile(fileOutL1.Data(), "RECREATE");
    // Build TH1D with the same Rbins_purity binning
    TH1D* hL1Eff = new TH1D("hL1Eff", Form("L1 Charge Cut Efficiency - Z=%s", ion),
                            (int)Rbins_purity.size()-1, Rbins_purity.data());
    hL1Eff->GetXaxis()->SetTitle("R (GV)");
    hL1Eff->GetYaxis()->SetTitle("L1 cut efficiency");

    // Fill histogram from the stored vectors
    for (size_t ib = 0; ib < L1eff.size(); ++ib) {
        auto low_edge = hL1Eff->GetBinLowEdge(ib+1);
        auto high_edge= hL1Eff->GetBinLowEdge(ib+1)+hL1Eff->GetBinWidth(ib+1);
        double x = get_x_from_flux(charge,low_edge,high_edge);
        hL1Eff->SetBinContent(hL1Eff->FindBin(x), L1eff[ib]);
        hL1Eff->SetBinError(hL1Eff->FindBin(x), (ib < L1effErr.size() ? L1effErr[ib] : 0.0));
    }

    // write to output file
    outputL1->cd();
    outputL1->WriteTObject(hL1Eff,"hL1Eff");
}
void fitHistograms(std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>>& templatesByCharge,
                   unsigned int charge, unsigned int n,TString option) {
    //binning based on charge
    std::vector<double> Rbins_purity;
    if (charge == 15) {
        Rbins_purity = {0.8, 3.64, 5.37, 7.76, 11.00, 15.3, 21.1, 38.9, 70.00, 150.00, 300.00, 500.00, 1000.00};
    } else {
        Rbins_purity = {
            2.15, 2.40, 2.67, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37,
            5.90, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0,
            14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8,
            33.5, 38.9, 45.1, 52.2, 60.3, 69.7, 80.5, 93.0, 108., 125.,
            147., 175., 211., 259., 330., 441., 660., 1200., 3000.
        };
    }
    int nRbins_purity = Rbins_purity.size();
    /*const int nP_rigBins = 13; 
	const double P_rigBins[nP_rigBins] = {0.8,3.64,5.37,7.76,11.00,15.3,21.1,38.9,70.00,150.00,300.00,500.00,1000.00};*/
    // --- 1. Check that the main charge exists and has the necessary templates.
    if (templatesByCharge.find(charge) == templatesByCharge.end()) {
        std::cerr << "Missing templates for charge " << charge << "\n";
        return;
    }
    const auto &mainTemplates = templatesByCharge.at(charge);
    if (mainTemplates.find("L1") == mainTemplates.end() || mainTemplates.at("L1").empty()) {
        std::cerr << "No L1 templates for charge " << charge << "\n";
        return;
    }
    if (mainTemplates.find("L2") == mainTemplates.end() || mainTemplates.at("L2").empty()) {
        std::cerr << "No L2 templates for charge " << charge << "\n";
        return;
    }
    
    // --- 2. Set up the templates for data and signal.
    auto dataTemplates   = templatesByCharge[charge]["L1"];
    auto signalTemplates = templatesByCharge[charge]["L2"];
    const int nBins = dataTemplates.size();

    // --- Set up the output
    const TString nucleus = getIonName(charge);
    const TString basePath = "../Fragmentation/BelowL1/Fractions/";
    const TString outPDF  = basePath + nucleus + "_purityFit.pdf";
    TCanvas canvas("canvasBin", "Bin", 800, 600);
    canvas.SetLogy();
    canvas.Print((outPDF+"[").Data());
    auto f_vs_r = new TH1D("f_vs_r","", nRbins_purity-1, Rbins_purity.data());
    f_vs_r->GetXaxis()->SetTitle("R (GV)");
    f_vs_r->GetYaxis()->SetTitle("Fraction");

    //debug
    std::vector<double> backgroundIntegral;
    std::vector<double> sig;
    // Loop over bins (each bin corresponds to one pair of data/signal templates)
    for (size_t i = 0; i < static_cast<size_t>(nBins); ++i) {
        // Define the observable for the fit.
        // (Adjust the range as needed.)
        RooRealVar x("x", "Q_{L1}", charge - 5, charge + 5);
        
        // --- Data and Signal RooFit objects.
        TH1D* histData   = dataTemplates.at(i);
        TH1D* histSignal = (i < signalTemplates.size()) ? signalTemplates.at(i) : signalTemplates.back();
        if (histSignal->GetEntries() == 0) {
            std::cerr << "Signal histogram is empty!" << std::endl;
            continue;
        }
        if (option == "y") {
            const double shift = -0.02;  // <-- adjust this value as needed
            TH1D* shifted = shiftHistogram(histSignal, shift);
            if (shifted) histSignal = shifted;
        }
        TString title = histData->GetTitle();
        RooDataHist dataRoo("dataRoo", "Data", x, Import(*histData));
        RooDataHist signalRoo("signalRoo", "Signal", x, Import(*histSignal));
        RooHistPdf signalPdf("signalPdf", "Signal PDF", x, signalRoo);
        // Define the asymmetric smearing model
        RooRealVar mean("mean", "mean", -0.1, -0.05, -0.2);
        RooRealVar sigma("sigma", "sigma", 0.05, 0.01, 0.1);  // Width of the central Gaussian
        RooRealVar alpha("alpha", "alpha", 6, 5, 10);        // Controls the tail on one side
        RooRealVar nn("n", "n", 0.5 , 0.1, 1);               // Controls the fall-off rate in the tail

        // Use RooCrystalBall for asymmetry
        RooCrystalBall asymGauss("asymGauss", "Asymmetrical Gaussian-like smearing", x, mean, sigma, alpha, nn);
        // Convolve the L2 PDFs with the asymmetrical function
        RooFFTConvPdf signalPdf_1("convL2_z1", "Convolution of L2_z1 with Asymmetrical Smearing", x, signalPdf_1, asymGauss);
        
        // --- 3. Build the background RooFit objects.
        // We will keep the individual background PDFs and their coefficients separated.
        RooArgList bgPdfList;    // List of individual background PDFs
        RooArgList bgCoeffList;  // Corresponding list of background coefficients
        
        // Loop over background charges (from charge+1 to charge+n).
        for (unsigned int j = charge + 1; j <= charge + n; j++) {
            // Check that background templates for this charge exist.
            if (templatesByCharge.find(j) == templatesByCharge.end())
                continue;
            const auto &temp = templatesByCharge.at(j);
            if (temp.find("L1") == temp.end() || temp.at("L1").empty())
                continue;
            
            // Select the appropriate histogram (i-th template if available; otherwise the last one).
            TH1D* histBkg = (i < temp.at("L1").size()) ? temp.at("L1").at(i) : temp.at("L1").back();
            
            // Build the RooDataHist and then the RooHistPdf for this background component.
            RooDataHist* bgData = new RooDataHist(Form("bgData_%u_%zu", j, i), "Background Data", x, Import(*histBkg));
            RooHistPdf* bgPdf  = new RooHistPdf(Form("bgPdf_%u_%zu", j, i), "Background PDF", x, *bgData);
            
            // Add this background PDF to the list.
            bgPdfList.add(*bgPdf);
            
            // Initialize the coefficient for this background component from the histogram integral.
            double initVal = histBkg->Integral();
            RooRealVar* coeff = new RooRealVar(Form("coeffBg_%u_%zu", j, i),
                                                "Background Coefficient", initVal, 0., 2 * initVal);
            // Background-specific ranges
            auto bkRange = getBkgRange(j);
            double bgRangeLow = bkRange.first;
            double bgRangeHigh= bkRange.second;
            /*double bgRangeLow = j-0.7;//This shoul be called by a function that based on the charge give me the correct range
            double bgRangeHigh = j+0.9;*/
            double dataIntegralBg = histData->Integral(histData->FindBin(bgRangeLow), histData->FindBin(bgRangeHigh));
            double templateIntegralBg = histBkg->Integral(histBkg->FindBin(bgRangeLow), histBkg->FindBin(bgRangeHigh));
            /*if (templateIntegralBg == 0) {
                coeff->setVal(0);
                bgCoeffList.add(*coeff);
                std::cerr << "Background template integral is 0, avoiding division for charge " << j << std::endl;
                continue;
            }*/
            //coeff->setVal(templateIntegralBg /dataIntegralBg);
            bgCoeffList.add(*coeff);
        }
        
        // Check that at least one background component was found.
        if (bgPdfList.getSize() == 0) {
            std::cerr << "No background components found for bin " << i << std::endl;
            continue;
        }
        
        // --- 4. Set up coefficient RooFit variables for the overall signal and for each background component.
        double initSig = histSignal->Integral();
        RooRealVar coeffSig("coeffSig", "Signal Coefficient", initSig, 0, 2 * initSig);
        // Example: Set coefficients using data integrals over specific ranges
        double dataIntegralSignal = histData->Integral(histData->FindBin(getSigRange(charge).first), histData->FindBin(getSigRange(charge).second));
        double templateIntegralSignal = histSignal->Integral(histSignal->FindBin(getSigRange(charge).first), histSignal->FindBin(getSigRange(charge).second));
        /*if (templateIntegralSignal == 0) {
            std::cerr << "templateIntegralSignal is 0, avoiding division" << std::endl;
            continue;
        }*/
        //coeffSig.setVal(templateIntegralSignal /dataIntegralSignal );
        
        // (Each background coefficient is already defined and stored in bgCoeffList.)
        
        // --- 5. Build the overall composite model as the sum of the signal PDF and each background PDF,
        // each multiplied by its own coefficient.
        RooArgList modelList;
        RooArgList modelCoeffs;
        
        // Add the signal PDF and its coefficient.
        modelList.add(signalPdf);
        modelCoeffs.add(coeffSig);
        
        // Loop over each background component and add it with its corresponding coefficient.
        for (int j = 0; j < bgPdfList.getSize(); j++) {
            // The j-th background PDF and its coefficient are already stored in the lists.
            modelList.add(*(dynamic_cast<RooAbsPdf*>(bgPdfList.at(j))));
            modelCoeffs.add(*(dynamic_cast<RooRealVar*>(bgCoeffList.at(j))));
        }
        
        // Build the composite model (sum of signal + individual backgrounds).
        RooRealSumPdf model("model", "Composite Model", modelList, modelCoeffs, kTRUE);
        
        // --- 6. Fit the model to the data.
        x.setRange("integrationRange", 0, charge+2);
        RooFitResult* fitResult = model.fitTo(dataRoo, Extended(kTRUE), Save(), Timer(true) );
        if (!fitResult) {
            std::cerr << "Fit failed for bin " << i << std::endl;
            continue;
        }
        
        // (Optional) Retrieve and print fit parameters.
        std::cout << "Fit results for bin " << i << ":\n";
        coeffSig.Print("v");
        for (int j = 0; j < bgCoeffList.getSize(); j++) {
            dynamic_cast<RooRealVar*>(bgCoeffList.at(j))->Print("v");
        }

        // --- 7. Compute the integral (yield) for each component over the full x range.
        // (Assuming the PDFs are normalized over x's full range.)
        double signalIntegral = coeffSig.getVal() *
            signalPdf.createIntegral(x)->getVal();
        sig.push_back(signalIntegral);
        std::cout << "Bin " << i << " - Signal integral: " << signalIntegral << std::endl;
        
        for (int j = 0; j < bgPdfList.getSize(); j++) {
            RooAbsPdf* bkgPdf = dynamic_cast<RooAbsPdf*>(bgPdfList.at(j));
            RooRealVar* bkgCoeff = dynamic_cast<RooRealVar*>(bgCoeffList.at(j));
            double bkgIntegral = bkgCoeff->getVal() *
                bkgPdf->createIntegral(x)->getVal();
            std::cout << "Bin " << i << " - Background component " << j
                      << " integral: " << bkgIntegral << std::endl;
            backgroundIntegral.push_back(bkgIntegral);
        }

        // (Optional) Plotting.
        RooPlot* frame = x.frame();
        frame->GetYaxis()->SetTitle("Counts");
        gStyle->SetLabelFont(62,"XYZ");
        gStyle->SetTitleFont(62,"XYZ");
        dataRoo.plotOn(frame);
        model.plotOn(frame, LineColor(kGreen) );  // Full model     
        model.plotOn(frame, Components(signalPdf), LineColor(kBlue));  // Signal
        // Loop over background components to plot them individually.
        for (int j = 0; j < bgPdfList.getSize(); j++) {
            model.plotOn(frame, Components(*(dynamic_cast<RooAbsPdf*>(bgPdfList.at(j)))),
                         LineColor(kRed + 2*j));
        }
        model.plotOn(frame, LineColor(kGreen));
        frame->SetTitle( Form("%s @ %s - #chi^{2}=%f",nucleus.Data(),title.Data(),frame->chiSquare() ) );
        frame->Draw();
        frame->GetYaxis()->SetRangeUser(1., frame->GetMaximum()*1.);
        canvas.Print(outPDF.Data());
        /*double h = computeFraction(coeffSig, bgCoeffList, fitResult, signalIntegral,backgroundIntegral);
        if (h>0) {
            f_vs_r->SetBinContent(Rbins_HighZ[i+1], 1-h);
            f_vs_r->SetBinError(Rbins_HighZ[i+1],0.000001);
        }*/
        // Cleanup: delete the fit result.
        //delete fitResult;
    }
    canvas.SetLogx();
    canvas.SetLogy(0);
    f_vs_r->SetMarkerStyle(20);
    f_vs_r->Draw("P");
    canvas.Print(outPDF.Data());
    //auto spline_purity = autospline(f_vs_r,4,30);
    //spline_purity->Draw();
    canvas.Print(outPDF.Data());
    canvas.Print((outPDF+"]").Data());
    TString fileOut = "../Fragmentation/BelowL1/Fractions/contamination_"+nucleus+".root";
    TFile *output = new TFile(fileOut.Data(), "RECREATE");
    output->WriteTObject(f_vs_r,"f_vs_r");
    //print fraction
    for (int i=0; i<f_vs_r->GetNbinsX(); i++) 
        std::cout << f_vs_r->GetBinContent(i) << std::endl;
    //print numerator
    std::cout << " ----------- " << std::endl;
    for (int i=0; i<backgroundIntegral.size(); i++ )
        std::cout << backgroundIntegral[i] << std::endl;
    //print sigVal
    std::cout << " ----------- " << std::endl;
    for (int i=0; i<sig.size(); i++)
        std::cout << sig[i] << std::endl;
    //output->WriteTObject(spline_purity,"spline_purity");
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
double computeFraction(const RooRealVar &coeffSig, const RooArgList &bgCoeffList,
                       std::vector<double> &den, std::vector<double> &bkg,
                       std::vector<double> &h, std::vector<double> &h_err,
                       double sigIntegral, std::vector<double> bkgIntegral,
                       double err_sig, std::vector<double> err_bkg,
                       RooFitResult* fitRes,
                       std::vector<std::vector<double>> &singleContribution,
                       std::vector<std::vector<double>> &singleContributionErr,
                       int uncertaintyOption, // 1 = fit error only, 2 = template stats, 3 = quadrature
                       int nToys) {
    // --- 1) Compute N and B from fit coefficients * integrals
    double nP = coeffSig.getVal();
    double N  = nP * sigIntegral;

    double B = 0;
    for (int i = 0; i < bgCoeffList.getSize(); ++i) {
        const RooRealVar& bi = static_cast<const RooRealVar&>(bgCoeffList[i]);
        B += bi.getVal() * bkgIntegral[i];
    }

    double Den = N + B;
    double fr  = (Den > 0) ? B / Den : 0.0;

    den.push_back(Den);
    bkg.push_back(B);
    h.push_back(fr);

    // --- 2) Collect parameter names and integrals
    std::vector<std::string> names;
    std::vector<double> Ivals;
    names.push_back(coeffSig.GetName());
    Ivals.push_back(sigIntegral);
    for (int i = 0; i < bgCoeffList.getSize(); ++i) {
        const RooRealVar& bi = static_cast<const RooRealVar&>(bgCoeffList[i]);
        names.push_back(bi.GetName());
        Ivals.push_back(bkgIntegral[i]);
    }

    // --- 3) Derivatives of total background fraction (df/dn)
    std::map<std::string, double> dfdn;
    for (size_t i = 0; i < names.size(); ++i) {
        if (names[i] == coeffSig.GetName())
            dfdn[names[i]] = -Ivals[i] * B / (Den * Den);
        else
            dfdn[names[i]] = Ivals[i] * N / (Den * Den);
    }

    // --- 4) Covariance-based uncertainty propagation for total fraction
    double sigma2_fit = 0;
    if ((uncertaintyOption == 1 || uncertaintyOption == 3) && fitRes) {
        const RooArgList& params = fitRes->floatParsFinal();
        const TMatrixDSym& covmat = fitRes->covarianceMatrix();

        for (size_t i = 0; i < names.size(); ++i) {
            int ii = params.index(params.find(names[i].c_str()));
            if (ii < 0) continue;
            for (size_t j = 0; j < names.size(); ++j) {
                int jj = params.index(params.find(names[j].c_str()));
                if (jj < 0) continue;
                sigma2_fit += dfdn[names[i]] * dfdn[names[j]] * covmat[ii][jj];
            }
        }
    }

    // --- 5) Template-statistical errors for total fraction
    double sigma2_stat = 0;
    if (uncertaintyOption == 2 || uncertaintyOption == 3) {
        sigma2_stat += std::pow(sigIntegral * err_sig * B / (Den * Den), 2);
        for (size_t i = 0; i < err_bkg.size(); ++i) {
            double I_i = bkgIntegral[i];
            sigma2_stat += std::pow(I_i * err_bkg[i] * N / (Den * Den), 2);
        }
    }

    double sigma2_total = 0;
    if      (uncertaintyOption == 1) sigma2_total = sigma2_fit;
    else if (uncertaintyOption == 2) sigma2_total = sigma2_stat;
    else if (uncertaintyOption == 3) sigma2_total = sigma2_fit + sigma2_stat;

    double sigma_fr = (sigma2_total > 0) ? std::sqrt(sigma2_total) : 0.;
    h_err.push_back(sigma_fr);

    // --- 6) Compute individual background fractions and uncertainties
    std::vector<double> thisFrac;
    std::vector<double> thisFracErr;

    const RooArgList& params = (fitRes ? fitRes->floatParsFinal() : RooArgList());
    const TMatrixDSym& covmat = (fitRes ? fitRes->covarianceMatrix() : TMatrixDSym());

    for (int k = 0; k < bgCoeffList.getSize(); ++k) {
        const RooRealVar& bk = static_cast<const RooRealVar&>(bgCoeffList[k]);
        double n_k = bk.getVal();
        double I_k = bkgIntegral[k];
        double f_k = (Den > 0) ? n_k * I_k / Den : 0.0;
        thisFrac.push_back(f_k);

        double sigma2_fk = 0;

        // --- Derivatives wrt each parameter n_j
        std::map<std::string, double> dfdn_k;
        for (size_t j = 0; j < names.size(); ++j) {
            const std::string& nj = names[j];
            double Ij = Ivals[j];
            if (nj == coeffSig.GetName()) {
                dfdn_k[nj] = -Ivals[0] * n_k * I_k / (Den * Den);
            } else if (nj == bk.GetName()) {
                dfdn_k[nj] = (Ij * Den - n_k * I_k * Ij) / (Den * Den);
            } else {
                dfdn_k[nj] = -n_k * I_k * Ij / (Den * Den);
            }
        }

        // --- Fit uncertainty via covariance
        if ((uncertaintyOption == 1 || uncertaintyOption == 3) && fitRes) {
            for (size_t i = 0; i < names.size(); ++i) {
                int ii = params.index(params.find(names[i].c_str()));
                if (ii < 0) continue;
                for (size_t j = 0; j < names.size(); ++j) {
                    int jj = params.index(params.find(names[j].c_str()));
                    if (jj < 0) continue;
                    sigma2_fk += dfdn_k[names[i]] * dfdn_k[names[j]] * covmat[ii][jj];
                }
            }
        }

        // --- Template statistical uncertainty (quadrature)
        if (uncertaintyOption == 2 || uncertaintyOption == 3) {
            double dfdnP = -sigIntegral * n_k * I_k / (Den * Den);
            sigma2_fk += std::pow(dfdnP * err_sig, 2);
            for (size_t j = 0; j < err_bkg.size(); ++j) {
                double Ij = bkgIntegral[j];
                double dfdnj = (j == (size_t)k)
                                   ? (Ij * Den - n_k * I_k * Ij) / (Den * Den)
                                   : -n_k * I_k * Ij / (Den * Den);
                sigma2_fk += std::pow(dfdnj * err_bkg[j], 2);
            }
        }

        thisFracErr.push_back((sigma2_fk > 0) ? std::sqrt(sigma2_fk) : 0.);
    }

    singleContribution.push_back(thisFrac);
    singleContributionErr.push_back(thisFracErr);

    return fr;
}
void fitHistograms(std::pair<std::vector<TH1D*>, std::vector<TH1D*>> z1, std::pair<std::vector<TH1D*>, std::vector<TH1D*>> z2, unsigned int charge) {
    // Get the correct binning
    int binL1z1, binL2z1, binL1z2, binL2z2;
    binL1z1 = z1.first.size();  //Z1, L1
    binL2z1 = z1.second.size(); //Z1, L2
    binL1z2 = z2.first.size();  //Z2, L1
    binL2z2 = z2.second.size(); //Z2, L2

    //vector to store the fractions:
    std::vector<double> frac;
    std::vector<double> frac_err;
    std::vector<double> sp_err;
    std::vector<double> par_b_err, par_c_err, par_b_val, par_c_val;
    std::vector<TString> par_b_stat, par_c_stat;
    std::vector<FitStatus> fitStatusList;

    TString path1 = "../Fragmentation/BelowL1/Fractions/";
    TString nucleus = getIonName(charge);
    TString out1 = path1+nucleus+"_purityFit.pdf[";
    TString out2 = path1+nucleus+"_purityFit.pdf";
    TString out3 = path1+nucleus+"_purityFit.pdf]";
    TCanvas *canvasBin = new TCanvas("canvasBin", "Bin", 800, 600);
    canvasBin->SetLogy();
    canvasBin->cd();
    canvasBin->Print(out1.Data() );
    auto f_vs_r = (TH1D*)hist_rig_highZ->Clone();
    f_vs_r->GetXaxis()->SetTitle("R (GV)");
    f_vs_r->GetYaxis()->SetTitle("Fraction");
    auto range = getRange(charge);

    for (int i=0; i<binL1z1; i++) { 
        // Define the observable variable (Q)
        RooRealVar x("x", "Q", range.first , range.second );

        // Create RooDataHist from TH1D ---- first=L1, second=L2
        TH1D* hist1 = z1.first[i];  // z1, L1 --- To be fitted
        TH1D* hist2 = nullptr;
        TH1D* hist3 = nullptr;
        if (i >= binL2z1 || i >= binL2z2) { //Signal from L2, background from L1
        hist2 = z1.second[binL2z1-1]; // z1, L2
        hist3 = z2.first[binL2z2-1]; // z2, L1
        } else {
        hist2 = z1.second[i]; // z1, L2
        hist3 = z2.first[i]; // z2, L1 
        }

        //hist2->Scale( 1./hist2->Integral(),"WIDTH");
        //hist3->Scale( 1./hist3->Integral(),"WIDTH");

        // Create RooDataHist objects
        RooDataHist dataHistL1_z1("dataHistL1_z1", "dataHistL1_z1", x, Import(*hist1) ); //z1.first[i]
        RooDataHist dataHistL2_z1("dataHistL2_z1", "dataHistL2_z1", x, Import(*hist2) ); //z1.second[i]
        RooDataHist dataHistL2_z2("dataHistL2_z2", "dataHistL2_z2", x, Import(*hist3) ); //z2.first[i]

        // Create RooHistPdf from RooDataHist for L2
        RooHistPdf pdfL2_z1("pdfL2_z1", "pdfL2_z1", x, dataHistL2_z1 );
        RooHistPdf pdfL2_z2("pdfL2_z2", "pdfL2_z2", x, dataHistL2_z2 );

        // Define the asymmetric smearing model
        RooRealVar mean("mean", "mean", -0.1, -0.05, -0.2);
        RooRealVar sigma("sigma", "sigma", 0.05, 0.01, 0.1);  // Width of the central Gaussian
        RooRealVar alpha("alpha", "alpha", 6, 5, 10);        // Controls the tail on one side
        RooRealVar n("n", "n", 0.5 , 0.1, 1);               // Controls the fall-off rate in the tail

        // Use RooCrystalBall for asymmetry
        RooCrystalBall asymGauss("asymGauss", "Asymmetrical Gaussian-like smearing", x, mean, sigma, alpha, n);

        // Convolve the L2 PDFs with the asymmetrical function
        RooFFTConvPdf convL2_z1("convL2_z1", "Convolution of L2_z1 with Asymmetrical Smearing", x, pdfL2_z1, asymGauss);

        //RooFFTConvPdf convL2_z2("convL2_z2", "Convolution of L2_z2 with Gaussian", x, pdfL2_z2, gauss);

        x.setRange("integrationRange", charge-0.57, charge+0.89);
        RooAbsReal* integral = convL2_z1.createIntegral(x, RooFit::Range("integrationRange"));
        RooAbsReal* conv_integral_den = convL2_z1.createIntegral(x, RooFit::Range(0, charge + 0.89));
        double integralValue = integral->getVal();

        // Normalize the coefficients based on the number of entries in the histograms
        double coeffL1_z1,coeffL1_z2;
        coeffL1_z1 = hist2->Integral(0, hist2->FindBin(charge + 0.89)); //OK, ORIGINAL
        //coeffL1_z1=integralValue;
        coeffL1_z2 = 5*hist3->Integral(0, hist3->FindBin(charge + 0.89)-coeffL1_z1); //OK

        RooRealVar coeffz1("coeffz1", "Coefficient for Z1", coeffL1_z1, coeffL1_z1, coeffL1_z1);
        RooRealVar coeffz2("coeffz2", "Coefficient for Z2", coeffL1_z2, coeffL1_z2, coeffL1_z2);

        // Create the model as a sum of the convolved PDFs
        RooRealSumPdf model("model", "model", RooArgList(pdfL2_z1, pdfL2_z2), RooArgList(coeffz1, coeffz2));

        // Fit the model to the L1 data
        RooFitResult* r = model.fitTo(dataHistL1_z1, 
                                    RooFit::Timer(true), 
                                    RooFit::Save(), 
                                    RooFit::NumCPU(8,kTRUE),
                                    RooFit::PrintLevel(-1));

        // Set the range for visualization
        x.setRange("range", 0, charge + 0.8);

        RooPlot* frame = x.frame();
        frame->GetYaxis()->SetTitle("Events");
        TString title = hist1->GetTitle();
        dataHistL1_z1.plotOn(frame);
        model.plotOn(frame,LineColor(kBlue));
        model.plotOn(frame,Components(pdfL2_z2), LineColor(kRed));
        model.plotOn(frame,LineColor(kGreen));    
        frame->SetTitle( Form("%s @ %s - #chi^{2}=%f",nucleus.Data(),title.Data(),frame->chiSquare() ) );
        frame->Draw();

        // Check if the fit result is valid
        if (r) {
            int fitStatus = r->status();      // Status of the fit (0 = OK, non-zero = issues)
            int covQuality = r->covQual();    // Covariance matrix quality (0, 1, 2, or 3)
            double chiSquare = frame->chiSquare(); // Chi-square value
            // Store in the vector
            fitStatusList.emplace_back(fitStatus, covQuality,chiSquare);
        } else {
            std::cerr << "Fit result 'r' is null!" << std::endl;
        }
        canvasBin->Print(out2.Data() );


        //integral della RooAddPdf
        //auto pdf= model.pdf("pdf");
        //x.setRange("range",0,charge+0.8);
        //auto integral_z1 = model.createIntegral(x,RooFit::NormSet(x), RooFit::Range("range"));
        //auto integral_z2 = pdfL2_z2.createIntegral(x,RooFit::NormSet(x),RooFit::Range("range"));

        const RooArgList& params = r->floatParsFinal();
        RooRealVar* param_b = (RooRealVar*)params.find("coeffz1");
        RooRealVar* param_c = (RooRealVar*)params.find("coeffz2");
        double param_b_value = param_b->getVal();
        double param_b_error = param_b->getError();
        double param_c_value = param_c->getVal();
        double param_c_error = param_c->getError();
        double den, num, error_num, error_den;


        if (param_b->isConstant()) {
            par_b_stat.push_back("param_b (coeffz1) is fixed and not floating in the fit.");
        } else {
            par_b_stat.push_back("param_b (coeffz1) is floating in the fit.");
        }

        if (param_c->isConstant()) {
            par_c_stat.push_back("param_c (coeffz2) is fixed and not floating in the fit.");
        } else {
            par_c_stat.push_back("param_c (coeffz2) is floating in the fit.");
        }


        //if (i >= binL2z1 || i >= binL2z2) { 

            //CONVOLUTUION//
            /*den = param_b_value*integralValue + 
                param_c_value*z2.first.at(binL1z2-1)->Integral(0, z2.first.at(binL1z2-1)->FindBin(charge+0.8));*/

            //NO CONVOLUTION//
            den = param_b_value*z1.second.at(binL2z1-1)->Integral(hist2->FindBin(charge -0.57),z1.second.at(binL2z1-1)->FindBin(charge+0.89)) + 
                param_c_value*z2.first.at(binL1z2-1)->Integral(hist2->FindBin(charge -0.57), z2.first.at(binL1z2-1)->FindBin(charge+0.89));

            num = param_c_value*z2.first.at(binL1z2-1)->Integral(hist2->FindBin(charge -0.57), z2.first.at(binL1z2-1)->FindBin(charge+0.89));

            error_num = param_c_value*z2.first.at(binL1z2-1)->GetBinError(hist2->FindBin(charge -0.57), z2.first.at(binL1z2-1)->FindBin(charge+0.89)) + 
                            param_c_error*z2.first.at(binL1z2-1)->Integral(hist2->FindBin(charge -0.57), z2.first.at(binL1z2-1)->FindBin(charge+0.89));

            error_den = param_b_value*z1.second.at(binL2z1-1)->GetBinError(hist2->FindBin(charge -0.57), z1.second.at(binL2z1-1)->FindBin(charge+0.89)) + 
                        param_b_error*z1.second.at(binL2z2-1)->Integral(hist2->FindBin(charge -0.57), z1.second.at(binL2z2-1)->FindBin(charge+0.89)) +
                            param_c_value*z2.first.at(binL1z2-1)->GetBinError(hist2->FindBin(charge -0.57), z2.first.at(binL1z2-1)->FindBin(charge+0.89)) + 
                            param_c_error*z2.first.at(binL1z2-1)->Integral(hist2->FindBin(charge -0.57), z2.first.at(binL1z2-1)->FindBin(charge+0.89));
        /*} else {
            den = param_b_value*z1.first.at(i)->Integral(0, z1.first.at(i)->FindBin(charge+0.8)) + 
                param_c_value*z2.first.at(i)->Integral(0, z2.first.at(i)->FindBin(charge+0.8));
            num = param_c_value*z2.second.at(i)->Integral(0, z2.second.at(i)->FindBin(charge+0.8));

            error_num = sqrt(pow(param_c_value*z2.first.at(i)->GetBinError(0, z2.first.at(i)->FindBin(charge+0.8)), 2) + 
                            pow(param_c_error*z2.first.at(i)->Integral(0, z2.first.at(i)->FindBin(charge+0.8)), 2));

            error_den = sqrt(pow(param_b_value*z1.second.at(i)->GetBinError(0, z1.second.at(i)->FindBin(charge+0.8)), 2) + 
                            pow(param_b_error*z1.second.at(i)->Integral(0, z1.second.at(i)->FindBin(charge+0.8)), 2) + 
                            pow(param_c_value*z2.first.at(i)->GetBinError(0, z2.first.at(i)->FindBin(charge+0.8)), 2) + 
                            pow(param_c_error*z2.first.at(i)->Integral(0, z2.first.at(i)->FindBin(charge+0.8)), 2));
        }*/

        double h = num / den;
        double h_error = h * sqrt(pow(error_num / num, 2) + pow(error_den / den, 2));  // Propagating the error on the ratio

        frac.push_back(h);
        frac_err.push_back(h_error);
        par_b_err.push_back(param_b_error);
        par_c_err.push_back(param_c_error);
        par_b_val.push_back(param_b_value);
        par_c_val.push_back(param_c_value);

        if (z1.first.at(i)->GetEntries() > 0) {
            if (frame->chiSquare() <= 200000) {
                f_vs_r->SetBinContent(Rbins_HighZ[i+1], h);
                f_vs_r->SetBinError(Rbins_HighZ[i+1],h_error);
                if (h==0) {
                    continue;
                }
            }
        } else {
            continue;
        }
    }  
    canvasBin->SetLogx();
    canvasBin->SetLogy(0);
    TAxis *axisY = f_vs_r->GetYaxis();
    TAxis *axisX = f_vs_r->GetXaxis();
    f_vs_r->SetMarkerStyle(20);
    f_vs_r->Draw("L");
    axisX->SetRangeUser(1.8,3100);
    axisY->SetRangeUser(0,0.003);
    canvasBin->Print(out2.Data() );
    canvasBin->Print(out3.Data() );
        for (const auto x : par_b_val) 
        std::cout << "Par b value is : " << x << std::endl;
        for (const auto x : par_c_val) 
        std::cout << "Par c value is : " << x << std::endl;
    auto spline_purity = autospline(f_vs_r,1,8,3,4);
    //get the spline errors
    for (int ii=0; ii<spline_purity->GetNbinsX(); ii++) 
        sp_err.push_back(spline_purity->GetBinError(ii));

    TString fileOut = "../Fragmentation/BelowL1/Fractions/contamination_"+nucleus+".root";

    TFile *output = new TFile(fileOut.Data(), "RECREATE");

    output->WriteTObject(f_vs_r,"f_vs_r");
    output->WriteTObject(spline_purity,"spline_purity");
    /*for (const auto& status : fitStatusList) {
        std::cout << "Fit Status: " << status.status 
                  << ", Covariance Quality: " << status.covQual 
                  << ", Chi Square/ndf : "<< status.chiSquare << std::endl;
    }*/
    output->Close();
}
std::pair<double,double> getBkgRange(unsigned int charge) {
    double xmin,xmax;
    switch(charge) {
        case 15:
            xmin=13.3;
            xmax=15.9;
            break;
        case 16:
            xmin=16.27;
            xmax=16.28;
            break;
        case 17:
            xmin=17.275;
            xmax=17.276;
            break;
    }
    return std::make_pair(xmin,xmax);
} 
std::pair<double,double> getSigRange(unsigned int charge) {
    double xmin,xmax;
    switch(charge) {
        case 14 :
            xmin=12.4;
            xmax=14.9;
            break;
        case 15:
            xmin=13.3;
            xmax=15.9;
            break;
        case 16:
            xmin=14.2;
            xmax=17.0;
            break;
    }
    return std::make_pair(xmin,xmax);
} 
TH1D* shiftHistogram(const TH1D* inputHist, double shiftAmount) {
    if (!inputHist) return nullptr;

    TString name = inputHist->GetName();
    name += "_shifted";

    auto shiftedHist = (TH1D*)inputHist->Clone(name);
    shiftedHist->Reset();

    for (int i = 1; i <= inputHist->GetNbinsX(); ++i) {
        double x = inputHist->GetBinCenter(i);
        double content = inputHist->GetBinContent(i);
        double error = inputHist->GetBinError(i);
        int newBin = shiftedHist->FindBin(x + shiftAmount);

        if (newBin >= 1 && newBin <= shiftedHist->GetNbinsX()) {
            shiftedHist->AddBinContent(newBin, content);
            // Add errors in quadrature if multiple entries land in same bin
            double prevError = shiftedHist->GetBinError(newBin);
            shiftedHist->SetBinError(newBin, std::sqrt(prevError * prevError + error * error));
        }  
    }
    return shiftedHist;
}
TH1D* shiftHistogramResample(const TH1D* inputHist, double shiftAmount) {
    if (!inputHist || inputHist->GetEntries() == 0) return nullptr;

    TString name = inputHist->GetName();
    name += "_resampled_shifted";

    // Clone and reset to get same binning
    TH1D* shiftedHist = (TH1D*)inputHist->Clone(name);
    shiftedHist->Reset();

    // Determine number of samples: use total histogram entries if not provided
    auto nSamples = inputHist->GetEntries();

    for (int i = 0; i < 10e+1*nSamples; ++i) {
        double x = inputHist->GetRandom(); // Sample from the original histogram
        shiftedHist->Fill(x + shiftAmount);
    }
    for (int b = 1; b <= shiftedHist->GetNbinsX(); ++b) {
        shiftedHist->SetBinError(b, sqrt(shiftedHist->GetBinContent(b)) );
    }
    shiftedHist->Scale(inputHist->Integral() / shiftedHist->Integral());
    return shiftedHist;
}
TH1D* shiftHistogramInterpolation(const TH1D* inputHist, double shiftAmount) {
    if (!inputHist || inputHist->GetEntries() == 0) return nullptr;

    TString name = inputHist->GetName();
    name += "_shifted";

    TH1D* shiftedHist = (TH1D*)inputHist->Clone(name);
    shiftedHist->Reset();

    for (int b = 1; b <= inputHist->GetNbinsX(); ++b) {
        double x = inputHist->GetBinCenter(b);
        double xShifted = x - shiftAmount;

        // --- Interpolate bin content (ROOT does this internally with TH1::Interpolate)
        double y = inputHist->Interpolate(xShifted);

        // --- Approximate error: take surrounding bins and interpolate
        int binLeft  = const_cast<TH1D*>(inputHist)->FindBin(xShifted);
        int binRight = binLeft + 1;

        // Clamp to valid bin range
        binLeft  = std::max(1, std::min(binLeft,  inputHist->GetNbinsX()));
        binRight = std::max(1, std::min(binRight, inputHist->GetNbinsX()));

        double xL = inputHist->GetBinCenter(binLeft);
        double xR = inputHist->GetBinCenter(binRight);
        double errL = inputHist->GetBinError(binLeft);
        double errR = inputHist->GetBinError(binRight);

        // Linear interpolation of error (safe, fast)
        double weight = (xShifted - xL) / (xR - xL + 1e-10); // avoid divide by zero
        double err = (1 - weight) * errL + weight * errR;

        // Assign content and error
        shiftedHist->SetBinContent(b, y);
        shiftedHist->SetBinError(b, err);
    }

    return shiftedHist;
}
TH1D* smoothHistogram(const TH1D* inputHist, int nSamples) {
    if (!inputHist || inputHist->GetEntries() == 0) return nullptr;

    TString name = inputHist->GetName();
    name += "_smoothed";

    TH1D* smoothed = (TH1D*)inputHist->Clone(name);
    smoothed->Reset();

    for (int i = 0; i < nSamples; ++i) {
        double x = inputHist->GetRandom();
        smoothed->Fill(x);
    }

    // Assign Poisson errors
    for (int b = 1; b <= smoothed->GetNbinsX(); ++b)
        smoothed->SetBinError(b, std::sqrt(smoothed->GetBinContent(b)));

    // Normalize to original
    double norm = inputHist->Integral();
    if (norm > 0) smoothed->Scale(norm / smoothed->Integral());

    return smoothed;
}
int findRigidityBin(double rigidity, const std::vector<double>& bins) {
    for (size_t i = 0; i < bins.size() - 1; ++i) {
        if (rigidity >= bins[i] && rigidity < bins[i+1]) return i;
    }
    return bins.size() - 2; // assign last bin if above highest edge
}
double initialFractionValue(unsigned int bkgCharge, unsigned int charge, double rig, std::vector<double> Rbins_purity) {
    double f = 0.1;
    int bin = findRigidityBin(rig, Rbins_purity);

    switch (charge) {
        case 14: {
            auto it = charge14_initialFractions.find(bkgCharge);
            if (it != charge14_initialFractions.end() && bin < it->second.size())
                f = it->second[bin];
            break;
        }
        case 15: {
            auto it = charge15_initialFractions.find(bkgCharge);
            if (it != charge15_initialFractions.end() && bin < it->second.size())
                f = it->second[bin];
            break;
        }
        case 16: {
            auto it = charge16_initialFractions.find(bkgCharge);
            if (it != charge16_initialFractions.end() && bin < it->second.size())
                f = it->second[bin];
            break;
        }
        case 18: {
            auto it = charge18_initialFractions.find(bkgCharge);
            if (it != charge18_initialFractions.end() && bin < it->second.size())
                f = it->second[bin];
            break;
        }
        case 20: {
            auto it = charge20_initialFractions.find(bkgCharge);
            if (it != charge20_initialFractions.end() && bin < it->second.size())
                f = it->second[bin];
            break;
        }
    }
    return f;
}
TH1D* removeZeroBins(const TH1D* originalHist) {
    if (!originalHist) return nullptr;

    // Create a new histogram with the same binning
    TH1D* cleanedHist = (TH1D*)originalHist->Clone();
    cleanedHist->SetName((TString(originalHist->GetName()) + "_cleaned").Data());
    cleanedHist->Reset();  // Clear all content

    for (int bin = 1; bin <= originalHist->GetNbinsX(); ++bin) {
        double content = originalHist->GetBinContent(bin);
        double error = originalHist->GetBinError(bin);
        if (content > 0) {
            cleanedHist->SetBinContent(bin, content);
            cleanedHist->SetBinError(bin, error);
        }
        if (content <=0)
            cleanedHist->SetBinContent(bin, 1e-40);
    }

    return cleanedHist;
}
std::tuple<double, int, double> computeChi2FromFrame(RooPlot* frame, const char* modelName, const char* dataHistName, int nFitParams) {
    if (!frame) {
        std::cerr << "[Chi2Calc] Frame is null!" << std::endl;
        return {-1.0, -1, -1.0};
    }

    // Get RooHist from frame
    RooHist* hist = dynamic_cast<RooHist*>(frame->findObject(dataHistName));
    if (!hist) {
        std::cerr << "[Chi2Calc] Data hist \"" << dataHistName << "\" not found in RooPlot!" << std::endl;
        return {-1.0, -1, -1.0};
    }

    // Count only valid data points (e.g., with non-zero error)
    int nPoints = 0;
    for (int i = 0; i < hist->GetN(); ++i) {
        double x, y;
        hist->GetPoint(i, x, y);
        if (hist->GetErrorY(i) > 0) ++nPoints;
    }
    int ndf = nPoints-nFitParams; 

    if (ndf <= 0) {
        std::cerr << "[Chi2Calc] Non-positive ndf (" << ndf << "). Check nFitParams or data!" << std::endl;
        return {-1.0, ndf, -1.0};
    }

    // Compute chi²/ndf from RooPlot
    double chi2_over_ndf = frame->chiSquare(modelName, dataHistName, nFitParams);
    double chi2 = chi2_over_ndf * ndf;

    return {chi2, ndf, chi2_over_ndf};
}
std::tuple<double, double, double, double, int> computeChi2WithBootstrap(
    RooRealVar x,
    TH1D* data,
    TH1D* histSignal,
    std::vector<TH1D*> histBkg,
    unsigned int charge,
    int binIndex,
    const std::vector<double>& Rbins_purity,
    int nFitParams,
    int nBootstrap) {

    std::vector<double> chi2Vals;
    std::vector<int> ndfVals;
    std::vector<double> chi2OverNdfVals;

    TRandom3 rand(0);
    int smt = smoothAmount(charge);

    //Bootstrap iterations
    for (int n = 0; n < nBootstrap; ++n) {  
        //Define the new signal and backgrounds
        TH1D* sigBoot = (TH1D*)histSignal->Clone(Form("sigBoot_%d", n));
        std::vector<TH1D*> bkgBoots;
        for (int ii=0; ii<histBkg.size(); ii++) {
            auto tempBkg = (TH1D*)histBkg[ii]->Clone( Form("bkgBoot_%u_%d",ii, n) );
            bkgBoots.push_back(tempBkg);
        }

        //Resampling signal
        for (int b = 1; b <= histSignal->GetNbinsX(); ++b) {
            if (histSignal->GetBinContent(b)==0) continue;
            double content = rand.Poisson(histSignal->GetBinContent(b));
            if (content>0) {
                sigBoot->SetBinContent(b, content);
                sigBoot->SetBinError(b, sqrt(sigBoot->GetBinContent(b)) );
            } else {
                sigBoot->SetBinContent(b, 10e-40);
                sigBoot->SetBinError(b, 0);
            }
        }

        //Remove 0 bins
        sigBoot->Smooth(smt);
        sigBoot = removeZeroBins(sigBoot);

        //Resampling bkg
        for (int j=0; j<bkgBoots.size(); j++) {
            for (int b = 1; b <= histBkg[j]->GetNbinsX(); ++b) {
                if (histBkg[j]->GetBinContent(b)==0) continue;
                double content = rand.Poisson(histBkg[j]->GetBinContent(b));
                if (content>0) {
                    bkgBoots[j]->SetBinContent(b, content);
                    bkgBoots[j]->SetBinError(b, sqrt(bkgBoots[j]->GetBinContent(b)) );
                } else {
                    bkgBoots[j]->SetBinContent(b, 10e-40);
                    bkgBoots[j]->SetBinError(b, 0);
                }
            }
        }

        //Remove 0 bins
        for (int j=0; j<bkgBoots.size(); j++) {
            bkgBoots[j]->Smooth(smt);
            bkgBoots[j] = removeZeroBins(bkgBoots[j]);
        }


        if (sigBoot->Integral() <= 0) {
            delete sigBoot;
            continue;
        }

        RooFitResult* fitResult = nullptr;
        auto [model,frame] = buildAndFitModel(x, data, sigBoot, bkgBoots, charge, binIndex, Rbins_purity, fitResult, binIndex, n);
        if (!model || !fitResult) {
            delete sigBoot;
            delete fitResult;
            delete model;
            continue;   
        }
        /*if (binIndex==0) {
            TCanvas* c = new TCanvas(Form("c_fit_%d", n), "Fit Result", 800, 600);
            frame->Draw();  // This shows data + model + components already
            c->SetLogy();
            c->SaveAs(Form("fit_result_bootstrap_%d.pdf", n));  // Optional: save to file
        }*/
        // Evaluate chi2
        
        auto [chi2, ndf, chi2ndf] = computeChi2FromFrame(frame, "model", "dataHist", nFitParams);
        
        if (ndf > 0) {  
            chi2Vals.push_back(chi2);
            ndfVals.push_back(ndf);
            chi2OverNdfVals.push_back(chi2ndf);
        }
    }

    int nValid = chi2Vals.size();
    if (nValid == 0)
        return {0., 0., 0., 0., 0};

    double meanChi2 = std::accumulate(chi2Vals.begin(), chi2Vals.end(), 0.0) / nValid;
    double meanChi2OverNdf = meanChi2 / ndfVals[0];

    double stdDevChi2OverNdf = 0.0;
    for (double v : chi2OverNdfVals)
        stdDevChi2OverNdf += (v - meanChi2OverNdf) * (v - meanChi2OverNdf);
    stdDevChi2OverNdf = std::sqrt(stdDevChi2OverNdf / nValid);
    return {meanChi2, ndfVals[0], meanChi2OverNdf, stdDevChi2OverNdf, nValid};
}
std::tuple<RooAddPdf*,RooPlot*> buildAndFitModel(
    RooRealVar x,
    TH1D* data,
    TH1D* histSignal,
    std::vector<TH1D*> histBkg,
    unsigned int charge,
    int i,
    const std::vector<double>& Rbins_purity,
    RooFitResult*& fitResult,
    int bin,
    int nBoot) {

    int smt = smoothAmount(charge);
        
    if (!histSignal || histSignal->GetEntries() == 0) return {nullptr,nullptr};

    histSignal->Smooth(smt);


    RooDataHist dataRoo("dataRoo", "Data", RooArgList(x), Import(*data));
    RooDataHist signalRoo("signalRoo", "Signal", RooArgList(x), Import(*histSignal));
    RooHistPdf* signalPdf = new RooHistPdf("signalPdf", "Signal PDF", RooArgSet(x), signalRoo);

    RooArgList bgPdfList, bgCoeffList;
    RooArgList modelList, modelCoeffs;



    // Loop over each background component and add it with its corresponding coefficient.
    for (int j = 0; j < histBkg.size(); j++) {
        // Build the RooDataHist and then the RooHistPdf for this background component.
        histBkg[j]->Smooth(smt);
        RooDataHist* bgData = new RooDataHist(Form("bgData_%u_%zu", j, i), "Background Data", RooArgList(x), Import(*histBkg[j]));
        RooHistPdf* bgPdf  = new RooHistPdf(Form("bgPdf_%u_%zu", j, i), "Background PDF", RooArgSet(x), *bgData);
        // Add this background PDF to the list.
        bgPdfList.add(*bgPdf);
                
        // Initialize the coefficient for this background component from the histogram integral.
        double initVal = data->Integral()/histBkg[j]->Integral();
        RooRealVar* coeff = new RooRealVar(Form("coeffBg_%u_%zu", j, i),
                                            "Background Coefficient", initialFractionValue(j,charge,Rbins_purity[i],Rbins_purity), 
                                            0., 10e+7*initialFractionValue(j,charge,Rbins_purity[i],Rbins_purity));
        bgCoeffList.add(*coeff);
    }

    //Adding signal
    double initSig = data->Integral() / histSignal->Integral();
    RooRealVar* coeffSig = new RooRealVar("coeffSig", "Signal Coeff",
                                          initialFractionValue(charge, charge, Rbins_purity[i], Rbins_purity),
                                          0., 10e+7 * initialFractionValue(charge, charge, Rbins_purity[i], Rbins_purity) );

    modelList.add(*signalPdf);
    modelCoeffs.add(*coeffSig);

    //Adding bkg
    for (int j = 0; j < histBkg.size(); j++) {
        // The j-th background PDF and its coefficient are already stored in the lists.
         modelList.add(*(dynamic_cast<RooAbsPdf*>(bgPdfList.at(j))));
         modelCoeffs.add(*(dynamic_cast<RooRealVar*>(bgCoeffList.at(j))));
    }

    RooAddPdf* model = new RooAddPdf("model", "Model", modelList, modelCoeffs);
    auto cutLim = cutLimits(charge);
    double chMin = cutLim.first-0.2; double chMax = cutLim.second+0.2;
    x.setRange("fitRange", chMin, chMax);
    fitResult = model->fitTo(dataRoo, Save(true), Range("fitRange"), PrintLevel(-1), Verbose(false));


    RooPlot* frame = x.frame( Name(Form("frame_%d_%d",bin,nBoot) ) );
            dataRoo.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.4),Name("dataHist"));
            model->plotOn(frame,
                        LineColor(kGreen+1),
                        //LineWidth(1.2),
                        FillColor(kGreen - 9),
                        FillStyle(1001),  // Solid fill
                        //DrawOption("F"),  // Fill area
                        Precision(1e-2),
                        Range("fitRange"), NormRange("fitRange"),
                        Name("model_filled"));
            model->plotOn(frame,
                        Components(*signalPdf),
                        LineColor(kAzure),
                        //LineWidth(1.2),
                        FillColor(kAzure - 9),
                        FillStyle(1001),
                        //DrawOption("F"),
                        Precision(1e-2),
                        Range("fitRange"), NormRange("fitRange"),
                        Name("signal_filled"));
            // Loop over background components to plot them individually.
            for (int j = 0; j < bgPdfList.getSize(); j++) {
                RooAbsPdf* bg = dynamic_cast<RooAbsPdf*>(bgPdfList.at(j));
                model->plotOn(frame,
                            Components(*bg),
                            LineColor(kRed + 1.2*j),
                            //LineWidth(1.2),
                            FillColor(kRed + j - 9),  // Slightly lighter fill
                            FillStyle(1001),
                            //DrawOption("F"),
                            Precision(1e-2),
                            Range("fitRange"), NormRange("fitRange"),
                            Name(Form("bg%d_filled", j)));
            }
            model->plotOn(frame,LineColor(kGreen+1),
            LineWidth(1.2),
            Range("fitRange"), NormRange("fitRange"),Name("model"));
    
    return {model,frame};
}
double stepByCharge(unsigned int charge) {
 double step = -1;
 switch(charge) {
    case 14:
        step = 0.015;
        break;
    case 15:
        step = 0.015;
        break;
    case 16:
        step = 0.015;
        break;
    case 18:
        step = 0.02;
        break;
    case 20:
        step = 0.01;
        break;
 }
 return step;
}
std::pair<double,double> rangeShiftByCharge(unsigned int charge) {
  double min,max;
  switch(charge) {
    case 14:
        min = -0.2;
        max = 0.1;
        break;
    case 15:
        min = -0.1;
        max = 0.2;
        break;
    case 16:
        min = -0.2;
        max = 0.1;
        break;
    case 18:
        min = -0.3;
        max = 0.1;
        break;
    case 20:
        min = -0.1;
        max = 0.1;
        break;
    break;
  }
  return make_pair(min,max);
}
void formatFrameSinglePage(std::unique_ptr<RooPlot> &frame, frameValues &FrameValues,
                           std::vector<double> Rbins_purity, TCanvas &c, unsigned int charge) {
    auto i = FrameValues.bin;

    float topHeight = 0.70;
    float bottomHeight = 0.30;

    // First draw the bottom pad (avoid clipping due to z-order)
    TPad* pad2 = new TPad("pad2", "Bottom Pad", 0, 0.0, 1, bottomHeight);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);  // Give enough space for labels
    pad2->Draw();

    // Now draw the top pad
    TPad* pad1 = new TPad("pad1", "Top Pad", 0, bottomHeight, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    // --- Top Pad ---
    pad1->cd();
    gPad->SetLogy(true);
    frame->SetName(Form("fit_bin%zu_shift%.2f", i, FrameValues.shift));
    frame->SetLabelFont(62,"XYZ");
    frame->SetTitleFont(62,"XYZ");
    frame->GetYaxis()->SetTitle("Counts");
    frame->GetXaxis()->SetTitle("Q_{L1} (c.u.)");
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetYaxis()->SetLabelSize(0.065);
    frame->GetYaxis()->SetLabelOffset(0.001);
    frame->GetXaxis()->SetLabelSize(0.046);
    frame->GetYaxis()->SetTitleSize(0.065);
    frame->GetYaxis()->SetTitleOffset(0.6);
    frame->GetXaxis()->SetTitleSize(0.051);
    frame->GetXaxis()->SetTitleOffset(0.86);
    frame->GetYaxis()->SetRangeUser(2., 5*frame->GetMaximum()/100.);
    
    frame->GetXaxis()->SetLabelSize(0);   // Hides X axis tick labels
    frame->GetXaxis()->SetTitleSize(0);   // Hides X axis title

    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTitleFont(62,"XYZ");
    frame->SetTitle("");
    TLegend* leg = new TLegend(0.1, 0.71, 0.4, 0.86);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);   
    leg->SetTextSize(0.045);
    leg->SetTextFont(62);
    leg->AddEntry((TObject*)nullptr, Form("#chi^{2}/ndf = %.2f/%d = %.2f", 
                  FrameValues.chi2,FrameValues.ndf,FrameValues.chi2ndf), "");
    leg->AddEntry((TObject*)nullptr, Form("Shift = %.3f %c %.3f c.u.",FrameValues.shift, 0xB1 , FrameValues.shift_err), "");
    leg->AddEntry((TObject*)nullptr, Form("[%.1f,%.1f] GV",Rbins_purity[i],Rbins_purity[i+1]),"");
    frame->Draw();
    leg->Draw();

    //L1 charge cut limits
    auto cut = cutLimits(charge);
    double x1 = cut.first;  // First x-position
    double x2 = cut.second;  // Second x-position

    TLine* vline1 = new TLine(x1, frame->GetYaxis()->GetXmin(), x1, frame->GetMaximum());
    vline1->SetLineColor(kGray+2);
    vline1->SetLineStyle(1);
    vline1->SetLineWidth(2);
    vline1->Draw("same");

    TLine* vline2 = new TLine(x2, frame->GetYaxis()->GetXmin(), x2, frame->GetMaximum());
    vline2->SetLineColor(kGray+2);
    vline2->SetLineStyle(1);
    vline2->SetLineWidth(2);
    vline2->Draw("same");

    // --- Legend with components ---
    TLegend* leg1 = new TLegend(0.60, 0.71, 0.9, 0.86);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);   
    leg1->SetTextSize(0.045);
    leg1->SetTextFont(62);
    leg1->SetNColumns(2);
    // --- Add styled entries ---
    TObject* dataHist  = frame->findObject("dataHist");
    TObject* model     = frame->findObject("model");
    TObject* signalPdf = frame->findObject("signal_filled");
    TObject* bg1       = frame->findObject( "bg0_filled" );
    TObject* bg2       = frame->findObject( "bg1_filled" );
    TObject* bg3       = frame->findObject( "bg2_filled" );

    if (dataHist) leg1->AddEntry(dataHist,  "Data", "lep");
    if (model)    leg1->AddEntry(model,     "Fit", "l");
    if (bg1)      leg1->AddEntry(bg1,       Form("Z=%u",charge-1), "l"); //charge-1
    if (signalPdf)leg1->AddEntry(signalPdf, Form("Z=%u",charge), "l");   //charge
    if (bg2)      leg1->AddEntry(bg2,       Form("Z=%u",charge+1), "l"); //charge+1
    if (bg3)      leg1->AddEntry(bg3,       Form("Z=%u",charge+2), "l"); //charge+2
    leg1->Draw();


    // --- Bottom Pad ---
    pad2->cd();
    pad2->SetGridy();  // Enable horizontal grid lines
	gStyle->SetGridColor(kGray+3);
    auto pull = frame->pullHist("dataHist", "model", false);
    TH1D* hPull = RooHistToTH1D(pull, "pullHist");
    if (!pull) {
        std::cerr << "[ERROR] pullHist returned nullptr.\n";
        return;
    }
    double mean = hPull->GetMean();
    double sigma = hPull->GetStdDev();

    // Create a small legend in the bottom right corner
    TLegend* legPull = new TLegend(0.1, 0.85, 0.4, 0.99);  // (x1, y1, x2, y2) in NDC
    legPull->SetBorderSize(0);
    legPull->SetFillStyle(0);  // Transparent background
    legPull->SetTextSize(0.1);  // Smaller text for small pad
    legPull->AddEntry((TObject*)nullptr, Form("#mu = %.3f, #sigma = %.3f",mean,sigma),"" );
    legPull->SetTextFont(62);

    pull->SetTitle("");
    pull->GetYaxis()->SetTitle("Pull");
    pull->GetXaxis()->CenterTitle(true);
    pull->GetYaxis()->CenterTitle(true);
    pull->GetXaxis()->SetTitle("Q_{L1} (c.u.)");
    pull->GetYaxis()->SetLabelSize(0.15);
    pull->GetYaxis()->SetTitleSize(0.15);
    pull->GetYaxis()->SetTitleOffset(0.25);
    pull->GetXaxis()->SetLabelSize(0.15);
    pull->GetXaxis()->SetTitleSize(0.15);
    pull->GetXaxis()->SetTitleOffset(1.04);

    pull->GetYaxis()->SetNdivisions(505);
    pull->GetXaxis()->SetNdivisions(505);

    pull->SetMarkerStyle(20);
    pull->SetMarkerSize(0.8);
    pull->Draw("AP");

    // Get axis limits from pull
    double xMin = pull->GetXaxis()->GetXmin();
    double xMax = pull->GetXaxis()->GetXmax();

    TGraph* band = new TGraph(4);
    band->SetPoint(0, xMin, -1);
    band->SetPoint(1, xMax, -1);
    band->SetPoint(2, xMax, 1);
    band->SetPoint(3, xMin, 1);
    band->SetFillColorAlpha(kGray+1, 0.3);
    band->SetLineColor(0);
    band->SetFillStyle(1001);
    band->SetName("sigmaBand");

    legPull->Draw();


    // Draw it *before* pull, so it's behind
    band->Draw("F");       
    pull->Draw("P same");  // Redraw pulls on top

    pull->Draw("P same");
    /*TF1* fitPull = new TF1("fitPull", "gaus", -5, 5);
    hPull->Draw();
    hPull->Fit(fitPull, "Q0");

    double mean = fitPull->GetParameter(1);
    double sigma = fitPull->GetParameter(2);

    TLatex lx;
    lx.SetNDC();
    lx.SetTextSize(0.10);
    lx.DrawLatex(0.15, 0.85, Form("#mu = %.2f", mean));
    lx.DrawLatex(0.15, 0.72, Form("#sigma = %.2f", sigma));*/

}
TH1D* RooHistToTH1D(const RooHist* rh, const TString& name) {
    if (!rh) return nullptr;

    int nPoints = rh->GetN();
    if (nPoints < 1) return nullptr;

    std::vector<double> xVals(nPoints), yVals(nPoints);
    double xLow = 0, xHigh = 0;

    for (int i = 0; i < nPoints; ++i) {
        double x, y;
        rh->GetPoint(i, x, y);
        xVals[i] = x;
        yVals[i] = y;

        double xerrLow = rh->GetErrorXlow(i);
        double xerrHigh = rh->GetErrorXhigh(i);

        if (i == 0)
            xLow = x - xerrLow;
        if (i == nPoints - 1)
            xHigh = x + xerrHigh;
    }
    const auto [ymin, ymax] = std::minmax_element(begin(yVals), end(yVals));
    TH1D* h = new TH1D(name, "", 15, *ymin, *ymax  );
    for (int i = 0; i < nPoints; ++i)
        h->Fill(yVals[i]);

    return h;
}
std::pair<double,double> cutLimits(unsigned int charge) {
    double x1,x2;
    switch(charge) {
        case 9:
            x1 = 7.9;
            x2 = 9.6;
            break;
        case 14:
            x1 = 12.4;
            x2 = 14.9;
            break;
        case 15:
            x1 = 13.3;
            x2 = 15.9;
            break;
        case 16:
            x1 = 14.2;
            x2 = 17.0;
            break;
        case 18:
            x1 = 16.0;
            x2 = 19.1;
            break;
        case 20:
            x1 = 17.8;
            x2 = 21.2;
            break;
    }
    return std::make_pair(x1,x2);
}
static std::pair<double,double> computeL1ChargeEfficiencyFromPdf(
        RooHistPdf* signalPdf,
        RooRealVar& x,
        RooDataHist dataRoo,   // histData convertito in RooDataHist
        double cutLow,
        double cutHigh) {
    if (!signalPdf) return {0.0, 0.0};
    // definisci l'intervallo sul RooRealVar
    x.setRange("cutRange", cutLow, cutHigh);

    // integrale del PDF normalizzato in [cutLow,cutHigh]
    std::unique_ptr<RooAbsReal> numInt(
        signalPdf->createIntegral(x, NormSet(x), Range("cutRange")));
    double numVal = numInt ? numInt->getVal() : 0.0;

    // integrale totale (PDF normalizzato → dovrebbe essere 1)
    std::unique_ptr<RooAbsReal> denomInt(
        signalPdf->createIntegral(x, NormSet(x)));
    double denomVal = denomInt ? denomInt->getVal() : 1.0;

    if (denomVal <= 0.0) return {0.0, 0.0};

    double eff = numVal / denomVal;

    // --- Stima errore statistico dai dati -------------------
    // Creo un TH1 dai dati per propagare in modo "binomial-like"
    std::unique_ptr<TH1> htmp(dataRoo.createHistogram("htmp", x));
    if (!htmp) return {eff, 0.0};

    int binLow  = htmp->GetXaxis()->FindBin(cutLow);
    int binHigh = htmp->GetXaxis()->FindBin(cutHigh);

    double numH   = htmp->Integral(binLow, binHigh);
    double denomH = htmp->Integral(1, htmp->GetNbinsX());

    if (denomH <= 0.0) return {eff, 0.0};

    // Propagazione errori come nel binomial
    double err2_num = 0.0;
    double err2_tot = 0.0;
    for (int b = 1; b <= htmp->GetNbinsX(); ++b) {
        double e = htmp->GetBinError(b);
        if (b >= binLow && b <= binHigh) err2_num += e*e;
        err2_tot += e*e;
    }

    double term1 = err2_num / (denomH * denomH);
    double term2 = (numH * numH) * err2_tot / std::pow(denomH,4);
    double var = term1 + term2;
    if (var < 0) var = 0;

    double effErr = std::sqrt(var);

    return {eff, effErr};
}
int smoothAmount(unsigned int charge) {
    int a = 1;
    switch(charge) {
        case 18:
        case 20:
            a = 4;
            break;
        default:
            a = 2;
    }
    return a;
}