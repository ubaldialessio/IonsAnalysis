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
//additional header
#include <string_view>
#include "TLatex.h"
#include "TPaveText.h"
static constexpr auto EXENAME = R"(ShiftTemplates)";
using namespace RooFit;

std::map<std::string, std::vector<TH1D*>> getTemplates(unsigned int charge); 
std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> loadTemplates(unsigned int charge);
std::pair<double, RooPlot*> chiFromFit(unsigned int charge, TH1D *data, TH1D *model, double shift); //this thing returns chi2 and fit
void shiftTemplates(unsigned int charge, std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> templates);
TH1D* removeZeroBins(const TH1D* originalHist);
TH1D* shiftHistogram(const TH1D* inputHist, double shiftAmount);

std::map<unsigned int, std::vector<double>> charge15_initialFractions = {
    {14, {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.035, 0.045, 0.03, 0.03, 0.03}},
    {15, {0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7}},
    {16, {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3}},
    {17, {0.11, 0.09, 0.09, 0.09, 0.04, 0.04, 0.04, 0.04, 0.06, 0.04, 0.04, 0.04}}
};

int main(int argc, char *argv[]) {
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    TH1::SetDefaultSumw2(true);
	if (argc < 2) {
		printf("Usage: \n");
	    printf("%s <charge> \n", argv[0]);
	    return 1;
	}
    unsigned int charge=atoi(argv[1]);
    auto templates = loadTemplates(charge);
    shiftTemplates(charge, templates);
}

std::map<std::string, std::vector<TH1D*>> getTemplates(unsigned int charge) {
    std::map<std::string, std::vector<TH1D*>> templates;
    // Costruisci il nome della lista basato sulla carica
    TString listName1 = "L1TemplateList_" + TString::Format("%d", charge);
    TString listName2 = "L2TemplateList_" + TString::Format("%d", charge);
    // Apri il file ROOT
    TString path1 = "../Fragmentation/BelowL1/";
    TString nucleus = getIonName(charge);
    TString fileName = path1 + nucleus + "/" + nucleus + ".root";
    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
        printf("Errore nell'aprire il file\n");
        return templates;
    }
    // Recupera le liste
    TList *list1 = dynamic_cast<TList*>(file->Get(listName1.Data()));
    TList *list2 = dynamic_cast<TList*>(file->Get(listName2.Data()));
    if (!list1 || !list2) {
        printf("Errore nel recuperare le liste\n");
        return templates;
    }
    // Get L1 Templates
    for (auto *obj : *list1) {
        TH1D *hist1 = dynamic_cast<TH1D*>(obj);
        if (hist1) {
            templates["L1"].push_back(hist1);
        } else {
            printf("Errore nel recuperare un istogramma da L1TemplateList per carica %u\n", charge);
        }
    }
    // Get L2 Templates
    for (auto *obj : *list2) {
        TH1D *hist2 = dynamic_cast<TH1D*>(obj);
        if (hist2) {
            templates["L2"].push_back(hist2);
        } else {
            printf("Errore nel recuperare un istogramma da L2TemplateList per carica %u\n", charge);
        }
    }
    return templates;
}
std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> loadTemplates(unsigned int charge) {
    std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> templatesByCharge;
    //Original with charge, charge+1, charge+2
    /*for (unsigned int i = 0; i <= n; i++) {
        unsigned int currentCharge = charge + i;
        templatesByCharge[currentCharge] = getTemplates(currentCharge); 
    }*/
    //Now including also charge -1 
    std::vector<unsigned int> allowedZ = {charge-1,charge,charge+1,charge+2};
    for (auto k : allowedZ)
        templatesByCharge[k] = getTemplates(k); 
    return templatesByCharge;
}
std::pair<double, RooPlot*> chiFromFit(unsigned int charge, TH1D *data, TH1D *model, double shift) {
    //this thing returns me the chi2 for a given bin (namely a given data and template)
    // --- 1. Set up the output
    double chi2 = 999;
    RooRealVar x("x", "Q_{L1}", charge - 3, charge + 3);
    // Shifting (need to check the function maybe)
    if (shift != 0) {
        TH1D* shiftedTemplate = shiftHistogram(model, shift);
        if (shiftedTemplate) model = shiftedTemplate;
    }
    // Remove zero bins
    model = removeZeroBins(model);
    // --- 2. Set up the RooFit data and template
    RooDataHist dataRoo("dataRoo", "Data", RooArgList(x), Import(*data));
    RooDataHist templateRoo("templateRoo", "Template", RooArgList(x), Import(*model));
    RooHistPdf  templatePdf("templatePdf", "Template PDF", RooArgSet(x), templateRoo);
    // --- 3. Fit the template to the data 
    x.setRange("integrationRange", charge-2, charge+2);
    double init = data->Integral()/model->Integral();
    RooRealVar coeff("coeff", "Coefficient",0.03*init,init,2*init);
    RooRealSumPdf modelToFitWith("model", "Model", templatePdf, coeff, kTRUE);
    RooFitResult* fitResult = modelToFitWith.fitTo(dataRoo, Save(),Timer(true),RecoverFromUndefinedRegions(10), 
                                                   EvalErrorWall(true),  Extended(true),Range("integrationRange"));
    // --- 4. Retrieving the chi2
    x.setRange("defaultRange", charge - 3, charge + 3);
    RooPlot *frame = x.frame();
    dataRoo.plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.5));
    modelToFitWith.plotOn(frame,LineColor(kRed),LineColor(kRed),LineWidth(1.5), Range("defaultRange"), NormRange("integrationRange"));
    //modelToFitWith.plotOn(frame,LineColor(kRed),Components("templatePdf"), Range("defaultRange"), NormRange("integrationRange"));
    chi2 = frame->chiSquare(1);
    return {chi2, frame};
}
void shiftTemplates(unsigned int charge, std::map<unsigned int, std::map<std::string, std::vector<TH1D*>>> templates) {
    //Actually I want two things:
    // 1. The 2D plots chi2 (of the fit of the under study L1 distribution with only signal) vs some shifts
    // 2. The shifted templates

    const int nP_rigBins = 13; 
	const double P_rigBins[nP_rigBins] = {0.8,3.64,5.37,7.76,11.00,15.3,21.1,38.9,70.00,150.00,300.00,500.00,1000.00};	

    // --- 0. Set up the output file
    std::vector<RooPlot*> allFrames;
    auto ion = getIonName(charge).Data();
    TString outPath = Form("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/%s/%s_chi.root",
                          ion,ion);
    auto outfile = new TFile(outPath.Data(), "recreate");

    // --- 1. Set up data and templataes
    auto datas  = templates[charge]["L1"];
    auto models = templates[charge]["L2"];
    const int nBins = datas.size();

    // --- 2. Define the shift scan
    double step = 0.01;
    double shift_min = -0.1, shift_max = 0.1;
    int nShiftSteps = static_cast<int>((shift_max - shift_min) / step) + 1;

    // Loop over each bins
    for (size_t i = 0; i < static_cast<size_t>(nBins); ++i) {
        TH1D* data  = datas[i];
        TH1D* model = models[i];
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
        for (int p = 1; p <= 20; ++p) {
            cAll->cd(p);
            gPad->SetLogy(1);  // <- Set logY for each sub-pad
        }
        // --- 3. Set up the output namely chi vs shift 
        TH2D* chi_vs_shift = new TH2D(Form("chi_vs_shift-bin %zu",i),
                             Form("%s (Z=%d) @ %s;Shift(c.u.);#chi^{2}", getIonName(charge).Data(), charge, data->GetTitle()), 
                             nShiftSteps, shift_min, shift_max, 100, 0, 100);
        //Loop over shift values
        for (int s = 0; s < nShiftSteps; ++s) {
            double currentShift = shift_min + s * step;
            auto [chi, frame] = chiFromFit(charge, data, model, currentShift);
            chi_vs_shift->Fill(currentShift, chi);
            if (frame) {
                cAll->Update();
                frame->SetName(Form("fit_bin%zu_shift%.2f", i, currentShift));
                frame->GetYaxis()->SetTitle("Counts");
                frame->GetYaxis()->SetLabelSize(0.062);
                frame->GetXaxis()->SetLabelSize(0.062);
                frame->GetYaxis()->SetTitleSize(0.062);
                frame->GetXaxis()->SetTitleSize(0.062);
                frame->GetYaxis()->SetTitleOffset(0.8);
                frame->GetXaxis()->SetTitleOffset(0.62);
                frame->GetYaxis()->SetRangeUser(1., frame->GetMaximum()*1.);
                gStyle->SetLabelFont(62,"XYZ");
                gStyle->SetTitleFont(62,"XYZ");
                frame->SetTitle("");
                allFrames.push_back(frame);
                cAll->cd(s+1);
                frame->Draw();
                TLegend* leg = new TLegend(0.65, 0.83, 0.88, 0.99);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.074);
                leg->AddEntry((TObject*)nullptr, Form("#chi^{2} = %.2f", chi), "");
                leg->AddEntry((TObject*)nullptr, Form("Shift = %.2f", currentShift), "");
                leg->Draw();
                if (s==0) {
                    TLegend* rigLeg = new TLegend(0.10, 0.83, 0.26, 0.99);
                    rigLeg->SetBorderSize(0);
                    rigLeg->SetFillStyle(0);
                    rigLeg->SetTextSize(0.074);
                    rigLeg->AddEntry((TObject*)nullptr, Form("[%.1f,%.1f]",P_rigBins[i],P_rigBins[i+1]),"");
                    rigLeg->Draw();
                }

            }
        }
        cAll->Print("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"+getIonName(charge)+
                    Form("/%s_shift.pdf",ion));
        outfile->WriteTObject(chi_vs_shift,Form("chi_vs_shift_bin%zu",i));
        if (i== nBins-1) cAll->Print("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"+getIonName(charge)+
                                    Form("/%s_shift.pdf]",ion));
    }
    outfile->Close();
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