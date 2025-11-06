#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"
#include "SplineUtility.h"
#include <TParameter.h>

#include <TGraphAsymmErrors.h>

#include "TCanvas.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <charge>\n", argv[0]);
        return 1;
    }
    unsigned int charge = atoi(argv[1]);
    auto pub = GetPublishedFlux(charge);
    TCanvas* c = new TCanvas();
    c->SaveAs( Form("../output/flux_fitTest_%s.pdf[",getIonName(charge).Data()) );

    //Fit
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetOptStat(0);

    std::vector<double> xvals, yerr68;
    for (int i = 0; i < nRbins_HighZ - 1; ++i) {
        double x_center = 0.5 * (Rbins_HighZ[i] + Rbins_HighZ[i+1]);
        xvals.push_back(x_center);
    }
    std::vector<double> knots = {1,3,5,8,10,30,50,80,100,300,900};
    auto fit = spfit(pub, 3, 1., 1000, knots, "", &xvals, &yerr68);
    std::cout << fit->GetTitle() << std::endl;

    pub->SetTitle(fit->GetTitle());
    pub->GetXaxis()->SetTitle("R (GV)");
    pub->GetYaxis()->SetTitle("#phi R^{2.7} (m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    pub->SetMarkerStyle(20);
    pub->GetXaxis()->SetRangeUser(2.,3000);
    fit->SetRange(0.1,3000);
    formatAxis(pub,1);
    c->SetLogx();
    pub->Draw();
    fit->Draw("SAME");
    c->SaveAs( Form("../output/flux_fitTest_%s.pdf",getIonName(charge).Data()) );

    //Recover
    auto recover = MultiplyTF1ByXPower(fit,-2.7);
    MultiplyByXPower(pub,-2.7);
    c->SetLogy();
    pub->Draw();
    recover->Draw("SAME");

    c->SaveAs( Form("../output/flux_fitTest_%s.pdf",getIonName(charge).Data()) );
    
    //Output
    c->SaveAs( Form("../output/flux_fitTest_%s.pdf]",getIonName(charge).Data()) );
}