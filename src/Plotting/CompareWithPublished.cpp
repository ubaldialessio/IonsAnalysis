#include "binning.h"
#include "utils.h"

#include <TGraphAsymmErrors.h>

#include "TCanvas.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

struct FluxComparison {
    TString label;          //name             
     TH1D* flux;  //flux
     TH1D* fluxOverAverage;
     TH1D* mineOverAverage;
     TH1D* errorContour;
};

void drawTwoPanel(TH1* mine, TH1* other, TH1* otherOverAverage, TH1* mineOverAverage, 
                  TH1* averageContour,  TString titleMine,  TString titleOther, Int_t color, 
                  TPad *pad1, TPad *pad2, TCanvas *c, unsigned int charge);
TH1D *DiffOverRatio( TH1D *pub,  TH1D *my);
short int GetFluxColor(unsigned int charge);
std::vector<FluxComparison> BuildOtherFluxes(unsigned int charge, TString timePeriod);
std::pair<double,double> fluxRange(unsigned int charge);

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s <time> <charge>\n", argv[0]);
        return 1;
    }
    TString timePeriod = argv[1];
    unsigned int charge = atoi(argv[2]);
    //Mine flux
     auto mineFlux = GetFlux(charge, timePeriod);

    //Others
    auto comps = BuildOtherFluxes(charge, timePeriod);

    // Create canvas and divide into two pads: top for flux, bottom for ratio
    TCanvas* c = new TCanvas();
    float topHeight = 0.70;
    float bottomHeight = 0.30;

    // First draw the bottom pad (avoid clipping due to z-order)
    TPad* pad2 = new TPad("pad2", "Bottom Pad", 0, 0.0, 1, bottomHeight);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);  // Give enough space for labels
    // Now draw the top pad
    TPad* pad1 = new TPad("pad1", "Top Pad", 0, bottomHeight, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad2->Draw();
    pad1->Draw();
    c->SaveAs( Form("../output/Comparison/%s_%s.pdf[",getIonName(charge).Data(), timePeriod.Data()) );

    // Loop over and draw
    for (auto& comp : comps) {
         TH1* otherFlux = comp.flux;
         TH1* contourComp = comp.errorContour;
         int color = GetFluxColor(charge);
         drawTwoPanel(mineFlux, otherFlux, comp.fluxOverAverage, comp.mineOverAverage, contourComp,
                     Form("%s_%s", getIonName(charge).Data(), timePeriod.Data()),
                     comp.label,
                     color,
                     pad1,pad2,c,charge);
    }
    c->SaveAs( Form("../output/Comparison/%s_%s.pdf]",getIonName(charge).Data(), timePeriod.Data()) );

}

std::vector<FluxComparison> BuildOtherFluxes(unsigned int charge, TString timePeriod) {
     auto mine = GetFlux(charge,timePeriod);
     auto jack = getJackFlux(charge);
     auto qi   = RebinHistogramWeighted(GetQYanFlux(charge));
     auto jose = RebinHistogramWeighted(GetJoseFlux(15));
     auto zhen = GetZhenFlux(charge); MultiplyByXPower(zhen,2.7);
     auto pub  = RebinHistogramWeighted(GetPublishedFlux(charge));

     auto jack_av = getAverage(mine,jack);
     auto qi_av   = getAverage(mine,qi);
     auto jose_av = getAverage(mine,jose);
     auto zhen_av = getAverage(mine,zhen);
     auto pub_av  = getAverage(mine,pub);

     auto jack_ov = FluxOverAverage(jack,jack_av,"");
     auto qi_ov   = FluxOverAverage(qi,qi_av,"");
     auto jose_ov = FluxOverAverage(jose,jose_av,"");
     auto zhen_ov = FluxOverAverage(zhen,zhen_av,"");
     auto pub_ov  = FluxOverAverage(pub,pub_av,"");

     auto mine_ov_jack = FluxOverAverage(mine,jack_av,"");
     auto mine_ov_qi  = FluxOverAverage(mine, qi_av,"");
     auto mine_ov_jose = FluxOverAverage(mine,jose_av,"");
     auto mine_ov_zhen = FluxOverAverage(mine,zhen_av,"");
     auto mine_ov_pub  = FluxOverAverage(mine,pub_av,"");

     auto jack_err = ErrorContour(jack,jack_av);
     auto qi_err =   ErrorContour(qi,qi_av);
     auto jose_err = ErrorContour(jose, jose_av);
     auto zhen_err = ErrorContour(zhen,zhen_av);
     auto pub_err  = ErrorContour(pub,pub_av);

    

    // Define all comparisons
    std::vector<FluxComparison> comps = {
        {"MIT (J. Bolster)",       jack, jack_ov, mine_ov_jack, jack_err},
        {"IHEP (Q. Yan)",         qi,   qi_ov,   mine_ov_qi,   qi_err},
        {"Jose",       jose, jose_ov, mine_ov_jose, jose_err},
        {"Rome (Z. Liu)",zhen, zhen_ov, mine_ov_zhen, qi_err},
        {"Pub",        pub,  pub_ov,  mine_ov_pub,  pub_err}
    };
    return comps;
}
/*int main(int argc, char *argv[]) {
	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <time> <charge> \n", argv[0]);
	    return 1;
	}
    setLogon();
    gStyle->SetLabelFont(62,"XYZ");
    TString timePeriod = argv[1];
    unsigned int charge=atoi(argv[2]);

    printComparison(charge,timePeriod);

    auto color = GetFluxColor(charge);
    // Vectors to store the published data
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> values;
    std::vector<double> syst;
    std::vector<double> errors; // To store the combined error
    
    readFluxTable(lower_bounds,upper_bounds,values,errors,charge,syst);
     auto published = buildPublished(lower_bounds,upper_bounds,values,errors);
     auto myflux = GetFlux(charge,timePeriod);
     auto ratio = DiffOverRatio(published,myflux);
     auto jackP = getJackFlux(charge);
     auto myflux4 = GetFlux(charge,timePeriod);

    std::vector<double> average = getAverage(published,myflux);
     auto h1 = FluxOverAverage(myflux,average,"mine");
     auto h2 = FluxOverAverage(published,average,"");
     auto h3 = FluxOverErrors(myflux,published);
     auto limits = ErrorContour(published,average);


     auto qyanOriginal = GetQYanFlux(charge);
     auto qyan = RebinHistogramWeighted(qyanOriginal);
     auto myflux2 = GetFlux(charge,timePeriod);
    std::vector<double> averageQYan = getAverage(qyan,myflux2);
     auto h4 = FluxOverAverage(myflux2,averageQYan,"mine");
     auto h5 = FluxOverAverage(qyan,averageQYan,"");
     auto h6 = FluxOverErrors(myflux2,qyan);
     auto lim = ErrorContour(qyan,averageQYan);

    //jose
     auto jose = GetJoseFlux(15);
     auto joseP = RebinHistogramWeighted(jose);
     auto myflux3 = GetFlux(charge,timePeriod);
    std::vector<double> averageJose = getAverage(joseP,myflux3);
     auto myOverAvJose = FluxOverAverage(myflux3,averageJose,"mine");
     auto JoseOverAvJose = FluxOverAverage(joseP,averageJose,"");
     auto lim2 = ErrorContour(joseP,averageJose);

    //MIT over average with published
     auto pubRebin = RebinHistogramWeighted(published);
    std::vector<double> averageMitPub = getAverage(qyan,pubRebin);
     auto h10 = FluxOverAverage(pubRebin,averageMitPub,"mine");
     auto h11 = FluxOverAverage(qyan,averageMitPub,"");
     auto lim10 = ErrorContour(qyan,averageMitPub);

    //Jack
     auto jackRebin = RebinHistogramWeighted(jackP);
    std::vector<double> averageJackMine = getAverage(jackRebin,myflux4);
     auto mineOverJack = FluxOverAverage(myflux4,averageJackMine,"");
     auto JackOverJack = FluxOverAverage(jackRebin,averageJackMine,"");
     auto errJack = ErrorContour(jackRebin,averageJackMine);

    double upper_range = 1;
    if (charge == 14 )
        upper_range = 35;
    if (charge == 16 )
        upper_range = 6;
    if (charge == 15)
        upper_range = 0.6;

    TCanvas* c1 = new TCanvas();
    TString out = "../output/Comparison/"+getIonName(charge)+Form("_%s.pdf",timePeriod.Data());
    std::cout << out << std::endl;
    setLogon();
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0.00001);
    setLogon();
    gStyle->SetLabelFont(62,"XYZ");
    c1->SaveAs( Form("%s[", out.Data() ) );
    c1->SetLogx();
	c1->SetGridy();
    myflux->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
    myflux->GetXaxis()->SetTitle("R (GV)");
    myflux->GetYaxis()->SetTitle("#phi R^{2.7} (m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    myflux->GetYaxis()->SetRangeUser(0.089,upper_range);
    myflux->GetXaxis()->SetRangeUser(2.9,1000);
    myflux->SetMarkerColor(color);
    myflux->SetLineColor(color);
    myflux->SetMarkerStyle(20); 
    myflux->SetMarkerSize(0.75);
    myflux->GetYaxis()->SetRangeUser(0,upper_range);
    formatAxis(myflux,1);
    myflux->Draw();
    c1->SaveAs(out.Data() );

    c1->Update();
    published->SetMarkerStyle(20);
    published->SetMarkerSize(0.75);
    published->Draw();
    published->SetTitle("Published flux");
    c1->Update();
    c1->SaveAs(out.Data() );

    //my vs pub
    myflux->Draw();
    published->Draw("SAME");
    myflux->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
    myflux->GetXaxis()->SetTitle("R (GV)");
    myflux->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    myflux->GetYaxis()->SetRangeUser(0.2,upper_range);
    myflux->GetXaxis()->SetRangeUser(2.15,1000);
    TLegend *a = new TLegend(0.1,0.9,0.25,0.75);
    a->AddEntry(myflux,"My flux");
    a->AddEntry(published, "Pub");
    a->SetTextSize(0.035); 
    a->Draw("SAME");
    c1->SetLogy();
    c1->Update();
    c1->SaveAs(out.Data() );

    //my vs pub over average
    c1->SetLogy(0);
    limits->SetFillColorAlpha(kBlue,0.15);
    h1->Draw("hist P");
    limits->Draw("P E3 SAME");
    h2->SetLineColor(kBlack);
    h2->Draw("hist P SAME");
    h1->SetMarkerColor(color);
    h1->SetLineColor(color);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(0.75);
    h1->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
    h1->GetXaxis()->SetTitle("R (GV)");
    h1->GetYaxis()->SetTitle("Flux/Average");
    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(0.75);
    h1->GetYaxis()->SetRangeUser(0.7,1.3);
    h1->GetXaxis()->SetRangeUser(2.15,1000);
    formatAxis(h1, 1);

    h1->SetTitleFont(62,"XYZ");
    h2->SetTitleFont(62,"XYZ");
    limits->SetTitleFont(62,"XYZ");
    h1->SetTitleSize(0.05,"XYZ");
    h2->SetTitleSize(0.05,"XYZ");
    limits->SetTitleSize(0.05,"XYZ");
    h1->SetLabelFont(62,"XYZ");
    h2->SetLabelFont(62,"XYZ");
    limits->SetLabelFont(62,"XYZ");
    h1->SetLabelSize(0.05,"XYZ");
    h2->SetLabelSize(0.05,"XYZ");
    limits->SetLabelSize(0.05,"XYZ");

    lim2->SetTitleFont(62,"XYZ");
    lim2->SetTitleSize(0.05,"XYZ");
    lim2->SetLabelFont(62,"XYZ");
    lim2->SetLabelSize(0.05,"XYZ");

    TLegend *b = new TLegend(0.70,0.9,0.9,0.7);
    b->AddEntry(h1,"My flux");
    b->AddEntry(h2, "Pub");
    b->AddEntry(limits, "Pub total errors");
    b->SetTextSize(0.035); 
    b->Draw("SAME");

    c1->Update();
    c1->SaveAs(out.Data() );
    h3->Draw();
    h3->SetMarkerStyle(20);
    h3->GetXaxis()->SetTitle("R (GV)");
    h3->SetTitleFont(62,"XYZ");
    h3->SetTitleSize(0.05,"XYZ");
    h3->SetLabelFont(62,"XYZ");
    h3->SetLabelSize(0.05,"XYZ");
    c1->Update();
    c1->SaveAs(out.Data() );

    //my vs qyan over average
    c1->SetLogy(0);
    lim->SetFillColorAlpha(kBlue,0.15);
    h4->Draw("hist P");
    lim->Draw("P E3 SAME");
    h5->SetLineColor(kBlack);
    h5->Draw("hist P SAME");
    h4->SetMarkerColor(color);
    h4->SetLineColor(color);
    h4->SetMarkerStyle(20);
    h4->SetMarkerSize(0.75);
    h4->SetTitle(Form("%s Flux", getIonName(charge).Data() ) );
    h4->GetXaxis()->SetTitle("R (GV)");
    h4->GetYaxis()->SetTitle("Flux/Average");
    h5->SetMarkerStyle(20);
    h5->SetMarkerSize(0.75);
    h4->GetYaxis()->SetRangeUser(0.7,1.32);
    h4->GetXaxis()->SetRangeUser(2.9,1000);
    formatAxis(h4, 1);

    h4->SetTitleFont(62,"XYZ");
    h5->SetTitleFont(62,"XYZ");
    lim->SetTitleFont(62,"XYZ");
    h4->SetTitleSize(0.05,"XYZ");
    h5->SetTitleSize(0.05,"XYZ");
    lim->SetTitleSize(0.05,"XYZ");
    h4->SetLabelFont(62,"XYZ");
    h5->SetLabelFont(62,"XYZ");
    lim->SetLabelFont(62,"XYZ");
    h4->SetLabelSize(0.05,"XYZ");
    h5->SetLabelSize(0.05,"XYZ");
    lim->SetLabelSize(0.05,"XYZ");

    TLegend *b1 = new TLegend(0.70,0.9,0.9,0.7);
    b1->AddEntry(h4,"My flux");
    b1->AddEntry(h5, "MIT");
    b1->AddEntry(lim, "MIT tot errors");
    b1->SetTextSize(0.035); 
    b1->Draw("SAME");
    c1->Update();
    c1->SaveAs(out.Data() );

    //my vs qyan
    myflux2->SetMarkerColor(color);
    myflux2->SetLineColor(color);
    myflux2->SetMarkerStyle(20); 
    myflux2->SetMarkerSize(0.75);
    myflux2->SetTitle("My flux");
    myflux2->GetYaxis()->SetRangeUser(0.01,upper_range);
    qyan->Draw();
    myflux2->Draw("SAME");
    qyan->SetMarkerStyle(20);
    qyan->SetMarkerSize(0.75);
    qyan->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
    qyan->GetXaxis()->SetTitle("R (GV)");
    qyan->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    //qyan->GetYaxis()->SetRangeUser(0.2,upper_range);
    qyan->GetXaxis()->SetRangeUser(2.9,1000);
    TLegend *a1 = new TLegend(0.1,0.9,0.25,0.75);
    a1->AddEntry(myflux,"My flux");
    a1->AddEntry(qyan, "MIT");
    a1->SetTextSize(0.035); 
    a1->Draw("SAME");
    c1->SetLogy(0);
    c1->Update();
    c1->SaveAs(out.Data() );

    //my vs jose
    joseP->SetMarkerStyle(20);
    joseP->SetMarkerSize(0.75);
    joseP->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
    joseP->GetXaxis()->SetTitle("R (GV)");
    joseP->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    myflux3->SetMarkerColor(color);
    myflux3->SetLineColor(color);
    myflux3->SetMarkerStyle(20); 
    myflux3->SetMarkerSize(0.75);
    myflux3->Draw();
    joseP->Draw("SAME");
    myflux3->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
    myflux3->GetXaxis()->SetTitle("R (GV)");
    myflux3->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    myflux3->GetYaxis()->SetRangeUser(0.2,upper_range);
    myflux3->GetXaxis()->SetRangeUser(2.15,1000);
    myflux3->SetMarkerColor(color);
    myflux3->SetLineColor(color);
    TLegend *aa = new TLegend(0.1,0.9,0.25,0.75);
    aa->AddEntry(myflux3,"My flux");
    aa->AddEntry(joseP, "Jose");
    aa->SetTextSize(0.035); 
    aa->Draw("SAME");
    c1->SetLogy(0);
    c1->Update();
    c1->SaveAs(out.Data() );

    c1->Divide(1, 2, 0, 0); // Divide into two sub-pads
    c1->cd(1);
        gPad->SetLogx(); // Use gPad instead of c1
        gPad->SetGridy();
        qyan->SetTitle("");
        qyan->GetYaxis()->SetTitleOffset(0.55);
        qyan->GetYaxis()->SetLabelSize(0.065);
        myflux2->GetXaxis()->SetRangeUser(2.67,1000);
        qyan->GetXaxis()->SetRangeUser(2.67,1000);
        myflux2->SetTitle("");
        qyan->Draw();
        myflux2->Draw("SAME");
        TLegend *a2 = new TLegend(0.1, 0.9, 0.25, 0.75);
        a2->AddEntry(myflux, "My flux");
        a2->AddEntry(qyan, "MIT");
        a2->SetTextSize(0.035); 
        //a2->Draw("SAME");

    c1->cd(2);
        gPad->SetLogx(); // Explicitly set log x for the second pad
        gPad->SetGridy();
        gPad->SetLogy(0); // Ensure log y is off if necessary
        h4->GetYaxis()->SetTitleOffset(0.55);
        h4->GetYaxis()->SetLabelSize(0.058);
        h4->GetXaxis()->SetLabelSize(0.06);
        h4->GetXaxis()->SetTitleOffset(0.95);
        h4->GetXaxis()->SetRangeUser(2.67,1000);
        h5->GetXaxis()->SetRangeUser(2.67,1000);
        h4->SetTitle("");
        h5->SetTitle("");
        h4->Draw("hist P");
        lim->Draw("P E3 SAME");
        h5->Draw("hist P SAME");

        TLegend *b2 = new TLegend(0.85, 0.85, 1., 1.);
        b2->AddEntry(h4, "My flux");
        b2->AddEntry(h5, "MIT");
        b2->AddEntry(lim, "MIT tot errors");
        b2->SetTextSize(0.035); 
        b2->Draw("SAME");

    c1->Update();
    c1->SaveAs(out.Data());

    //my vs jose over average
    myflux3->GetXaxis()->SetTitle("R (GV)");
    myflux3->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    c1->Update();
    c1->cd(1);
    c1->Update();
        gPad->SetLogx(); // Use gPad instead of c1
        gPad->SetGridy();
        joseP->SetTitle("");
        joseP->GetYaxis()->SetTitleOffset(0.55);
        joseP->GetYaxis()->SetLabelSize(0.065);
        joseP->GetXaxis()->SetRangeUser(2.9,1000);
        myflux3->GetXaxis()->SetRangeUser(2.9,1000);
        myflux3->SetTitle("");
        joseP->Draw();
        myflux3->Draw("SAME");
        TLegend *a3 = new TLegend(0.1, 0.9, 0.25, 0.75);
        a3->AddEntry(myflux3, "My flux");
        a3->AddEntry(joseP, "Jose");
        a3->SetTextSize(0.035); 
        a3->Draw("SAME");

    c1->cd(2);
        gPad->SetLogx(); // Explicitly set log x for the second pad
        gPad->SetGridy();
        gPad->SetLogy(0); // Ensure log y is off if necessary
        myOverAvJose->SetMarkerColor(color);
        myOverAvJose->GetYaxis()->SetTitleOffset(0.55);
        myOverAvJose->GetYaxis()->SetLabelSize(0.058);
        myOverAvJose->GetXaxis()->SetLabelSize(0.06);
        myOverAvJose->GetXaxis()->SetTitleOffset(0.95);
        myOverAvJose->GetYaxis()->SetRangeUser(0.7,1.3);
        myOverAvJose->GetXaxis()->SetRangeUser(2.67,1000);
        myOverAvJose->GetXaxis()->SetTitle("R (GV)");
        myOverAvJose->GetYaxis()->SetTitle("Flux/Average");
        myOverAvJose->SetTitle("");
        JoseOverAvJose->SetTitle("");
        myOverAvJose->Draw("hist P");
        lim2->SetFillColorAlpha(kBlue,0.15);
        lim2->Draw("P E3 SAME");
        JoseOverAvJose->Draw("hist P SAME");

        TLegend *b3 = new TLegend(0.85, 0.85, 1., 1.);
        b3->AddEntry(myOverAvJose, "My flux");
        b3->AddEntry(JoseOverAvJose, "Jose");
        b3->AddEntry(lim2, "Jose tot errors");
        b3->SetTextSize(0.035); 
        b3->Draw("SAME");

    c1->Update();
    c1->SaveAs(out.Data());

    //my vs pub over average
    c1->cd(1);
        gPad->SetLogx(); // Use gPad instead of c1
        gPad->SetGridy();
        myflux->SetMarkerColor(color);
        myflux->SetLineColor(color);
        myflux->SetMarkerStyle(20); 
        myflux->SetMarkerSize(0.75);
        myflux->SetTitle("My flux");
        myflux->GetYaxis()->SetRangeUser(0.01,upper_range);
        published->SetMarkerStyle(20);
        published->SetMarkerSize(0.75);
        published->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
        published->GetXaxis()->SetTitle("R (GV)");
        published->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
        published->GetYaxis()->SetTitleOffset(0.55);
        published->GetYaxis()->SetLabelSize(0.058);
        published->GetXaxis()->SetLabelSize(0.06);
        published->GetXaxis()->SetTitleOffset(0.95);
        myflux->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
        myflux->GetXaxis()->SetTitle("R (GV)");
        myflux->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
        myflux->GetYaxis()->SetRangeUser(0.2,upper_range);
        myflux->GetXaxis()->SetRangeUser(2.67,1300);
        published->GetXaxis()->SetRangeUser(2.67,1300);
        published->Draw();
        myflux->Draw("SAME");
        TLegend *a4 = new TLegend(0.1, 0.9, 0.25, 0.75);
        a4->AddEntry(myflux, "My flux");
        a4->AddEntry(published, "Pub");
        a4->SetTextSize(0.035); 
        a4->Draw("SAME");

    c1->cd(2);
        gPad->SetLogx(); // Explicitly set log x for the second pad
        gPad->SetGridy();
        gPad->SetLogy(0); // Ensure log y is off if necessary
        h1->SetMarkerColor(color);
        h1->GetYaxis()->SetTitleOffset(0.55);
        h1->GetYaxis()->SetLabelSize(0.058);
        h1->GetXaxis()->SetLabelSize(0.06);
        h1->GetXaxis()->SetTitleOffset(0.95);
        h1->GetXaxis()->SetRangeUser(2.67,1000);
        h2->GetXaxis()->SetRangeUser(2.67,1000);
        h1->SetTitle("");
        h2->SetTitle("");
        h1->Draw("hist P");
        limits->SetFillColorAlpha(kBlue,0.15);
        limits->Draw("P E3 SAME");
        h2->Draw("hist P SAME");

        TLegend *b4 = new TLegend(0.85, 0.85, 1., 1.);
        b4->AddEntry(h1, "My flux");
        b4->AddEntry(h2, "Pub");
        b4->AddEntry(limits, "Pub tot errors");
        b4->SetTextSize(0.035); 
        b4->Draw("SAME");

    c1->Update();
    c1->SaveAs(out.Data());

    //MIT vs published
    c1->cd(1);
        gPad->SetLogx(); // Use gPad instead of c1
        gPad->SetGridy();
        qyan->SetMarkerColor(color);
        qyan->SetLineColor(color);
        qyan->SetMarkerStyle(20); 
        qyan->SetMarkerSize(0.75);
        qyan->SetTitle("My flux");
        qyan->GetYaxis()->SetRangeUser(0.01,upper_range);
        pubRebin->SetMarkerStyle(20);
        pubRebin->SetMarkerSize(0.75);
        pubRebin->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
        pubRebin->GetXaxis()->SetTitle("R (GV)");
        pubRebin->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
        pubRebin->GetYaxis()->SetTitleOffset(0.55);
        pubRebin->GetYaxis()->SetLabelSize(0.058);
        pubRebin->GetXaxis()->SetLabelSize(0.06);
        pubRebin->GetXaxis()->SetTitleOffset(0.95);
        qyan->SetTitle( Form("%s Flux", getIonName(charge).Data() ) );
        qyan->GetXaxis()->SetTitle("R (GV)");
        qyan->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
        qyan->GetYaxis()->SetRangeUser(0.2,upper_range);
        qyan->GetXaxis()->SetRangeUser(2.67,1300);
        pubRebin->GetXaxis()->SetRangeUser(2.67,1300);
        pubRebin->Draw();
        qyan->Draw("SAME");
        TLegend *a5 = new TLegend(0.1, 0.9, 0.25, 0.75);
        a5->AddEntry(qyan, "MIT");
        a5->AddEntry(pubRebin, "Pub");
        a5->SetTextSize(0.035); 
        a5->Draw("SAME");
    c1->cd(2);
        gPad->SetLogx(); // Explicitly set log x for the second pad
        gPad->SetGridy();
        gPad->SetLogy(0); // Ensure log y is off if necessary
        h11->SetMarkerColor(color);
        h11->GetYaxis()->SetTitleOffset(0.55);
        h11->GetYaxis()->SetLabelSize(0.058);
        h11->GetXaxis()->SetLabelSize(0.06);
        h11->GetXaxis()->SetTitleOffset(0.95);
        h11->GetXaxis()->SetRangeUser(2.67,1000);
        h10->GetXaxis()->SetRangeUser(2.67,1000);
        h11->GetXaxis()->SetTitle("R (GV)");
        h11->SetTitle("");
        h10->SetTitle("");
        h11->GetYaxis()->SetRangeUser(0.7,1.3);
        h10->GetYaxis()->SetRangeUser(0.7,1.3);
        h11->Draw("hist P");
        lim10->SetFillColorAlpha(kBlue,0.15);
        lim10->Draw("P E3 SAME");
        h10->Draw("hist P SAME");

        TLegend *b5 = new TLegend(0.85, 0.85, 1., 1.);
        b5->AddEntry(h11, "MIT");
        b5->AddEntry(h10, "Pub");
        b5->AddEntry(limits, "MIT tot errors");
        b5->SetTextSize(0.035); 
        b5->Draw("SAME");


    c1->SaveAs(out.Data());
    //my vs jack over average
    myflux4->GetXaxis()->SetTitle("R (GV)");
    myflux4->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    jackRebin->GetYaxis()->SetTitle("#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    c1->cd(1);
        gPad->SetLogx(); // Use gPad instead of c1
        gPad->SetGridy();
        jackRebin->SetTitle("");
        jackRebin->GetYaxis()->SetTitleOffset(0.7);
        jackRebin->GetYaxis()->SetLabelSize(0.065);
        jackRebin->GetXaxis()->SetRangeUser(2.9,1000);
        jackRebin->SetMarkerSize(0.75);
        myflux4->GetXaxis()->SetRangeUser(2.9,1000);
        myflux4->SetTitle("");
        myflux4->SetMarkerColor(color);
        myflux4->SetLineColor(color);
        myflux4->SetMarkerStyle(20); 
        myflux4->SetMarkerSize(0.75);
        myflux4->SetTitle("My flux");
        myflux4->GetYaxis()->SetRangeUser(0.01,upper_range);
        jackRebin->Draw();
        myflux4->Draw("SAME");
        TLegend *a6 = new TLegend(0.1, 0.9, 0.25, 0.75);
        a6->AddEntry(myflux4, "My flux");
        a6->AddEntry(jackRebin, "MIT(Jack)");
        a6->SetTextSize(0.035); 
        a6->Draw("SAME");

    c1->cd(2);
        gPad->SetLogx(); // Explicitly set log x for the second pad
        gPad->SetGridy();
        gPad->SetLogy(0); // Ensure log y is off if necessary
        mineOverJack->SetMarkerColor(color);
        mineOverJack->SetMarkerSize(0.75);
        JackOverJack->SetMarkerSize(0.75);
        mineOverJack->GetYaxis()->SetTitleOffset(0.7);
        mineOverJack->GetYaxis()->SetLabelSize(0.058);
        mineOverJack->GetXaxis()->SetLabelSize(0.06);
        mineOverJack->GetXaxis()->SetLabelOffset(0.0001);
        mineOverJack->GetXaxis()->SetTitleOffset(0.95);
        mineOverJack->GetYaxis()->SetRangeUser(0.75,1.35);
        mineOverJack->GetXaxis()->SetRangeUser(2.67,1000);
        mineOverJack->GetXaxis()->SetTitle("R (GV)");
        mineOverJack->GetYaxis()->SetTitle("Flux/Average");
        mineOverJack->SetTitle("");
        JackOverJack->SetTitle("");
        mineOverJack->Draw("hist P");
        errJack->SetFillColorAlpha(kBlue,0.15);
        errJack->Draw("P E3 SAME");
        JackOverJack->Draw("hist P SAME");

        TLegend *b6= new TLegend(0.85, 0.85, 1., 1.);
        b6->AddEntry(mineOverJack, "My flux");
        b6->AddEntry(JackOverJack, "MIT(Jack)");
        b6->AddEntry(errJack, "MIT tot errors");
        b6->SetTextSize(0.035); 
        b6->Draw("SAME");
        c1->Update();
        c1->SaveAs(out.Data());

    c1->SaveAs(Form("%s]", out.Data()));

    return 1;
}*/
void drawTwoPanel(TH1* mine, TH1* other, TH1* otherOverAverage, TH1* mineOverAverage, 
                  TH1* averageContour,  TString titleMine,  TString titleOther, Int_t color, 
                  TPad *pad1, TPad *pad2, TCanvas *c, unsigned int charge) {
    double labelSizeBase = 0.035;
    double titleSizeBase = 0.03;
    double legendSizeBase = 0.03;
    double yTitleOffsetBase = 0.7;
    double xmin=2.15;
    setLogon();
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetOptStat(0);
    // Top pad: flux overlay
    pad1->cd();
    gPad->SetGridy();
	gStyle->SetGridColor(kGray+3);
    gPad->SetLogx();
    mine->GetXaxis()->SetRangeUser(xmin,1000);
    mine->SetMarkerColor(color); mine->SetLineColor(color);
    mine->SetMarkerStyle(20); mine->SetMarkerSize(0.75);
    mine->SetTitle(titleMine);
    mine->GetXaxis()->SetTitle("R (GV)");
    mine->GetYaxis()->SetTitle("#phi R^{2.7} (m^{-2}s^{-1}sr^{-1}GV^{1.7})");
    mine->GetYaxis()->SetTitleOffset(1.);
    mine->GetXaxis()->SetLabelSize(0);   // Hides X axis tick labels
    mine->GetXaxis()->SetTitleSize(0);   // Hides X axis title
    mine->GetYaxis()->CenterTitle(true);
    mine->SetLabelFont(62,"XYZ");
    mine->SetTitleFont(62,"XYZ");
    mine->GetYaxis()->SetLabelSize(labelSizeBase/0.7);
    mine->GetYaxis()->SetTitleSize(titleSizeBase/0.7);
    mine->Draw();
    other->SetMarkerStyle(20); other->SetMarkerSize(0.75); other->SetLineColor(kBlack);
    mine->SetTitle(titleMine);
    other->GetXaxis()->SetRangeUser(xmin,1000);
    other->Draw("SAME");
    TLegend* leg1 = new TLegend(0.1, 0.9, 0.25, 0.75);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);   
    leg1->SetTextSize(legendSizeBase/0.7);
    leg1->SetTextFont(62);
    leg1->AddEntry(mine, "Perugia (A. Ubaldi)");
    leg1->AddEntry(other, titleOther);
    leg1->Draw("SAME");

    // Bottom pad: flux / average
    pad2->cd();
    gPad->SetGridy();
    mineOverAverage->SetMarkerColor(color); mineOverAverage->SetLineColor(color);
    mineOverAverage->SetMarkerStyle(20); mineOverAverage->SetMarkerSize(0.75);
    otherOverAverage->SetLineColor(kBlack);
    otherOverAverage->SetMarkerStyle(20); otherOverAverage->SetMarkerSize(0.75);
    averageContour->SetLabelFont(62,"XYZ");
    averageContour->SetTitleFont(62,"XYZ");
    averageContour->GetYaxis()->SetTitle("Flux/Average");
    averageContour->GetXaxis()->SetTitle("R (GV)");
    averageContour->GetXaxis()->CenterTitle(true);
    averageContour->GetYaxis()->CenterTitle(true);
    averageContour->GetYaxis()->SetTitleOffset(0.42);
    averageContour->GetYaxis()->SetLabelSize(labelSizeBase/0.3);
    averageContour->GetXaxis()->SetLabelSize(labelSizeBase/0.3);
    averageContour->GetYaxis()->SetTitleSize(titleSizeBase/0.3);
    averageContour->GetXaxis()->SetTitleSize(titleSizeBase/0.3);

    averageContour->SetFillColorAlpha(kBlue, 0.15);
    averageContour->GetXaxis()->SetRangeUser(xmin,1000);
    auto rng = fluxRange(charge); 
    averageContour->GetYaxis()->SetRangeUser(rng.first,rng.second);
    mineOverAverage->GetXaxis()->SetRangeUser(xmin,1000);
    otherOverAverage->GetXaxis()->SetRangeUser(xmin,1000);
    gPad->SetLogx();
    gPad->SetGridy();
	gStyle->SetGridColor(kGray+3);
    averageContour->Draw("P E3");
    mineOverAverage->Draw("hist P SAME");
    otherOverAverage->Draw("hist P SAME");
    TLegend* leg2 = new TLegend(0.1, 0.9, 0.65, 0.85);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);   
    leg2->SetTextSize(legendSizeBase/0.3);
    leg2->SetTextFont(62);
    leg2->SetNColumns(3);
    leg2->AddEntry(mineOverAverage, "PERUGIA");
    leg2->AddEntry(otherOverAverage, titleOther);
    leg2->AddEntry(averageContour, "Total errors");
    leg2->Draw("SAME");

    // Save output
    TString outName = Form("../output/Comparison/%s.pdf", titleMine.Data());
    c->SaveAs(outName);
}
TH1D *DiffOverRatio( TH1D *pub,  TH1D *my) {
    TH1D *h = new TH1D();
    for (int ii=0; ii<pub->GetNbinsX(); ii++) {
        h->SetBinContent(ii, (pub->GetBinContent(ii) - my->GetBinContent(ii)) /my->GetBinError(ii));
    }
    return h;
}
short int GetFluxColor(unsigned int charge) {
    short int c;
    switch (charge) {
        case 14:
            c = kGreen-2;
            break;
        case 15:
            c = kViolet;
            break;      
        case 16:
            c = kOrange+2;
            break;
        case 18:
            c = kRed+1;
            break;
        case 20:
            c = kAzure+2;
            break;
    }
    return c;
}
std::pair<double,double> fluxRange(unsigned int charge) {
    double ymin,ymax;
    switch(charge) {
        case 14:
            ymin=0.94;
            ymax=1.06;
            break;
        case 15:
            ymin=0.85;
            ymax=1.15;
            break;
        case 16:
            ymin=0.85;
            ymax=1.15;
            break;
        case 18:
            ymin=0.85;
            ymax=1.15;
            break;
        case 20:
            ymin=0.85;
            ymax=1.15;
            break;
    }
    return std::make_pair(ymin,ymax);
}
