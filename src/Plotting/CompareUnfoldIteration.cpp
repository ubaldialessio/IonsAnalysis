#include "binning.h"
#include "utils.h"
#include <TParameter.h>

TFile *getFile(unsigned int charge, TString timePeriod);
std::vector<TH1D*> getHistogram(TFile *f, TString name, int iter);
std::vector<TF1*> getHistogramFit(TFile *f, TString name, int iter);
std::vector<TH1D*> getAcc(TFile *f, TString name, int iter);
void print(std::vector<TH1D*> v, unsigned int charge, TH1D *pub, TString timePeriod);
void printRelDiff(std::vector<TH1D*> v, unsigned int charge, TH1D *pub, TString timePeriod);
void print_fit_model(std::vector<TH1D*> flux, std::vector<TF1*> flux_model, unsigned int charge, TString timePeriod);
void print_acc_fit(std::vector<TH1D*> flux, std::vector<TH1D*> flux_model, unsigned int charge, TString timePeriod);
TH1D* ComputePull(TH1* hData, TF1* fitFunc, const std::string& pullName = "pullHist");

//-------------------------MAIN----------------------------
// Step 1. compare, for every iteration, the flux
int main(int argc, char **argv) {
	if (argc < 3) {
		printf("Usage: \n");
		printf("%s <time> <charge> \n", argv[0]);
		return 1;
	}
    TString timePeriod = argv[1];
    TString exename = argv[0];
    unsigned int charge = atoi(argv[2]);
    auto file = getFile(charge, timePeriod);
    int iter = 50;
    TFile *f= new TFile("../IonsSelected/"+getIonPath(charge)+"/Unfolding/"+getIonPath(charge)+
                        Form("_UnfoldingFactor_%s.root", timePeriod.Data()) );
    auto pub = (TH1D*)f->Get("pub_flux");
    MultiplyByXPower(pub,2.7);
    auto flux = getHistogram(file,"flux",iter);
    auto flux_model = getHistogramFit(file,"flux_model",iter);
    auto acc = getAcc(file,"acc",iter);
    auto spline_acc = getAcc(file,"spline_acc",iter);
    print(flux, charge,pub, timePeriod);
    //flux.insert(flux.begin(),pub);
    std::cout << flux.size() << std::endl;
    std::cout << flux_model.size() << std::endl;
    print_fit_model(flux, flux_model,charge,timePeriod);
    print_acc_fit(acc, spline_acc, charge, timePeriod);
    printRelDiff(flux, charge,pub, timePeriod);
}

TFile *getFile(unsigned int charge, TString timePeriod) {
    TString ionPath = getIonPath(charge);
	TFile *file = new TFile("../IonsSelected/"+ionPath+"/Unfolding/"+ionPath+Form("_UnfoldingFactor_%s.root",timePeriod.Data()) );
    if (file && !file->IsZombie()) {
        return file;
    } else std::cerr << "Error" << std::endl;
}
std::vector<TH1D*> getHistogram(TFile *f, TString name, int iter) {
    std::vector<TH1D*> v;
    for (int i=0; i<=iter; i++) {
        auto h = (TH1D*)f->Get(Form("%s_it%i",name.Data(),i) );
        if (h) {
            std::cout << "OK for iter " << i << std::endl;
            MultiplyByXPower(h,2.7);
            v.push_back(h);
        } 
    }
    return v;
}
std::vector<TF1*> getHistogramFit(TFile *f, TString name, int iter) {
    std::vector<TF1*> v;
    for (int i=0; i<=iter; i++) {
        auto h = (TF1*)f->Get(Form("%s_it%i",name.Data(),i) );
        if (h) {
            std::cout << "OK for iter " << i << std::endl;
            auto vv = MultiplyTF1ByXPower(h, 2.7);
                double chi2 = h->GetChisquare();
                int npoints = h->GetNumberFitPoints();
                int npar = h->GetNumberFreeParameters();
                int ndf = npoints - npar;
                auto title = Form("Chi2/ndf: %f", chi2 / ndf);
                vv->SetTitle(title);
                vv->SetName(title);
            v.push_back(vv);
        } 
    }
    return v;   
}
std::vector<TH1D*> getAcc(TFile *f, TString name, int iter) {
    std::vector<TH1D*> v;
    for (int i=0; i<=iter; i++) {
        auto h = (TH1D*)f->Get(Form("%s_it%i",name.Data(),i) );
        if (h) {
            std::cout << "OK for iter " << i << std::endl;
            v.push_back(h);
        } 
    }
    return v;   
}
void print(std::vector<TH1D*> v, unsigned int charge, TH1D *pub, TString timePeriod) {
    TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TCanvas *b = new TCanvas("","",2048,1280);
    TString prefix = getIonPath(charge);
	b->SaveAs( Form("../output/%s_unf_%s.pdf[", prefix.Data(),timePeriod.Data()) ,"RECREATE");		
    auto a = new TLegend(0.7, 0.80, 0.8, 0.9);
    //MultiplyByXPower(pub,2.7);
    pub->SetMarkerSize(1.6);
    pub->SetMarkerStyle(20);
    getColor(pub,charge,"dat");
    double max_rel_diff=0;
    for (int i=0; i<v.size(); i++) {
        a->Clear();
        a->AddEntry(pub,"Pub");
        b->Update();
        formatAxis(v[i], 1);		
		formatTitle(v[i],1);				
        setTitle(v[i], "Flux", "R (GV)", "#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})", charge);
        //MultiplyByXPower(v[i],2.7);
		adjustZoomX(v[i],"Flux");
        v[i]->GetYaxis()->SetLabelFont(62);
        v[i]->GetXaxis()->SetLabelFont(62);
        v[i]->GetYaxis()->SetTitleFont(62);
        v[i]->GetXaxis()->SetTitleFont(62);
        gPad->SetGridy();
        gStyle->SetGridColor(kGray);
        gStyle->SetGridStyle(1);
        gStyle->SetGridWidth(1);
        getColor(v[i], charge, "dat");
        setFriendColor(v[i],i);	
        v[i]->SetMarkerSize(1.6);
        if (i==0)
            a->AddEntry(v[i],Form("It %i",i+1));
        v[i]->Draw("SAME");
        gPad->SetLogx();
        //gPad->SetLogy();
        pub->Draw("SAME");
        if (i>0) {
            //auto relativeDiffHist = relDiff(v[i],v[i-1],0.9957, 2001);
            max_rel_diff = findMaxRelDiff(v[i],v[i-1],0.9957, 2001);
            a->AddEntry(v[i],Form("It %i - Diff %.2f",i+1, max_rel_diff) );
            //relativeDiffHist->Draw("SAME");
        }
        a->Draw();
        b->SaveAs( Form("../output/%s_unf_%s.pdf", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
    }
    b->SaveAs( Form("../output/%s_unf_%s.pdf]", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
}
void printRelDiff(std::vector<TH1D*> v, unsigned int charge, TH1D *pub, TString timePeriod) {
 TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TCanvas *b = new TCanvas("","",2048,1280);
    TString prefix = getIonPath(charge);
	b->SaveAs( Form("../output/%s_unfRel_%s.pdf[", prefix.Data(),timePeriod.Data()) ,"RECREATE");		
    auto a = new TLegend(0.7, 0.80, 0.8, 0.9);
    //MultiplyByXPower(pub,2.7);
    pub->SetMarkerSize(1.6);
    pub->SetMarkerStyle(20);
    getColor(pub,charge,"dat");
    double max_rel_diff=0;
    for (int i=0; i<v.size(); i++) {
        a->Clear();
        a->AddEntry(pub,"Pub");
        b->Update();
        //MultiplyByXPower(v[i],2.7);
		adjustZoomX(v[i],"Flux");
        gPad->SetGridy();
        gStyle->SetGridColor(kGray);
        gStyle->SetGridStyle(1);
        gStyle->SetGridWidth(1);
        setFriendColor(v[i],i);	
        v[i]->SetMarkerSize(1.6);
        gPad->SetLogx();
        if (i>0) {
            auto relativeDiffHist = relDiff(v[i],v[i-1],0.9957, 2001);
            relativeDiffHist->GetYaxis()->SetLabelFont(62);
            relativeDiffHist->GetXaxis()->SetLabelFont(62);
            relativeDiffHist->GetYaxis()->SetTitleFont(62);
            relativeDiffHist->GetXaxis()->SetTitleFont(62);
            formatAxis(relativeDiffHist, 1);		
            formatTitle(relativeDiffHist,1);				
            setTitle(relativeDiffHist, "Flux", "R (GV)", "Relative Difference %", charge);
            getColor(relativeDiffHist, charge, "dat");
            relativeDiffHist->GetYaxis()->SetTitleOffset(1.22);
            max_rel_diff = findMaxRelDiff(v[i],v[i-1],0.9957, 2001);
            a->AddEntry(v[i],Form("It %i - Diff %.2f",i+1, max_rel_diff) );
            relativeDiffHist->Draw();
        }
        b->SaveAs( Form("../output/%s_unfRel_%s.pdf", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
    }
    b->SaveAs( Form("../output/%s_unfRel_%s.pdf]", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
}
void print_fit_model(std::vector<TH1D*> flux, std::vector<TF1*> flux_model, unsigned int charge, TString timePeriod) {
    TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TCanvas *b = new TCanvas("","",2048,1280);
    TString prefix = getIonPath(charge);
	b->SaveAs( Form("../output/%s_fitUnf_%s.pdf[", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
    auto a = new TLegend(0.1, 0.9, 0.35, 0.75);
    a->SetBorderSize(0);
    a->SetFillStyle(0);  
    auto pads = buildPadsForFitAndPull(0.01);
    pads.second->Draw();
    pads.first->Draw("SAME");
    double labelSizeBase = 0.035;
    double titleSizeBase = 0.03;
    double legendSizeBase = 0.03;

    for (int i=0; i<flux.size(); i++) {

        //top pad    
        pads.first->cd();
        a->Clear();
        b->Update();
        setTitle(flux[i], "Flux", "R (GV)", "#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})", charge);
        formatAxis(flux[i], 1);		
		formatTitle(flux[i],1);	
        adjustZoomX(flux[i],"Flux");
        getColor(flux[i], charge, "dat");
        pads.first->SetGridy();
        gStyle->SetGridColor(kGray+3);
        flux[i]->GetXaxis()->CenterTitle(true);
        flux[i]->GetYaxis()->CenterTitle(true);
        flux[i]->GetXaxis()->SetRangeUser(2.,1000);
        flux[i]->GetYaxis()->SetLabelSize(labelSizeBase/0.7);
        flux[i]->GetYaxis()->SetTitleSize(titleSizeBase/0.7);
        //flux[i]->GetYaxis()->SetTitleOffset(0.25);
        flux[i]->GetXaxis()->SetLabelSize(labelSizeBase/0.7);
        flux[i]->GetXaxis()->SetTitleSize(titleSizeBase/0.7);
        //flux[i]->GetXaxis()->SetTitleOffset(1.04);
        //flux[i]->GetXaxis()->SetLabelSize(0);   // Hides X axis tick labels
        //flux[i]->GetXaxis()->SetTitleSize(0);
        //flux_model[i]->GetXaxis()->SetLabelSize(0);   // Hides X axis tick labels
        //flux_model[i]->GetXaxis()->SetTitleSize(0);
        flux[i]->Draw();
        flux_model[i]->Draw("SAME");
        a->AddEntry(flux[i],Form("Flux_it_%i",i) );
        //chi2
        const char* input_cstr = flux_model[i]->GetName();
        std::cout << flux_model[i]->GetName() << std::endl;
        TString input(input_cstr);
        TString result;
        int chi2Pos = input.Index("Chi2/ndf:");
        if (chi2Pos != kNPOS) {
            TString chi2Part = input(chi2Pos + 9, input.Length() - chi2Pos - 9);
            chi2Part.Strip(TString::kBoth);  // Remove leading/trailing spaces
            double chi2Value = chi2Part.Atof();
            result.Form("Spline: #chi^{2}/ndf = %.3f", chi2Value);
            a->AddEntry(flux_model[i], Form("Fit: Chi2/ndf %.3f", chi2Value));
        }

        a->SetTextSize(0.036);
        a->Draw();
        gPad->SetLogx();

        //bottom pad, pull
        auto pull_temp = ComputePull(flux[i],flux_model[i]);
        TH1D* pull = new TH1D("pull", "; Pull; Entries", 6, -3, 3);
        for (int i = 1; i <= pull_temp->GetNbinsX(); i++) {
            pull->Fill(pull_temp->GetBinContent(i));
        }

        pads.second->cd();
        pads.second->SetGridy();
        formatAxis(pull, 1);		
        pull->SetTitle("");
        pull->GetXaxis()->CenterTitle(true);
        pull->GetYaxis()->CenterTitle(true);
        pull->GetYaxis()->SetLabelSize(labelSizeBase/0.3);
        pull->GetYaxis()->SetTitleSize(titleSizeBase/0.3);
        pull->GetYaxis()->SetTitleOffset(0.45);
        pull->GetXaxis()->SetLabelSize(labelSizeBase/0.3);
        pull->GetXaxis()->SetTitleSize(titleSizeBase/0.3);
        pull->GetXaxis()->SetLabelOffset(0.005);
        //pull->GetXaxis()->SetTitleOffset(1.04);
        pull->SetMarkerSize(0.8);
        pull->SetMarkerColor(kBlack);
        pull->SetLineColor(kBlack);

        pull->GetYaxis()->SetNdivisions(505);
        //pull->GetXaxis()->SetNdivisions(505);

        pull->SetMarkerStyle(20);
        pull->SetMarkerSize(1.8);

        // Get axis limits from pull
        double xMin = pull->GetXaxis()->GetXmin();
        double xMax = pull->GetXaxis()->GetXmax();

        TF1* gausFit = new TF1("gausFit", "gaus", -3, 3);
        pull->Fit(gausFit, "Q");

        double mean = gausFit->GetParameter(1);
        double sigma = gausFit->GetParameter(2);
        double mean_err = gausFit->GetParError(1);
        double sigma_err = gausFit->GetParError(2);

        TLegend* legPull = new TLegend(0.05, 0.65, 0.4, 0.9);  // (x1, y1, x2, y2) in NDC
        legPull->SetBorderSize(0);
        legPull->SetFillStyle(0);  // Transparent background
        legPull->SetTextSize(0.1);  // Smaller text for small pad
        legPull->AddEntry((TObject*)nullptr, Form("#mu = %.3f",
                                            mean),"" );
        legPull->AddEntry((TObject*)nullptr, Form("#sigma = %.3f",sigma),"");                                   
        legPull->SetTextFont(62);

        pull->GetYaxis()->SetRangeUser(0.,58.);

        TGraph* band = new TGraph(4);
        band->SetPoint(0, xMin, -1);
        band->SetPoint(1, xMax, -1);
        band->SetPoint(2, xMax, 1);
        band->SetPoint(3, xMin, 1);
        band->SetFillColorAlpha(kGray+1, 0.3);
        band->SetLineColor(0);
        band->SetFillStyle(1001);
        band->SetName("sigmaBand");    
        pull->Draw("hist");
        band->Draw("F SAME");     
        gausFit->Draw("SAME");
        legPull->Draw();
        b->SaveAs( Form("../output/%s_fitUnf_%s.pdf", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
    }
    b->SaveAs( Form("../output/%s_fitUnf_%s.pdf]", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
}
void print_acc_fit(std::vector<TH1D*> flux, std::vector<TH1D*> flux_model, unsigned int charge, TString timePeriod) {
    TStyle *style = effHistStyle(0);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	setLogon();
	style->cd();
	TCanvas *b = new TCanvas("","",2048,1280);
    TString prefix = getIonPath(charge);
	b->SaveAs( Form("../output/%s_fitAcc_%s.pdf[", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
    auto a = new TLegend(0.1, 0.9, 0.35, 0.75);
    auto pads = buildPadsForFitAndPull(0.01);
    pads.second->Draw();
    pads.first->Draw("SAME");
    double labelSizeBase = 0.035;
    double titleSizeBase = 0.03;
    double legendSizeBase = 0.03;

    for (int i=0; i<flux.size(); i++) {

        //top pad    
        pads.first->cd();
        a->Clear();
        b->Update();
        setTitle(flux[i], "Flux", "R (GV)", "#phi R^{2.7}(m^{-2}s^{-1}sr^{-1}GV^{1.7})", charge);
        formatAxis(flux[i], 1);		
		formatTitle(flux[i],1);	
        adjustZoomX(flux[i],"Flux");
        getColor(flux[i], charge, "dat");
        pads.first->SetGridy();
        gStyle->SetGridColor(kGray+3);

        flux[i]->Rebin(9);
        flux[i]->Scale(1./9.);

        flux[i]->GetXaxis()->CenterTitle(true);
        flux[i]->GetYaxis()->CenterTitle(true);
        flux[i]->GetXaxis()->SetRangeUser(1.9,2000);
        flux[i]->GetYaxis()->SetLabelSize(labelSizeBase/0.7);
        flux[i]->GetYaxis()->SetTitleSize(titleSizeBase/0.7);
        //flux[i]->GetYaxis()->SetTitleOffset(0.25);
        flux[i]->GetXaxis()->SetLabelSize(labelSizeBase/0.7);
        flux[i]->GetXaxis()->SetTitleSize(titleSizeBase/0.7);
        //flux[i]->GetXaxis()->SetTitleOffset(1.04);
        //flux[i]->GetXaxis()->SetLabelSize(0);   // Hides X axis tick labels
        //flux[i]->GetXaxis()->SetTitleSize(0);
        //flux_model[i]->GetXaxis()->SetLabelSize(0);   // Hides X axis tick labels
        //flux_model[i]->GetXaxis()->SetTitleSize(0);
        flux[i]->Draw();
        getColor(flux_model[i], charge, "mc");
        flux_model[i]->SetFillStyle(0);
		flux_model[i]->DrawCopy("hist L SAME");
		flux_model[i]->SetFillStyle(1001);
		flux_model[i]->SetMarkerSize(0);
		flux_model[i]->Draw("E3 L SAME");
        //flux_model[i]->Draw("SAME");
        a->AddEntry(flux[i],Form("Acc_it_%i",i) );
        //chi2
        const char* input_cstr = flux_model[i]->GetName();
        TString input(input_cstr);
        TString result;
        int chi2Pos = input.Index("Chi2/ndf:");
        if (chi2Pos != kNPOS) {
            TString chi2Part = input(chi2Pos + 9, input.Length() - chi2Pos - 9);
            chi2Part.Strip(TString::kBoth);  // Remove leading/trailing spaces
            double chi2Value = chi2Part.Atof();
            result.Form("Spline: #chi^{2}/ndf = %.3f", chi2Value);
            a->AddEntry(flux_model[i], Form("Fit: Chi2/ndf %.3f", chi2Value));
        }
        a->SetBorderSize(0);
 		a->SetFillStyle(0);  // Transparent background
        a->SetTextSize(0.036);
        a->Draw();
        gPad->SetLogx();

        //bottom pad, pull
        /*auto pull_temp = ComputePull(flux[i],flux_model[i]);
        TH1D* pull = new TH1D("pull", "; Pull; Entries", 6, -3, 3);
        for (int i = 1; i <= pull_temp->GetNbinsX(); i++) {
            pull->Fill(pull_temp->GetBinContent(i));
        }

        pads.second->cd();
        pads.second->SetGridy();
        formatAxis(pull, 1);		
        pull->SetTitle("");
        pull->GetXaxis()->CenterTitle(true);
        pull->GetYaxis()->CenterTitle(true);
        pull->GetYaxis()->SetLabelSize(labelSizeBase/0.3);
        pull->GetYaxis()->SetTitleSize(titleSizeBase/0.3);
        pull->GetYaxis()->SetTitleOffset(0.45);
        pull->GetXaxis()->SetLabelSize(labelSizeBase/0.3);
        pull->GetXaxis()->SetTitleSize(titleSizeBase/0.3);
        pull->GetXaxis()->SetLabelOffset(0.005);
        //pull->GetXaxis()->SetTitleOffset(1.04);
        pull->SetMarkerSize(0.8);
        pull->SetMarkerColor(kBlack);
        pull->SetLineColor(kBlack);

        pull->GetYaxis()->SetNdivisions(505);
        //pull->GetXaxis()->SetNdivisions(505);

        pull->SetMarkerStyle(20);
        pull->SetMarkerSize(1.8);

        // Get axis limits from pull
        double xMin = pull->GetXaxis()->GetXmin();
        double xMax = pull->GetXaxis()->GetXmax();

        TF1* gausFit = new TF1("gausFit", "gaus", -3, 3);
        pull->Fit(gausFit, "Q");

        double mean = gausFit->GetParameter(1);
        double sigma = gausFit->GetParameter(2);
        double mean_err = gausFit->GetParError(1);
        double sigma_err = gausFit->GetParError(2);

        TLegend* legPull = new TLegend(0.05, 0.65, 0.4, 0.9);  // (x1, y1, x2, y2) in NDC
        legPull->SetBorderSize(0);
        legPull->SetFillStyle(0);  // Transparent background
        legPull->SetTextSize(0.1);  // Smaller text for small pad
        legPull->AddEntry((TObject*)nullptr, Form("#mu = %.3f",
                                            mean),"" );
        legPull->AddEntry((TObject*)nullptr, Form("#sigma = %.3f",sigma),"");                                   
        legPull->SetTextFont(62);

        pull->GetYaxis()->SetRangeUser(0.,58.);

        TGraph* band = new TGraph(4);
        band->SetPoint(0, xMin, -1);
        band->SetPoint(1, xMax, -1);
        band->SetPoint(2, xMax, 1);
        band->SetPoint(3, xMin, 1);
        band->SetFillColorAlpha(kGray+1, 0.3);
        band->SetLineColor(0);
        band->SetFillStyle(1001);
        band->SetName("sigmaBand");    
        pull->Draw("hist");
        band->Draw("F SAME");     
        gausFit->Draw("SAME");
        legPull->Draw();*/
        b->SaveAs( Form("../output/%s_fitAcc_%s.pdf", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
    }
    b->SaveAs( Form("../output/%s_fitAcc_%s.pdf]", prefix.Data(),timePeriod.Data()) ,"RECREATE");	
}
TH1D* ComputePull(TH1* hData, TF1* fitFunc, const std::string& pullName) {
    if (!hData || !fitFunc) return nullptr;

    int nbins = hData->GetNbinsX();
    auto pullHist = (TH1D*)hData->Clone(pullName.c_str());
    pullHist->Reset("ICES"); // Clears contents but keeps axis etc.

    for (int i = 1; i <= nbins; ++i) {
        double x = hData->GetBinCenter(i);
        double y_data = hData->GetBinContent(i);
        double y_fit = fitFunc->Eval(x);
        double err = hData->GetBinError(i);

        double pull = (err > 0) ? (y_data - y_fit) / err : 0;

        pullHist->SetBinContent(i, pull);

        printf("For x = %.2f: pull = %.2f\n",x, pull);
    }

    pullHist->SetTitle("Pull Distribution;Bin center;Pull");
    return pullHist;
}
TH1D* ComputePull(TH1* hData, TH1* fitFunc, const std::string& pullName) {
    if (!hData || !fitFunc) return nullptr;

    int nbins = hData->GetNbinsX();
    auto pullHist = (TH1D*)hData->Clone(pullName.c_str());
    pullHist->Reset("ICES"); // Clears contents but keeps axis etc.

    for (int i = 1; i <= nbins; ++i) {
        double x = hData->GetBinCenter(i);
        double y_data = hData->GetBinContent(i);
        double y_fit = fitFunc->GetBinCenter(i);
        double err = hData->GetBinError(i);

        double pull = (err > 0) ? (y_data - y_fit) / err : 0;

        pullHist->SetBinContent(i, pull);

        printf("For x = %.2f: pull = %.2f\n",x, pull);
    }

    pullHist->SetTitle("Pull Distribution;Bin center;Pull");
    return pullHist;
}
