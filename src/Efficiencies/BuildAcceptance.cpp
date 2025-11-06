#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "SplineUtility.h"

TH1D *finalAcceptance(unsigned int charge);
TH1D* GetEff(TH1* hpass, TH1* htotal, TString option);

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: \n");
	    printf("%s <charge>  \n", argv[0]);
	    return 1;
	}
    unsigned int charge=atoi(argv[1]);
    //auto final_acc_7 = finalAcceptance(charge);

    TString ionPath = getIonPath(charge);
    TString inpPath = "../IonsSelected/"+ionPath+"/"+ionPath+"_output.root";
    TFile *file = TFile::Open(inpPath.Data());
    TH1D * damc_l1_7 = (TH1D*)file->Get("damc_l1_7");
    TH1D * damc_tf_7 = (TH1D*)file->Get("damc_tf_7");
    TH1D * damc_tr_7 = (TH1D*)file->Get("damc_tr_7");
    TH1D * damc_tk_7 = (TH1D*)file->Get("damc_tk_7");
		auto final_damc_l1 = (TH1D*)file->Get("final_damc_l1_7");
		auto final_damc_tf = (TH1D*)file->Get("final_damc_tf_7");
		auto final_damc_tr = (TH1D*)file->Get("final_damc_tr_7");
		auto final_damc_tk = (TH1D*)file->Get("final_damc_tk_7");
    auto final_damc_tot = (TH1D*)file->Get("final_damc_tot_7");
	  TString mcPath = "../IonsSelected/"+ionPath+"/mc/mc.root";
    auto mc = new TFile(mcPath.Data() );
    TH1D *mc_pass_gen  = (TH1D *)mc->Get("mc_pass_gen");
    TH1D *mc_samp	   = (TH1D *)mc->Get("mc_samp");
    TH1D *mc_pass      = (TH1D *)mc->Get("mc_pass");
    mc_samp->Rebin(9);
    mc_pass->Rebin(9);
    mc_pass_gen->Rebin(9);
    auto true_acc  = GetEff(mc_pass_gen, mc_samp, "simple");
    true_acc->Scale(TMath::Pi() * 3.9 * 3.9, "nosw2");
    auto acc = (TH1D*)true_acc->Clone();
    acc = GetEff(mc_pass, mc_samp, "simple");
	  acc->Scale(TMath::Pi() * 3.9 * 3.9, "nosw2");
	  auto spline_acc = autospline(acc, 1, 500);
    auto final_acc = new TH1D();
    final_acc = (TH1D*)spline_acc->Clone();
    final_acc->Multiply(final_damc_tot);
    file->Close();

    TString output= "../IonsSelected/"+getIonPath(charge)+"/"+getIonPath(charge)+"_output.root";
    auto outfile = new TFile(output.Data() ,"UPDATE");
    outfile->Delete("final_acc;*");
    outfile->WriteTObject(final_acc,"final_acc_7");
    return 1;
}

TH1D *finalAcceptance(unsigned int charge) {
    TString ionPath = getIonPath(charge);
    TString inpPath = "../IonsSelected/"+ionPath+"/"+ionPath+"_output.root";
    TFile *file = TFile::Open(inpPath.Data());
    TH1D * damc_l1_7 = (TH1D*)file->Get("damc_l1_7");
    TH1D * damc_tf_7 = (TH1D*)file->Get("damc_tf_7");
    TH1D * damc_tr_7 = (TH1D*)file->Get("damc_tr_7");
    TH1D * damc_tk_7 = (TH1D*)file->Get("damc_tk_7");
		auto final_damc_l1 = autospline(damc_l1_7, 0.5, 30);
		auto final_damc_tf = autospline(damc_tf_7, 0.5, 15);
		auto final_damc_tr = autospline(damc_tr_7, 0.5, 65);
		auto final_damc_tk = autospline(damc_tk_7, 0.5, 7);
    auto final_damc_tot = new TH1D();
    final_damc_tot = (TH1D*)final_damc_l1->Clone();
    final_damc_tot->Multiply(final_damc_tf);
    final_damc_tot->Multiply(final_damc_tr);
    final_damc_tot->Multiply(final_damc_tk);
	  TString mcPath = "../IonsSelected/"+ionPath+"/mc/mc.root";
    auto mc = new TFile(mcPath.Data() );
    TH1D *mc_pass_gen  = (TH1D *)mc->Get("mc_pass_gen");
    TH1D *mc_samp	   = (TH1D *)mc->Get("mc_samp");
    TH1D *mc_pass      = (TH1D *)mc->Get("mc_pass");
    mc_samp->Rebin(9);
    mc_pass->Rebin(9);
    mc_pass_gen->Rebin(9);
    auto true_acc  = GetEff(mc_pass_gen, mc_samp, "simple");
    true_acc->Scale(TMath::Pi() * 3.9 * 3.9, "nosw2");
    auto acc = (TH1D*)true_acc->Clone();
    acc = GetEff(mc_pass, mc_samp, "simple");
	  acc->Scale(TMath::Pi() * 3.9 * 3.9, "nosw2");
	  auto spline_acc = autospline(acc, 1, 500);
    auto final_acc = new TH1D();
    final_acc = (TH1D*)spline_acc->Clone();
    final_acc->Multiply(final_damc_tot);
    file->Close();
    return final_acc;
}

TH1D* GetEff(TH1* hpass, TH1* htotal, TString option) {
  hpass->Sumw2();
  htotal->Sumw2();
  auto ratio = new TH1D("", "", nRbins_HighZ-1, Rbins_HighZ);
  if(option.Contains("efficiency")){
    auto eff = new TEfficiency(*hpass, *htotal);
    eff->SetStatisticOption(TEfficiency::kBJeffrey);
    for(int i=1; i<ratio->GetNbinsX(); ++i){
	  ratio->SetBinContent(i, eff->GetEfficiency(i));
      double hierr = eff->GetEfficiencyErrorUp(i);
      double loerr = eff->GetEfficiencyErrorLow(i);
      double err = hierr>loerr ? hierr : loerr;
      ratio->SetBinError(i, err);
    }
  }
  else if(option.Contains("correction")){
    for (int i=1; i<=hpass->GetNbinsX(); ++i){
      if (htotal->GetBinContent(i)!=0) {
      	ratio->SetBinContent(i, hpass->GetBinContent(i)/htotal->GetBinContent(i));
      	double err = sqrt(hpass->GetBinError(i)*hpass->GetBinError(i) + htotal->GetBinError(i)*htotal->GetBinError(i));
      	ratio->SetBinError(i, err); 
      }
    }
  }
  else if(option.Contains("simple")){
    ratio = (TH1D*)hpass->Clone();
    ratio->Divide(htotal);
  }
  return ratio;
}