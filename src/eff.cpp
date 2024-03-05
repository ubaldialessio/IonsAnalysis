#include "selection.h"

TH1D* GetEff(TH1* hpass, TH1* htotal, TString option) {
  hpass->Sumw2();
  htotal->Sumw2();
  auto ratio = new TH1D("", "", nRbins-1, Rbins);
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
      ratio->SetBinContent(i, hpass->GetBinContent(i)/htotal->GetBinContent(i));
      double err = sqrt(hpass->GetBinError(i)*hpass->GetBinError(i) + htotal->GetBinError(i)*htotal->GetBinError(i));
      ratio->SetBinError(i, err);
    }
  }
  else if(option.Contains("simple")){
    ratio = (TH1D*)hpass->Clone();
    ratio->Divide(htotal);
  }
  return ratio;
}

int main(int argc, char **argv) {
		if (argc < 2) {
			printf("Usage: n");
			printf("%s <input.root> <output.root> n", argv[0]);
			return 1;
		}
	TString inp = argv[1], out = argv[2];
	TFile *f = new TFile(inp.Data() );
	auto outfile = new TFile(out.Data() ,"recreate");
	TH2F *h1 = (TH2F *)f->Get("pass_tof");
	TH2F *h2 = (TH2F *)f->Get("sample_tof");
	TH2F *h3 = (TH2F *)f->Get("pass_l1");
	TH2F *h4 = (TH2F *)f->Get("sample_l1");
	TH2F *h5 = (TH2F *)f->Get("pass_tr");
	TH2F *h6 = (TH2F *)f->Get("sample_tr");
	TH2F *h7 = (TH2F *)f->Get("sample_trig");
	TH2F *h8 = (TH2F *)f->Get("phys_trig");
	TH2F *h9 = (TH2F *)f->Get("unbiased_trig");
	
	TH1D *pass_tof	 = h1->ProjectionY("a",start_bin,stop_bin);
	TH1D *sample_tof = h2->ProjectionY("b",start_bin,stop_bin);
	TH1D *pass_l1	 = h3->ProjectionY("c",start_bin,stop_bin);
	TH1D *sample_l1  = h4->ProjectionY("d",start_bin,stop_bin);
	TH1D *pass_tr    = h5->ProjectionY("e",start_bin,stop_bin);
	TH1D *sample_tr  = h6->ProjectionY("f",start_bin,stop_bin);
	TH1D *sample_trig= h7->ProjectionY("g",start_bin,stop_bin);
	TH1D *phys_trig  = h8->ProjectionY("h",start_bin,stop_bin);
	TH1D *unb_trig   = h9->ProjectionY("i",start_bin,stop_bin);
	TH1D *mc_pass_gen= (TH1D *)f->Get("mc_pass_gen");
	TH1D *mc_samp	 = (TH1D *)f->Get("mc_samp");	

	auto tof_eff = GetEff(pass_tof,sample_tof);
	auto l1_eff  = GetEff(pass_l1, sample_l1);
	auto tr_eff  = GetEff(pass_tr, sample_tr);
	auto trig_eff= GetEff(phys_trig,sample_trig);
	auto true_acc= GetEff(mc_pass_gen, mc_samp,"simple");
	true_acc->Scale(TMath::Pi() * 3.9*3.9, "nosw2");
	
	tof_eff->SetLineWidth(2);
	l1_eff->SetLineWidth(2);
	tr_eff->SetLineWidth(2);
	trig_eff->SetLineWidth(2);
	true_acc->SetLineWidth(2);
	
	outfile->WriteTObject(tof_eff,"tof_eff");
	outfile->WriteTObject(l1_eff,"l1_eff");
	outfile->WriteTObject(tr_eff,"tr_eff");
	outfile->WriteTObject(trig_eff,"trig_eff");
	outfile->WriteTObject(true_acc,"true_acc");
}
