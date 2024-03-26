 #include "binning.h"

TH1D* GetEff(TH1* hpass, TH1* htotal, TString option);


int main() {
	TFile *dat = new TFile("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/dat.root");
	TFile *mc    = new TFile("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/mc.root");
	auto outfile = new TFile("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/output.root","recreate");
////MONTECARLO////
	TH1D *m_pass_tof   = (TH1D*)((TH2F *)mc->Get("pass_tof"))->ProjectionY("a",start_bin,stop_bin);
	TH1D *m_sample_tof = (TH1D*)((TH2F *)mc->Get("sample_tof"))->ProjectionY("b",start_bin,stop_bin);
	TH1D *m_pass_l1	   = (TH1D*)((TH2F *)mc->Get("pass_l1"))->ProjectionY("c",start_bin,stop_bin);
	TH1D *m_sample_l1  = (TH1D*)((TH2F *)mc->Get("sample_l1"))->ProjectionY("d",start_bin,stop_bin);
	TH1D *m_pass_tr    = (TH1D*)((TH2F *)mc->Get("pass_tr"))->ProjectionY("e",start_bin,stop_bin);
	TH1D *m_sample_tr  = (TH1D*)((TH2F *)mc->Get("sample_tr"))->ProjectionY("f",start_bin,stop_bin);
	TH1D *m_sample_trig= (TH1D*)((TH2F *)mc->Get("sample_trig"))->ProjectionY("g",start_bin,stop_bin);
	TH1D *m_phys_trig  = (TH1D*)((TH2F *)mc->Get("phys_trig"))->ProjectionY("h",start_bin,stop_bin);
	auto m_tof_eff = GetEff(m_pass_tof,m_sample_tof,"simple");
	auto m_l1_eff  = GetEff(m_pass_l1, m_sample_l1,"simple");
	auto m_tr_eff  = GetEff(m_pass_tr, m_sample_tr,"simple");
	auto m_trig_eff= GetEff(m_phys_trig,m_sample_trig,"simple");
	TH1D *mc_pass_gen= (TH1D *)mc->Get("mc_pass_gen");
	TH1D *mc_samp	 = (TH1D *)mc->Get("mc_samp");
	TH1D *mc_pass    = (TH1D *)mc->Get("mc_pass");
	mc_samp->Rebin(9);
	mc_pass->Rebin(9);
	mc_pass_gen->Rebin(9);
	auto true_acc= GetEff(mc_pass_gen, mc_samp,"simple");
	true_acc->Scale(TMath::Pi() * 3.9*3.9, "nosw2");
	true_acc->SetName("true_acc");
////DATA///////	
	TH1D *d_pass_tof   = (TH1D*)((TH2F *)dat->Get("pass_tof"))->ProjectionY("i",start_bin,stop_bin);
	TH1D *d_sample_tof = (TH1D*)((TH2F *)dat->Get("sample_tof"))->ProjectionY("j",start_bin,stop_bin);
	TH1D *d_pass_l1	   = (TH1D*)((TH2F *)dat->Get("pass_l1"))->ProjectionY("k",start_bin,stop_bin);
	TH1D *d_sample_l1  = (TH1D*)((TH2F *)dat->Get("sample_l1"))->ProjectionY("l",start_bin,stop_bin);
	TH1D *d_pass_tr    = (TH1D*)((TH2F *)dat->Get("pass_tr"))->ProjectionY("m",start_bin,stop_bin);
	TH1D *d_sample_tr  = (TH1D*)((TH2F *)dat->Get("sample_tr"))->ProjectionY("n",start_bin,stop_bin);
	TH1D *d_sample_trig= (TH1D*)((TH2F *)dat->Get("sample_trig"))->ProjectionY("o",start_bin,stop_bin);
	TH1D *d_phys_trig  = (TH1D*)((TH2F *)dat->Get("phys_trig"))->ProjectionY("p",start_bin,stop_bin);
	auto d_tof_eff = GetEff(d_pass_tof,d_sample_tof,"simple");
	auto d_l1_eff  = GetEff(d_pass_l1, d_sample_l1,"simple");
	auto d_tr_eff  = GetEff(d_pass_tr, d_sample_tr,"simple");
	auto d_trig_eff= GetEff(d_phys_trig,d_sample_trig,"simple");
	auto hwidth = (TH1D*)hist_rig->Clone();
	for(int i=1; i<=hist_rig->GetNbinsX(); ++i) {
		hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
		hwidth->SetBinError(i, 0);
	}
	auto lt     = (TH1D *)((TH2D *)dat->Get("lvt_25"))->ProjectionY("j", start_bin, stop_bin);
	auto counts = (TH1D *)((TH2D *)dat->Get("rigidity"))->ProjectionY("k", start_bin, stop_bin);
	auto oldflux = (TH1D*)counts->Clone();
	oldflux->Divide(lt);
	oldflux->Divide(hwidth);
////CORRECTIONS///
	auto damc_l1     = GetEff(d_l1_eff,   m_l1_eff,   "correction");
    auto damc_tof    = GetEff(d_tof_eff,  m_tof_eff,  "correction");
    auto damc_track  = GetEff(d_tr_eff,   m_tr_eff,   "correction");
    auto damc_trig   = GetEff(d_trig_eff, m_trig_eff, "correction");
	auto final_damc_tot = (TH1D*)damc_l1->Clone();
    final_damc_tot->Multiply(damc_tof);
    final_damc_tot->Multiply(damc_track);
    final_damc_tot->Multiply(damc_trig);
	auto final_acc = (TH1D*)true_acc->Clone();
	final_acc ->Multiply(final_damc_tot);
    auto flux = (TH1D*)counts->Clone();
    flux->Divide(lt);
    flux->Divide(hwidth);
    flux->Divide(final_acc);
////WRITE////
	outfile->WriteTObject(d_tof_eff,"d_tof_eff");
	outfile->WriteTObject(d_l1_eff,"d_l1_eff");
	outfile->WriteTObject(d_tr_eff,"d_tr_eff");
	outfile->WriteTObject(d_trig_eff,"d_trig_eff");
	outfile->WriteTObject(m_tof_eff,"m_tof_eff");
	outfile->WriteTObject(m_l1_eff,"m_l1_eff");
	outfile->WriteTObject(m_tr_eff,"m_tr_eff");
	outfile->WriteTObject(m_trig_eff,"m_trig_eff");
	outfile->WriteTObject(true_acc,"true_acc");
	outfile->WriteTObject(lt,"lt");
	outfile->WriteTObject(counts,"counts");
	outfile->WriteTObject(flux,"flux");
	outfile->WriteTObject(oldflux,"oldflux");
	outfile->Close();
}


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
