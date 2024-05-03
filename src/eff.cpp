#include "binning.h"
#include "utils.h"

TH1D* GetEff(TH1* hpass, TH1* htotal, TString option);

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: \n");
	    printf("%s <charge> <mc/dat/both> \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString opt = argv[2];
	TString ionPath=getIonPath(charge);
	TString outPath = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/"+ionPath+"_output.root";
	auto outfile = new TFile(outPath.Data() ,"recreate");
////DATA///////	
if (opt.Contains("dat") ) {
	TString datPath = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/dat/dat.root";
	TFile *dat = new TFile(datPath.Data() );
 	auto counts_2d = (TH2D *)dat->Get("rigidity");
 	int start_bin=1, stop_bin=counts_2d->GetNbinsX();
	TH1D *d_pass_tof   = (TH1D*)((TH2F *)dat->Get("pass_tof"))->ProjectionY("i",start_bin,stop_bin);
	TH1D *d_sample_tof = (TH1D*)((TH2F *)dat->Get("sample_tof"))->ProjectionY("j",start_bin,stop_bin);
	TH1D *d_pass_l1	   = (TH1D*)((TH2F *)dat->Get("pass_l1"))->ProjectionY("k",start_bin,stop_bin);
	TH1D *d_sample_l1  = (TH1D*)((TH2F *)dat->Get("sample_l1"))->ProjectionY("l",start_bin,stop_bin);
	TH1D *d_pass_tr    = (TH1D*)((TH2F *)dat->Get("pass_tr"))->ProjectionY("m",start_bin,stop_bin);
	TH1D *d_sample_tr  = (TH1D*)((TH2F *)dat->Get("sample_tr"))->ProjectionY("n",start_bin,stop_bin);
	TH1D *d_sample_trig= (TH1D*)((TH2F *)dat->Get("sample_trig"))->ProjectionY("o",start_bin,stop_bin);
	TH1D *d_phys_trig  = (TH1D*)((TH2F *)dat->Get("phys_trig"))->ProjectionY("p",start_bin,stop_bin);
	auto d_tof_eff = GetEff(d_pass_tof,d_sample_tof,"efficiency");
	auto d_l1_eff  = GetEff(d_pass_l1, d_sample_l1,"efficiency");
	auto d_tr_eff  = GetEff(d_pass_tr, d_sample_tr,"efficiency");
	auto d_trig_eff= GetEff(d_phys_trig,d_sample_trig,"efficiency");
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
	outfile->WriteTObject(d_tof_eff,"d_tof_eff");
	outfile->WriteTObject(d_l1_eff,"d_l1_eff");
	outfile->WriteTObject(d_tr_eff,"d_tr_eff");
	outfile->WriteTObject(d_trig_eff,"d_trig_eff");
	outfile->WriteTObject(lt,"lt");
	outfile->WriteTObject(counts,"counts");
	outfile->WriteTObject(oldflux,"oldflux");
}
////MONTECARLO////
if (opt.Contains("mc") ) {
	TString mcPath = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/mc/mc.root";
	TFile *mc    = new TFile(mcPath.Data() );
	auto onlyToGetBin  = (TH2F *)mc->Get("pass_tof");
	int start_bin=1, stop_bin=onlyToGetBin->GetNbinsX();
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
	outfile->WriteTObject(m_tof_eff,"m_tof_eff");
	outfile->WriteTObject(m_l1_eff,"m_l1_eff");
	outfile->WriteTObject(m_tr_eff,"m_tr_eff");
	outfile->WriteTObject(m_trig_eff,"m_trig_eff");
	outfile->WriteTObject(true_acc,"true_acc"); 
}
////BOTH WITH SPLINE////
if (opt.Contains("both")) {
	TString datPath = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/dat/dat.root";
	TFile *dat = new TFile(datPath.Data() );
	TString mcPath = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/mc/mc.root";
	TFile *mc    = new TFile(mcPath.Data() );
	auto counts_2d = (TH2D *)dat->Get("rigidity");
	int start_bin=1, stop_bin=counts_2d->GetNbinsX();
	TH1D *d_pass_tof   = (TH1D*)((TH2F *)dat->Get("pass_tof"))->ProjectionY("i",start_bin,stop_bin);
	TH1D *d_sample_tof = (TH1D*)((TH2F *)dat->Get("sample_tof"))->ProjectionY("j",start_bin,stop_bin);
	TH1D *d_pass_l1	   = (TH1D*)((TH2F *)dat->Get("pass_l1"))->ProjectionY("k",start_bin,stop_bin);
	TH1D *d_sample_l1  = (TH1D*)((TH2F *)dat->Get("sample_l1"))->ProjectionY("l",start_bin,stop_bin);
	TH1D *d_pass_tr    = (TH1D*)((TH2F *)dat->Get("pass_tr"))->ProjectionY("m",start_bin,stop_bin);
	TH1D *d_sample_tr  = (TH1D*)((TH2F *)dat->Get("sample_tr"))->ProjectionY("n",start_bin,stop_bin);
	TH1D *d_sample_trig= (TH1D*)((TH2F *)dat->Get("sample_trig"))->ProjectionY("o",start_bin,stop_bin);
	TH1D *d_phys_trig  = (TH1D*)((TH2F *)dat->Get("phys_trig"))->ProjectionY("p",start_bin,stop_bin);
	TH1D *m_pass_tof   = (TH1D*)((TH2F *)mc->Get("pass_tof"))->ProjectionY("a",start_bin,stop_bin);
	TH1D *m_sample_tof = (TH1D*)((TH2F *)mc->Get("sample_tof"))->ProjectionY("b",start_bin,stop_bin);
	TH1D *m_pass_l1	   = (TH1D*)((TH2F *)mc->Get("pass_l1"))->ProjectionY("c",start_bin,stop_bin);
	TH1D *m_sample_l1  = (TH1D*)((TH2F *)mc->Get("sample_l1"))->ProjectionY("d",start_bin,stop_bin);
	TH1D *m_pass_tr    = (TH1D*)((TH2F *)mc->Get("pass_tr"))->ProjectionY("e",start_bin,stop_bin);
	TH1D *m_sample_tr  = (TH1D*)((TH2F *)mc->Get("sample_tr"))->ProjectionY("f",start_bin,stop_bin);
	TH1D *m_sample_trig= (TH1D*)((TH2F *)mc->Get("sample_trig"))->ProjectionY("g",start_bin,stop_bin);
	TH1D *m_phys_trig  = (TH1D*)((TH2F *)mc->Get("phys_trig"))->ProjectionY("h",start_bin,stop_bin);
	TH1D *mc_pass_gen= (TH1D *)mc->Get("mc_pass_gen");
	TH1D *mc_samp	 = (TH1D *)mc->Get("mc_samp");
	TH1D *mc_pass    = (TH1D *)mc->Get("mc_pass");
	mc_samp->Rebin(9);
	mc_pass->Rebin(9);
	mc_pass_gen->Rebin(9);
	auto hwidth = (TH1D*)hist_rig->Clone();
	for(int i=1; i<=hist_rig->GetNbinsX(); ++i) {
		hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
		hwidth->SetBinError(i, 0);
	}
	auto lvt_25     = (TH1D *)((TH2D *)dat->Get("lvt_25"))->ProjectionY("q", start_bin, stop_bin);
	auto counts    = (TH1D *)((TH2D *)dat->Get("rigidity"))->ProjectionY("r", start_bin, stop_bin);
	auto oldflux   = (TH1D*)counts->Clone();
	auto true_acc = GetEff(mc_pass_gen, mc_samp, "simple");
	true_acc->Scale(TMath::Pi() * 3.9 * 3.9, "nosw2");
	auto acc = (TH1D*)true_acc->Clone();
	auto final_damc_tot = new TH1D(); // damc = data/montecarlo corr
	auto final_acc = new TH1D();
	auto flux = (TH1D *)counts->Clone();
	flux->Divide(lvt_25);
	flux->Divide(hwidth);
	auto dataEff_l1 = GetEff(d_pass_l1  , d_sample_l1,"simple");
	auto dataEff_tf = GetEff(d_pass_tof , d_sample_tof,"simple");
	auto dataEff_tr = GetEff(d_phys_trig, d_sample_trig,"simple");
	auto dataEff_tk = GetEff(d_pass_tr  , d_sample_tr,"simple");
	auto rate = (TH1D *)counts->Clone();
	rate->Divide(lvt_25);
	rate->Divide(hwidth);
	outfile->WriteTObject(rate,"rate");
	outfile->WriteTObject(lvt_25, "lvt_25");
	outfile->WriteTObject(counts, "counts");
	
	for(int iter=0; iter<8; ++iter) {
		    auto mcEff_l1 = GetEff(m_pass_l1  , m_sample_l1, "simple");
		    auto mcEff_tf = GetEff(m_pass_tof , m_sample_tof, "simple");
		    auto mcEff_tr = GetEff(m_phys_trig, m_sample_trig, "simple");
		    auto mcEff_tk = GetEff(m_pass_tr  , m_sample_tr, "simple");
		    auto damc_l1 = GetEff(dataEff_l1, mcEff_l1, "correction");
		    auto damc_tf = GetEff(dataEff_tf, mcEff_tf, "correction");
		    auto damc_tr = GetEff(dataEff_tr, mcEff_tr, "correction");
		    auto damc_tk = GetEff(dataEff_tk, mcEff_tk, "correction");
		    auto final_damc_l1 = autospline(damc_l1, 0.5, 30);
		    auto final_damc_tf = autospline(damc_tf, 0.5, 15);
		    auto final_damc_tr = autospline(damc_tr, 0.5, 30);
		    auto final_damc_tk = autospline(damc_tk, 0.5, 7);
		    final_damc_tot = (TH1D*)final_damc_l1->Clone();
		    final_damc_tot->Multiply(final_damc_tf);
		    final_damc_tot->Multiply(final_damc_tr);
		    final_damc_tot->Multiply(final_damc_tk);
		    if(iter) {
		      acc = GetEff(mc_pass, mc_samp, "simple");
		      acc->Scale(TMath::Pi() * 3.9 * 3.9, "nosw2");
		    }
		    auto spline_acc = autospline(acc, 0.65, 45);
		    final_acc = (TH1D*)spline_acc->Clone();
		    final_acc->Multiply(final_damc_tot);
		    auto oldflux = (TH1D*)flux->Clone();
		    flux = (TH1D*)counts->Clone();
		    flux->Divide(lvt_25);
		    flux->Divide(hwidth);
		    flux->Divide(final_acc);
		    if (iter ==7 ) {
				outfile->WriteTObject(flux, Form("flux_%i", iter));
				outfile->WriteTObject(damc_l1, Form("damc_l1_%i", iter));
				outfile->WriteTObject(damc_tf, Form("damc_tf_%i", iter));
				outfile->WriteTObject(damc_tr, Form("damc_tr_%i", iter));
				outfile->WriteTObject(damc_tk, Form("damc_tk_%i", iter));
				outfile->WriteTObject(final_damc_l1, Form("final_damc_l1_%i", iter));
				outfile->WriteTObject(final_damc_tf, Form("final_damc_tf_%i", iter));
				outfile->WriteTObject(final_damc_tr, Form("final_damc_tr_%i", iter));
				outfile->WriteTObject(final_damc_tk, Form("final_damc_tk_%i", iter));
				outfile->WriteTObject(final_damc_tot, Form("final_damc_tot_%i", iter));
				outfile->WriteTObject(acc, Form("acc_%i", iter));
				outfile->WriteTObject(spline_acc, Form("spline_acc_%i", iter));
				outfile->WriteTObject(final_acc, Form("final_acc_%i", iter));
				outfile->WriteTObject(dataEff_l1,"dataEff_l1");
				outfile->WriteTObject(dataEff_tf,"dataEff_tf");
				outfile->WriteTObject(dataEff_tr,"dataEff_tr");
				outfile->WriteTObject(dataEff_tk,"dataEff_tk");
				outfile->WriteTObject(mcEff_l1,"mcEff_l1");
				outfile->WriteTObject(mcEff_tf,"mcEff_tf");
				outfile->WriteTObject(mcEff_tr,"mcEff_tr");
				outfile->WriteTObject(mcEff_tk,"mcEff_tk");
			}
	}
}
////WRITE////
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
