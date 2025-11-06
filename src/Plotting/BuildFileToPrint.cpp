#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"
#include "SplineUtility.h"

void printUsage(const char* exeName) {
    std::cout << "Usage: \n"
              << exeName << " <time> <charge> <all/single> <rebin/no> [nRebin] \n"
              << "  <all/single>:\n"
              << "    all    - efficiencies from general file (output of selectEv.cpp)\n"
              << "    single - efficiencies from single files (output of XEff.cpp)\n"
              << "  <rebin/no>:\n"
              << "    rebin  - rebin all efficiencies\n"
              << "    no     - no rebinning\n"
              << "    [nRebin]: (optional), required if rebin is selected. Suggested: 4\n";
}

int main(int argc, char *argv[]) {
  if (argc < 5 || (TString(argv[3]) == "rebin" && argc < 5)) {
    printUsage(argv[0]);
    return 1;
  }
  TString timePeriod = argv[1];
  TString exename = argv[0];
	unsigned int charge=atoi(argv[2]);
	TString opt2 = argv[3];
  TString rebin =argv[4];
  int nRebin = (rebin == "rebin") ? atoi(argv[4]) : 1;
  auto outfile = new TFile("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+
                             getIonPath(charge)+"/"+getIonPath(charge)+"_output.root","recreate");
//--Build the efficiencies
  auto r = StreamUtility::loadFile(charge,opt2,"","");
	TFile *dal1 = r.dal1;
	TFile *datf = r.datf;
	TFile *datr = r.datr;
	TFile *datk = r.datk;
	TFile *mcl1 = r.mcl1;
	TFile *mctf = r.mctf;
	TFile *mctr = r.mctr;
	TFile *mctk = r.mctk;
  TFile *livetime = r.livetime;
  TFile *rawacc   = r.rawacc;
	TFile *Counts   = r.counts;
	auto counts_2d = (TH2D *)Counts->Get("rigidity");
	int start_bin=1, stop_bin=counts_2d->GetNbinsX();
	TH1D *d_pass_tof   = (TH1D*)((TH2F *)datf->Get("pass_tof_den_ltotf_05"))->ProjectionY("i",start_bin,stop_bin);
	TH1D *d_sample_tof = (TH1D*)((TH2F *)datf->Get("sample_tof_den_ltof_05"))->ProjectionY("j",start_bin,stop_bin);
	TH1D *d_pass_l1	   = (TH1D*)((TH2F *)dal1->Get("pass_l1"))->ProjectionY("k",start_bin,stop_bin);
	TH1D *d_sample_l1  = (TH1D*)((TH2F *)dal1->Get("sample_l1"))->ProjectionY("l",start_bin,stop_bin);
	TH1D *d_pass_tr    = (TH1D*)((TH2F *)datk->Get("pass_tr"))->ProjectionY("m",start_bin,stop_bin);
	TH1D *d_sample_tr  = (TH1D*)((TH2F *)datk->Get("sample_tr"))->ProjectionY("n",start_bin,stop_bin);
	TH1D *d_sample_trig= (TH1D*)((TH2F *)datr->Get("unbiased_trig"))->ProjectionY("o",start_bin,stop_bin);
	TH1D *d_phys_trig  = (TH1D*)((TH2F *)datr->Get("phys_trig"))->ProjectionY("p",start_bin,stop_bin);
	TH1D *m_pass_tof   = (TH1D*)((TH2F *)mctf->Get("pass_tof"))->ProjectionY("a",start_bin,stop_bin);
	TH1D *m_sample_tof = (TH1D*)((TH2F *)mctf->Get("sample_tof"))->ProjectionY("b",start_bin,stop_bin);
	TH1D *m_pass_l1	   = (TH1D*)((TH2F *)mcl1->Get("pass_l1"))->ProjectionY("c",start_bin,stop_bin);
	TH1D *m_sample_l1  = (TH1D*)((TH2F *)mcl1->Get("sample_l1"))->ProjectionY("d",start_bin,stop_bin);
	TH1D *m_pass_tr    = (TH1D*)((TH2F *)mctk->Get("pass_tr"))->ProjectionY("e",start_bin,stop_bin);
	TH1D *m_sample_tr  = (TH1D*)((TH2F *)mctk->Get("sample_tr"))->ProjectionY("f",start_bin,stop_bin);
	TH1D *m_sample_trig= (TH1D*)((TH2F *)mctr->Get("unbiased_trig"))->ProjectionY("g",start_bin,stop_bin);
	TH1D *m_phys_trig  = (TH1D*)((TH2F *)mctr->Get("phys_trig"))->ProjectionY("h",start_bin,stop_bin);
	TH1D *mc_pass_gen  = (TH1D *)rawacc->Get("mc_pass_gen");
	TH1D *mc_samp	   = (TH1D *)rawacc->Get("mc_samp");
	TH1D *mc_pass      = (TH1D *)rawacc->Get("mc_pass");
  mc_samp->Rebin(9);
  mc_pass->Rebin(9);
  mc_pass_gen->Rebin(9);
  auto hwidth = (TH1D*)hist_rig_highZ->Clone();
	for(int i=1; i<=hist_rig_highZ->GetNbinsX(); ++i) {
		hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
		hwidth->SetBinError(i, 0);
	}
	auto lvt_25    = (TH1D *)((TH2D *)livetime->Get("lvt_25"))->ProjectionY("q", start_bin, stop_bin);
	auto counts    = (TH1D *)((TH2D *)Counts->Get("rigidity"))->ProjectionY("r", start_bin, stop_bin);
  auto rate = (TH1D *)counts->Clone();
	rate->Divide(lvt_25);
	rate->Divide(hwidth);
  auto rateMulti = (TH1D*)rate->Clone();
  MultiplyByXPower(rateMulti,2.7);
//-------Rebin efficiencies
  int tRebin = 4;
  if (rebin=="rebin") {
    printf("--- Rebinning tof efficiency.. \n");
    Rebin(d_sample_tof,nRebin);
    Rebin(d_pass_tof,  nRebin);
    Rebin(m_pass_tof,  nRebin);
    Rebin(m_sample_tof,nRebin);
    printf("--- Rebinning l1 efficiency.. \n");
    Rebin(d_sample_l1, nRebin);
    Rebin(d_pass_l1,   nRebin);
    Rebin(m_pass_l1,   nRebin);
    Rebin(m_sample_l1, nRebin);
    printf("--- Rebinning track efficiency.. \n");
    Rebin(d_sample_tr, nRebin);
    Rebin(d_pass_tr,   nRebin);
    Rebin(m_pass_tr,   nRebin);
    Rebin(m_sample_tr, nRebin);
    printf("--- Rebinning trigger efficiency.. \n");
    Rebin(d_sample_trig, nRebin);
    Rebin(d_phys_trig,   nRebin);
    Rebin(m_sample_trig, nRebin);
    Rebin(m_phys_trig,   nRebin);
  }
  auto l1_da = divide(d_pass_l1, d_sample_l1, "efficiency");
  auto tf_da = divide(d_pass_tof, d_sample_tof, "efficiency");
  auto tr_da = divide(d_phys_trig, d_sample_trig, "efficiency");
  auto tk_da = divide(d_pass_tr, d_sample_tr, "efficiency");
  auto l1_mc = divide(m_pass_l1, m_sample_l1, "efficiency");
  auto tf_mc = divide(m_pass_tof, m_sample_tof, "efficiency");
  auto tr_mc = divide(m_phys_trig, m_sample_trig, "efficiency");
  auto tk_mc = divide(m_pass_tr, m_sample_tr, "efficiency");
//--Build the ratio
  auto damc_l1 = divide(l1_da, l1_mc, "ratio");
  auto damc_tf = divide(tf_da, tf_mc, "ratio");
  auto damc_tr = divide(tr_da, tr_mc, "ratio");
  auto damc_tk = divide(tk_da, tk_mc, "ratio");
  auto l1_range = SetFitLimits(SplineUtility::Efficiency::L1Eff,charge);
  auto tf_range = SetFitLimits(SplineUtility::Efficiency::TofEff,charge);
  auto tr_range = SetFitLimits(SplineUtility::Efficiency::TriggerEff,charge);
  auto tk_range = SetFitLimits(SplineUtility::Efficiency::TrackEff,charge);
//--Build the spline
  auto final_damc_l1 = autospline(damc_l1, l1_range.first, l1_range.second);
  auto final_damc_tf = autospline(damc_tf, tf_range.first, tf_range.second);
  auto final_damc_tr = autospline(damc_tr, tr_range.first, tr_range.second,3,7);
  auto final_damc_tk = autospline(damc_tk, tk_range.first, tk_range.second);
//--Final correction
  auto final_damc_tot = (TH1D*)final_damc_l1->Clone();
  final_damc_tot->Multiply(final_damc_tf);
  final_damc_tot->Multiply(final_damc_tr);
  final_damc_tot->Multiply(final_damc_tk);
  auto acc = divide(mc_pass_gen, mc_samp, "efficiency");
  acc->Scale(TMath::Pi() * 3.9 * 3.9);
  auto spline_acc = autospline(acc, 2.15, 2000,5,20);
	auto final_acc = (TH1D*)spline_acc->Clone();
	final_acc->Multiply(final_damc_tot);
//--Folded lux    
  auto flux = GetFlux(charge,timePeriod);
//--Flux model for unfolding and its fit
  TString fileName = "../IonsSelected/"+getIonPath(charge)+"/Unfolding/"+getIonPath(charge)+"_UnfoldingFactor.root";
  TFile *file = TFile::Open(fileName.Data());
  TH1D *pub_flux,*unf,*sp_unf;
  TF1 *flux_model;
  if (file) {
    printf("--- Opening file: %s \n",fileName.Data() );
    pub_flux = (TH1D*)file->Get("pub_flux");
    flux_model=(TF1*)file->Get("flux_model_0");
  //--Unfolding factor
    unf = (TH1D*)file->Get("unf_factor");
    sp_unf = autospline(unf,2.7,1000);
  }
  else {
    printf("--- Unable to open file: %s \n",fileName.Data() );
  }
    outfile->WriteTObject(counts,"counts");
    outfile->WriteTObject(lvt_25,"lvt_25");
    outfile->WriteTObject(rate,"rate");
    outfile->WriteTObject(rateMulti,"rateMulti");
    outfile->WriteTObject(l1_da,"l1_da");
    outfile->WriteTObject(tf_da,"tf_da");
    outfile->WriteTObject(tr_da,"tr_da");
    outfile->WriteTObject(tk_da,"tk_da");
    outfile->WriteTObject(l1_mc,"l1_mc");
    outfile->WriteTObject(tf_mc,"tf_mc");
    outfile->WriteTObject(tr_mc,"tr_mc");
    outfile->WriteTObject(tk_mc,"tk_mc");
    outfile->WriteTObject(damc_l1, "damc_l1");
    outfile->WriteTObject(damc_tf, "damc_tf" );
    outfile->WriteTObject(damc_tr, "damc_tr" );
    outfile->WriteTObject(damc_tk, "damc_tk" );
    outfile->WriteTObject(final_damc_l1, "final_damc_l1");
    outfile->WriteTObject(final_damc_tf, "final_damc_tf");
    outfile->WriteTObject(final_damc_tr, "final_damc_tr");
    outfile->WriteTObject(final_damc_tk, "final_damc_tk");
    outfile->WriteTObject(final_damc_tot, "final_damc_tot");
    outfile->WriteTObject(acc,"acc");
    outfile->WriteTObject(spline_acc,"spline_acc");
    outfile->WriteTObject(final_acc,"final_acc");
    outfile->WriteTObject(pub_flux,"pub_flux");  
    outfile->WriteTObject(flux_model,"flux_model_0");   
    outfile->WriteTObject(unf,"unf_factor");
    outfile->WriteTObject(sp_unf,"sp_unf");
    outfile->WriteTObject(flux,"flux");
    
    outfile->Close();
}