#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"
#include "SplineUtility.h"
#include <TParameter.h>

TH1D *BuildPublishedFlux(unsigned int charge);
TGraph *build_histo_model(TH1D *histo);
double integrate(const TGraph* fun, double min, double max);
double integrate(const TF1* fun, double min, double max);
void RigidityMap(unsigned int charge, TH1D *h);
int getEntriesBySample(int sampleNum, TTree *t1, TTree *t2, TTree *t3, TTree *t4);
TString SamplePrefix(int sampleNum);
double maxDaqFit(unsigned int charge);
TH1D *buildDummyFlux(unsigned int charge, TString timePeriod, TString inp_sec_track);
int Z;
double A,mass;
const double prod_min = 0.9957, prod_max = 2000;
double RigToEkn(double Rig, double Z, double A, double M) {
  return (sqrt((Rig * Z) * (Rig * Z) + M * M) - M) / A;
}
double EknToRig(double Ekn, double Z, double A, double M) {
  // return sqrt((Ekn * A + M) * (Ekn * A + M) - (M * M)) / Z;
  return sqrt(A * Ekn * (A * Ekn + 2. * M)) / Z;
}
//original
double force_field(double* x, double* p) {
  double Rig      = x[0];
  double Norm     = p[0];
  double phi      = p[1];
  double exponent = p[2];

  double Etot = RigToEkn(Rig, Z, A, mass) + mass;

  double SM_factor = (Etot * Etot - mass * mass) / ((Etot + phi) * (Etot + phi) - mass * mass);

  double flux = pow(EknToRig(Etot - mass + phi, Z, A, mass), exponent);

  return Norm * flux * SM_factor;
}
//original and multiplied
double force_field_multiplied(double* x, double* p) {
  double Rig      = x[0];
  double Norm     = p[0];
  double phi      = p[1];
  double exponent = p[2];

  double Etot = RigToEkn(Rig, Z, A, mass) + mass;

  double SM_factor = (Etot * Etot - mass * mass) / ((Etot + phi) * (Etot + phi) - mass * mass);

  double flux = pow(EknToRig(Etot - mass + phi, Z, A, mass), exponent);

  return Norm * flux * SM_factor * pow(Rig,2.7);
}
//with broken power law
double broken_force_field(double* x, double* p) {
  double Rig    = x[0];
  double Norm   = p[0];  // Normalization
  double phi    = p[1];  // Modulation potential (GV)
  double gamma1 = p[2];  // Index below break
  double gamma2 = p[3];  // Index above break
  double rBreak = p[4];  // Break rigidity

  if (Rig <= 0) return 0;

  // Total energy
  double Etot = RigToEkn(Rig, Z, A, mass) + mass;

  // Solar modulation factor
  double SM_factor = (Etot * Etot - mass * mass) /
                     ((Etot + phi) * (Etot + phi) - mass * mass);

  // Broken power law with continuity at break
  double flux;
  if (Rig < rBreak) {
    flux = pow(Rig, gamma1);
  } else {
    double continuity = pow(rBreak, gamma1 - gamma2);
    flux = continuity * pow(Rig, gamma2);
  }

  return Norm * flux * SM_factor;
}
TF1* build_flux_model(unsigned int Z,
                      TH1* flux_hist,
                      double prod_min,
                      double prod_max) {
    // aggiorno le globali usate in force_field
    ::Z    = Z;
    ::A    = 2*Z;
    ::mass = 0.938272075 * ::A;

    auto f = new TF1(Form("flux_model_Z%d", Z),
                     force_field,
                     prod_min,
                     prod_max,
                     3);
    f->SetParameters(20, 0.6, -2.7);
    flux_hist->Fit(f);
    f->SetNpx(100000);

    double chi2    = f->GetChisquare();
    int npoints    = f->GetNumberFitPoints();
    int npar       = f->GetNumberFreeParameters();
    int ndf        = npoints - npar;
    printf("[Z=%d] Chi2/ndf = %.2f / %d = %.2f\n",
           Z, chi2, ndf, chi2/ndf);

    return f;
} 
void printUsage(const char* exeName) {
    std::cout << "\nUsage:\n"
              << "  " << exeName << " <charge> <all/single> <reweigh/no> <rebin/no> <rigmap/no> <mc> <normal/emulate> <period> -sec_track <y/n> <removeContamination> <geom> <treeType>\n"
              << "Arguments:\n"
              << "  <all/single>:\n"
              << "    all          Efficiencies from general file (output of selectEv.cpp)\n"
              << "    single       Efficiencies from single files (output of XEff.cpp)\n\n"
              << "  <reweigh/no>:\n"
              << "    r            Reweigh the Monte Carlo\n"
              << "    no           No reweigh\n\n"
              << "  <rebin/no>:\n"
              << "    rebin        Rebin track data and Monte Carlo efficiency\n"
              << "    no           No rebinning\n\n"
              << "  <rigmap/no>:\n"
              << "    rigmap       Reweigh track efficiency using rigidity map\n"
              << "    no           No reweigh\n\n"
              << "  <mc>:\n"
              << "    s            Use single MC\n"
              << "    g            Use global MC\n\n"
              << "  <normal/emulate>:\n"
              << "    normal  (-n) Trigger efficiency with unbiased method\n"
              << "    emulate (-e) Trigger efficiency with emulation\n\n"
              << "  <period>:\n"
              << "    10y          Use 10-year period\n"
              << "    12.5y        Use 12.5-year period\n"
              << "    13.5y        Use 13.5-year period\n\n"
              << "  -sec_track:\n"
              << "    y            Enable secondary track selection\n"
              << "    n            Disable secondary track selection\n\n"
              << "  <removeContamination>:\n"
              << "    y            Remove contamination\n"
              << "    n            Keep contamination\n\n"
              << "  <geom>\n"
              << "    IL1          Inner L1\n"
              << "    FS           Full Span\n"
              << "  <treeType>\n"
              << "    old          reco_gen1, ..\n";
}


int main(int argc, char *argv[]) {
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TH1::SetDefaultSumw2(true);
  if (argc < 13) {
    printUsage(argv[0]);
    return 1;
  }
  TString exename = argv[0];
  unsigned int charge = atoi(argv[1]);
  TString opt2 = argv[2];
  TString rew = argv[3];
  TString rebin = argv[4];
  TString rigmap = argv[5];
  TString mcType = argv[6];
  TString trgType = argv[7];
  TString timePeriod = argv[8];
  TString sec_track = argv[9];
  TString optBkg = argv[10];
  TString geom = argv[11];
  TString treeType = argv[12];

  TString ionPath = getIonPath(charge);
  TString outPath;
  TString outname;
  TString sec_track_out;
  if (sec_track=="y") sec_track_out = "_sec_track";
  else sec_track_out = "";

  //outname = StreamUtility::formatOutnameOptions(charge,getIonPath(charge)+"_UnfoldingFactor",sec_track,"",mcType, timePeriod);
  if (mcType=="s")
    outPath = StreamUtility::getOutputDir(charge, exename, ionPath + 
                                          Form("_UnfoldingFactor%s_%s.root",sec_track_out.Data(),timePeriod.Data()));
  if (mcType=="g")
    outPath = StreamUtility::getOutputDir(charge, exename, ionPath + 
                                          Form("_UnfoldingFactor_global%s_%s.root",sec_track_out.Data(),timePeriod.Data()));
  std::cout << outPath << std::endl;
  auto outfile = new TFile(outPath.Data(), "recreate");
//-------Define the corrections and the efficiencies
  TH1D  *damc_l1,*damc_tf,*damc_tr,*damc_tk,*damc_l1u,*damc_tkch,
        *final_damc_l1,*final_damc_tf,*final_damc_tr,*final_damc_tk,*final_damc_l1u,*final_damc_tkch,*final_daq,*final_l1ch,
        *final_damc_tot,
        *l1_mc,*tf_mc,*tr_mc,*tk_mc,*l1u_mc,*tkch_mc,
        *tr_da, *daq_da, *l1ch,
        *acc,*spline_acc;
//-------Load the histograms------------  
  auto r = StreamUtility::loadFile(charge,opt2,mcType,sec_track);
//------Get the selected time period
	int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
  auto hwidth = (TH1D*)hist_rig_highZ->Clone();
	for(int i=1; i<=hist_rig_highZ->GetNbinsX(); ++i) {
		hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
		hwidth->SetBinError(i, 0);
	}
  auto hwidth_mc = (TH1D*)hist_log->Clone();
	for(int i=1; i<=hist_log->GetNbinsX(); ++i) {
		hwidth_mc->SetBinContent(i, hwidth_mc->GetBinWidth(i));
		hwidth_mc->SetBinError(i, 0);
	}
  outfile->WriteTObject(hwidth, "hwidth");
  outfile->WriteTObject(hwidth_mc, "hwidth_mc" );
  auto lvt_25 = (TH1D *)((TH2D *)r.livetime->Get("lvt_25;1"))->ProjectionY("q", start_bin, stop_bin);
  auto counts = (TH1D*) GetHist2D(r.counts, geom, "rigidity")->ProjectionY("r", start_bin, stop_bin);
  outfile->WriteTObject(counts, "countsRaw");
//-------Removing contamination below L1 and above L1
 if (optBkg=="y") {
  TString contName= "../Fragmentation/BelowL1/Fractions/contamination_"+getIonName(charge)+".root";
  TFile *contFile = TFile::Open(contName.Data());
  TString above = "../IonsSelected/"+getIonPath(charge)+"/Fragments/Result/fragment.root";
  TFile *aboveFile = TFile::Open(above.Data());
  auto fragm  = (TH1D*)aboveFile->Get("fraction");
  if (!contFile || contFile->IsZombie()) {
    printf("Errore nell'aprire il file per la purity\n");
  } if (!aboveFile || aboveFile->IsZombie()) {
    printf("Errore nell'aprire aboveL1\n");
  } else {
    auto purity = (TF1*)contFile->Get("purityFit");
    auto fragm  = (TH1D*)aboveFile->Get("fraction");
    //Removing contamination BelowL1 and AboveL1 from counts
    for (int ibin = 0; ibin < counts->GetNbinsX(); ibin++) {
      double lowEdge  = counts->GetBinLowEdge(ibin+1);
      double highEdge = counts->GetBinLowEdge(ibin+1) + counts->GetBinWidth(ibin+1);
      counts->SetBinContent(ibin + 1, counts->GetBinContent(ibin + 1)
                             *(1-purity->Eval(counts->GetBinCenter(ibin+1))
                                            - fragm->GetBinContent((ibin+1) )) ) ;
            double N = counts->GetBinContent(ibin + 1);
            double P = purity->Eval(counts->GetBinCenter(ibin+1));
            double F = fragm->GetBinContent(ibin + 1);
            double sigma_N = counts->GetBinError(ibin + 1);
            //double sigma_P = purity->GetBinError(ibin + 1);
            double sigma_F = fragm->GetBinError(ibin + 1);
            double sigma_hcounts = sqrt(pow((1 - P) * (1 - F) * sigma_N, 2) +
                                        //pow(N * (1 - F) * sigma_P, 2) +
                                        pow(N * (1 - P) * sigma_F, 2));
            counts->SetBinError(ibin + 1, 0.);
      }
  }
 }
  auto rateMulti = (TH1D *)counts->Clone();
	rateMulti->Divide(lvt_25);
	rateMulti->Divide(hwidth);
  MultiplyByXPower(rateMulti,2.7);
  outfile->WriteTObject(rateMulti,"rateMulti");
  outfile->WriteTObject(counts, "counts" );
  outfile->WriteTObject(lvt_25, "lvt_25" );
  auto d_pass_tof   = (TH1D*) GetHist2D(r.datf,  geom, "pass_tof_den_ltof")
                          ->ProjectionY("i", start_bin, stop_bin);
  auto d_sample_tof = (TH1D*) GetHist2D(r.datf,  geom, "sample_tof_den_ltof")
                          ->ProjectionY("j", start_bin, stop_bin);
  auto d_pass_l1    = (TH1D*) GetHist2D(r.dal1,  geom, "pass_l1_den")
                          ->ProjectionY("k", start_bin, stop_bin);
  auto d_sample_l1  = (TH1D*) GetHist2D(r.dal1,  geom, "sample_l1_den")
                          ->ProjectionY("l", start_bin, stop_bin);
  auto d_pass_tr    = (TH1D*) GetHist2D(r.datk,  geom, "pass_tr")
                          ->ProjectionY("m", start_bin, stop_bin);
  auto d_sample_tr  = (TH1D*) GetHist2D(r.datk,  geom, "sample_tr")
                          ->ProjectionY("n", start_bin, stop_bin);
  auto d_unb_trig   = (TH1D*) GetHist2D(r.datr,  geom, "unbiased_trig")
                          ->ProjectionY("o", start_bin, stop_bin);
  auto d_phys_trig  = (TH1D*) GetHist2D(r.datr,  geom, "phys_trig")
                          ->ProjectionY("p", start_bin, stop_bin);
  auto d_pass_l1u   = (TH1D*) GetHist2D(r.dal1u, geom, "pass_l1Unb_den")
                          ->ProjectionY("s", start_bin, stop_bin);
  auto d_sample_l1u = (TH1D*) GetHist2D(r.dal1u, geom, "sample_l1Unb_den")
                          ->ProjectionY("t", start_bin, stop_bin);
  auto d_pass_trch  = (TH1D*) GetHist2D(r.datkch,geom, "pass_trch_den_ltof_l9_1")
                          ->ProjectionY("u", start_bin, stop_bin);
  auto d_sample_trch= (TH1D*) GetHist2D(r.datkch,geom, "sample_trch_den_ltof_l9_1")
                          ->ProjectionY("v", start_bin, stop_bin);
  l1ch = getL1chEff(charge);
//-------Reweigh data track using rigidity map
if (rigmap=="rigmap") {
  printf("[INFO] Reweighing track data with rigidity map.. \n");
  RigidityMap(charge,d_pass_tr);
  RigidityMap(charge,d_sample_tr);
}
//-------Rebin efficiencies
    printf("[INFO] Rebinning trigger data efficiency.. \n");
    RebinByEfficiencyOverall(d_unb_trig,SplineUtility::Efficiency::TriggerEff,charge);
    RebinByEfficiencyOverall(d_phys_trig,SplineUtility::Efficiency::TriggerEff,charge);
  if (rebin=="rebin") {
    printf("[INFO] Rebinning l1 detect data efficiency.. \n");
    RebinByEfficiencyOverall(d_sample_l1u,SplineUtility::Efficiency::L1UnbEff,charge);
    RebinByEfficiencyOverall(d_pass_l1u,SplineUtility::Efficiency::L1UnbEff,charge);
    printf("[INFO] Rebinning track charge data efficiency.. \n");
    RebinByEfficiencyOverall(d_sample_trch,SplineUtility::Efficiency::TrackChEff,charge);
    RebinByEfficiencyOverall(d_pass_trch,SplineUtility::Efficiency::TrackChEff,charge);
    printf("[INFO] Rebinning track data efficiency.. \n");
    RebinByEfficiencyOverall(d_sample_tr,SplineUtility::Efficiency::TrackEff,charge);
    RebinByEfficiencyOverall(d_pass_tr,SplineUtility::Efficiency::TrackEff,charge);
    printf("[INFO] Rebinning tof data efficiency.. \n");
    RebinByEfficiencyOverall(d_sample_tof,SplineUtility::Efficiency::TofEff,charge);
    RebinByEfficiencyOverall(d_pass_tof,SplineUtility::Efficiency::TofEff,charge);
    printf("[INFO] Rebinning l1 data efficiency.. \n");
    RebinByEfficiencyOverall(d_sample_l1,SplineUtility::Efficiency::L1Eff,charge);
    RebinByEfficiencyOverall(d_pass_l1,SplineUtility::Efficiency::L1Eff,charge);
    /*RebinHistAbove(d_sample_l1,100.);
    RebinHistAbove(d_pass_l1,100.);*/
  }
  auto l1_da = divide(d_pass_l1, d_sample_l1, "efficiency");
  auto tf_da = divide(d_pass_tof, d_sample_tof, "efficiency");
  if (trgType=="n") {
    tr_da = divide(d_phys_trig,d_unb_trig , "trigger");
  }
  if (trgType=="e") tr_da = emulateTrigger(charge,start_bin,stop_bin,geom,rebin);
  daq_da = emulateDAQ(charge,start_bin,stop_bin,geom,rebin);
  auto tk_da = divide(d_pass_tr, d_sample_tr, "efficiency");
  auto l1u_da= divide(d_pass_l1u,d_sample_l1u,"efficiency");
  auto tkch_da=divide(d_pass_trch,d_sample_trch,"efficiency");
  outfile->WriteTObject(l1_da,"l1_da");
  outfile->WriteTObject(tf_da,"tf_da");
  outfile->WriteTObject(tr_da,"tr_da");
  outfile->WriteTObject(tk_da,"tk_da");
  outfile->WriteTObject(l1u_da,"l1u_da");
  outfile->WriteTObject(tkch_da,"tkch_da");
  outfile->WriteTObject(daq_da,"daq_da");
  outfile->WriteTObject(l1ch,"l1ch");
//-------Get the TTree with passed events
  double gen;
  /*float den_reco_il1, den_reco_inn, den_reco_beta;
  float den_ltof_reco_il1, den_ltof_reco_inn, den_ltof_reco_beta;
  float den_ltof_l9Unb_05_reco_il1, den_ltof_l9Unb_05_reco_inn, den_ltof_l9Unb_05_reco_beta;
  float den_ltof_l9Unb_1_reco_inn, den_ltof_l9Unb_1_reco_il1, den_ltof_l9Unb_1_reco_beta;
  short den_passed;
  short den_ltof_passed, den_ltof_l9Unb_05_passed, den_ltof_l9Unb_1_passed;
  unsigned int den_species,den_ltof_species,den_ltof_l9Unb_05_species,den_ltof_l9Unb_1_species;*/
  Sample den, den_ltof, den_ltof_l9Unb_05, den_ltof_l9Unb_1;
  double ngen;
  double w;
  TString t1,t2,t3,t4, rwcc;
  if (mcType=="s") {t1 = Form("IL1/tree1_Z%u_IL1",charge);
                 t2 = Form("IL1/tree2_Z%u_IL1",charge);
                 t3 = Form("IL1/tree3_Z%u_IL1",charge);
                 t4 = Form("IL1/tree4_Z%u_IL1",charge);
                 rwcc = "IL1/mc_samp";} if (mcType=="g")
                 {t1= "reco_gen1"; t2 = "reco_gen2"; t3 = "reco_gen3"; t4 = "reco_gen4";}
                 if (treeType=="old") {t1= "reco_gen1"; t2 = "reco_gen2"; t3 = "reco_gen3"; t4 = "reco_gen4"; rwcc = "mc_samp";}
  std::cout << "[INFO] : trees stuff:\n"
            << "tree1 name: "<< t1 <<"\n"
            << "tree2 name: "<< t2 <<"\n"
            << "tree3 name: "<< t3 << "\n"
            << "tree4 name: "<< t4 << "\n"
            << std::endl;
  auto tree1 = (TTree*)r.passed->Get( t1.Data() );
  auto tree2 = (TTree*)r.passed->Get( t2.Data() );
  auto tree3 = (TTree*)r.passed->Get( t3.Data() );
  auto tree4 = (TTree*)r.passed->Get( t4.Data() );
  ngen= ((TH1D*)r.rawacc->Get(rwcc.Data()))->GetSumOfWeights();
  tree1->SetBranchAddress("gen", &gen);
  tree1->SetBranchAddress("den_species",&den.species);
  tree1->SetBranchAddress("den_p", &den.passed);
  tree1->SetBranchAddress("den_inn", &den.reco_inn);
  tree1->SetBranchAddress("den_il1", &den.reco_il1);
  tree1->SetBranchAddress("den_b", &den.reco_beta);

  tree2->SetBranchAddress("den_ltof_species",&den_ltof.species);
  tree2->SetBranchAddress("den_ltof_p", &den_ltof.passed);
  tree2->SetBranchAddress("den_ltof_inn", &den_ltof.reco_inn);
  tree2->SetBranchAddress("den_ltof_il1", &den_ltof.reco_il1);
  tree2->SetBranchAddress("den_ltof_b", &den_ltof.reco_beta);

  tree3->SetBranchAddress("den_ltof_l9Unb_05_species",&den_ltof_l9Unb_05.species);
  tree3->SetBranchAddress("den_ltof_l9Unb_05_p", &den_ltof_l9Unb_05.passed);
  tree3->SetBranchAddress("den_ltof_l9Unb_05_inn", &den_ltof_l9Unb_05.reco_inn);
  tree3->SetBranchAddress("den_ltof_l9Unb_05_il1", &den_ltof_l9Unb_05.reco_il1);
  tree3->SetBranchAddress("den_ltof_l9Unb_05_b", &den_ltof_l9Unb_05.reco_beta);

  tree4->SetBranchAddress("den_ltof_l9Unb_1_species",&den_ltof_l9Unb_1.species);
  tree4->SetBranchAddress("den_ltof_l9Unb_1_p", &den_ltof_l9Unb_1.passed);
  tree4->SetBranchAddress("den_ltof_l9Unb_1_inn", &den_ltof_l9Unb_1.reco_inn);
  tree4->SetBranchAddress("den_ltof_l9Unb_1_il1", &den_ltof_l9Unb_1.reco_il1);
  tree4->SetBranchAddress("den_ltof_l9Unb_1_b", &den_ltof_l9Unb_1.reco_beta);

  std::cout << "Entries tree 1 for L1, L1Unb, Track, Triggerr : " << tree1->GetEntries() << std::endl;
  std::cout << "Entries tree 2 for UToF :                       " << tree2->GetEntries() << std::endl;
  std::cout << "Entries tree 3 for TrackCharge :                " << tree3->GetEntries() << std::endl;
  std::cout << "Entries tree 4 :                                " << tree4->GetEntries() << std::endl;

  //-------Read published flux-----
  TH1D *pub_flux;
  pub_flux = GetQYanFlux(charge);
  auto flux = (TH1D*)pub_flux->Clone("flux");
  MultiplyByXPower(flux,-2.7);
  outfile->WriteTObject(pub_flux,"pub_flux");
//-------Unfolding procedure------------
  auto unf_factor = new TH1D();
  auto final_acc  = new TH1D();
//-------Creating model-----------------
  Z    = charge;
  A    = charge*2;
  mass = 0.938272075*A;

//---------Flux maps for global MC reweighing
  // supponiamo di avere già un istogramma flux_B, flux_C, flux_O
  // (uno per specie: B=5, C=6, O=8, ecc.)
  std::map<unsigned int, TF1*> flux_models;
  std::map<unsigned int, double> norms;

  std::map<unsigned int, TH1*> flux_hists;
  std::vector<short> allowedZ = {14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  if (mcType=="g") {
    for (const auto f : allowedZ) {
      if (f == 21 || f == 22 || f == 23 || f == 24 || f == 25)
        flux_hists[f] = buildDummyFlux(23,"13.5y","");
      else {
        auto flux = GetQYanFlux(f);
        MultiplyByXPower(flux,-2.7);
        flux_hists[f] = flux;
      }
      for (auto const& kv : flux_hists) {
        short Z = kv.first;
        TH1* flux_hist = kv.second;
        TF1* f = build_flux_model(Z, flux_hist, prod_min, prod_max);
        flux_models[Z] = f;
        norms[Z] = f->Integral(prod_min, prod_max);
      }
    }
  }

  //-------Set fit limits
  auto l1_range = SetFitLimits(SplineUtility::Efficiency::L1Eff,charge);
  auto tf_range = SetFitLimits(SplineUtility::Efficiency::TofEff,charge);
  auto tr_range = SetFitLimits(SplineUtility::Efficiency::TriggerEff,charge);
  auto tk_range = SetFitLimits(SplineUtility::Efficiency::TrackEff,charge);
  auto l1u_range= SetFitLimits(SplineUtility::Efficiency::L1UnbEff,charge);
  auto tkch_range=SetFitLimits(SplineUtility::Efficiency::TrackChEff,charge);
  auto daq_range =SetFitLimits(SplineUtility::Efficiency::DaqEff,charge);
  auto acc_range= SetFitLimits(SplineUtility::Efficiency::Acc,charge);
  auto l1ch_range=SetFitLimits(SplineUtility::Efficiency::L1ChargeCut,charge);

  auto l1_knots= SetKnots(SplineUtility::Efficiency::L1Eff,charge);
  auto tf_knots= SetKnots(SplineUtility::Efficiency::TofEff,charge);
  auto tr_knots= SetKnots(SplineUtility::Efficiency::TriggerEff,charge);
  auto tk_knots= SetKnots(SplineUtility::Efficiency::TrackEff,charge);
  auto l1u_knots= SetKnots(SplineUtility::Efficiency::L1UnbEff,charge);
  auto tkch_knots=SetKnots(SplineUtility::Efficiency::TrackChEff,charge);
  auto daq_knots= SetKnots(SplineUtility::Efficiency::DaqEff,charge);
  auto acc_knots =SetKnots(SplineUtility::Efficiency::Acc,charge);
  auto l1ch_knots=SetKnots(SplineUtility::Efficiency::L1ChargeCut,charge);

  printf("[INFO] Knots for L1 pickup fitting : %i - %i \n",l1_knots.first, l1_knots.second);
  printf("[INFO] Knots for Tof fitting : %i - %i \n",tf_knots.first, tf_knots.second);
  printf("[INFO] Knots for Trigger fitting : %i - %i \n",tr_knots.first, tr_knots.second);
  printf("[INFO] Knots for Track fitting : %i - %i \n",tk_knots.first,tk_knots.second);
  printf("[INFO] Knots for L1 detect fitting : %i - %i \n",l1u_knots.first, l1u_knots.second);
  printf("[INFO] Knots for Track charge fitting : %i - %i \n",tkch_knots.first,tkch_knots.second);
  printf("[INFO] Knots for Daq fitting : %i - %i \n",daq_knots.first,daq_knots.second);
  printf("[INFO] Knots for Acceptance fitting : %i - %i \n",acc_knots.first,acc_knots.second);

  printf("[INFO] Range for Daq fitting : %f - %f GV\n",daq_range.first,daq_range.second);
  final_daq  = autospline(daq_da,daq_range.first,daq_range.second,daq_knots.first,daq_knots.second,maxDaqFit(charge));   
  final_l1ch = autospline(l1ch,l1ch_range.first,l1ch_range.second,l1ch_knots.first,l1ch_knots.second);
  outfile->WriteTObject(final_daq,"final_daq");
  outfile->WriteTObject(final_l1ch,"final_l1ch");


  std::vector<double> xvals, yerr68;
  for (int i = 0; i < nRbins_HighZ - 1; ++i) {
       double x_center = 0.5 * (Rbins_HighZ[i] + Rbins_HighZ[i+1]);
    xvals.push_back(x_center);
  }
  std::vector<double> knots = {1,3,4,6,8,10,30,100,300};

  for (int iter=0; iter<50; iter++) {
        printf("[UNFOLDING] Currently doing unfolding iteration %i --- \n", iter);
        auto mc_samp     = (TH1D*)hist_log->Clone("mc_samp");
        auto mc_pass     = (TH1D*)hist_log->Clone("mc_pass");
        auto mc_pass_gen = (TH1D*)hist_log->Clone("mc_pass_gen");
        /*mc_samp->Rebin(9);
        mc_pass->Rebin(9);
        mc_pass_gen->Rebin(9);*/
        auto samp_l1_mc  = (TH1D*)hist_rig_highZ->Clone("samp_l1_mc");
        auto samp_tf_mc  = (TH1D*)hist_rig_highZ->Clone("samp_tf_mc");
        auto samp_tr_mc  = (TH1D*)hist_rig_highZ->Clone("samp_tr_mc");
        auto samp_tk_mc  = (TH1D*)hist_rig_highZ->Clone("samp_tk_mc");
        auto samp_l1u_mc = (TH1D*)hist_rig_highZ->Clone("samp_l1u_mc");
        auto samp_tkch_mc= (TH1D*)hist_rig_highZ->Clone("samp_tkch_mc");
        auto pass_l1_mc  = (TH1D*)hist_rig_highZ->Clone("pass_l1_mc");
        auto pass_tf_mc  = (TH1D*)hist_rig_highZ->Clone("pass_tf_mc");
        auto pass_tr_mc  = (TH1D*)hist_rig_highZ->Clone("pass_tr_mc");
        auto pass_tk_mc  = (TH1D*)hist_rig_highZ->Clone("pass_tk_mc");
        auto pass_l1u_mc = (TH1D*)hist_rig_highZ->Clone("pass_l1u_mc");
        auto pass_tkch_mc= (TH1D*)hist_rig_highZ->Clone("pass_tkch_mc");

//------Build Flux model to reweigh mc
        //auto flux_model_temp = spfit(flux, 3, prod_min, prod_max, knots, "", &xvals, &yerr68);
        //auto flux_model = MultiplyTF1ByXPower(flux_model_temp, -2.7);
        auto flux_model = new TF1("flux_model",force_field,0.8,prod_max,3);
        flux_model->SetParameters(20,0.6,-2.7);
        flux->Fit(flux_model);
        flux_model->SetNpx(100000);
        double chi2 = flux_model->GetChisquare();
        int npoints = flux_model->GetNumberFitPoints();
        int npar = flux_model->GetNumberFreeParameters();
        int ndf = npoints - npar;
        auto title = Form("Chi2/ndf: %f", chi2 / ndf);
        flux_model->SetTitle(title);
        flux_model->SetName(title);
        //MultiplyByXPower(flux,-2.7);      
        outfile->WriteTObject(flux_model, Form("flux_model_it%i", iter));
        outfile->WriteTObject(flux, Form("flux_it%i", iter));
//------Printing fit result
        printf("Chi2 / ndf = %.2f / %d = %.2f\n",chi2,ndf,chi2/ndf);
        //double norm = integrate(flux_model,prod_min,prod_min); //ORIGINALE
        double norm = flux_model->Integral(prod_min,prod_max);
        for(int i=1; i<=mc_samp->GetNbinsX(); i++) {
            //mc_samp->SetBinContent(i, integrate(flux_model, mc_samp->GetBinLowEdge(i), mc_samp->GetBinLowEdge(i+1))); //ORIGINALE
            mc_samp->SetBinContent(i, flux_model->Integral(mc_samp->GetBinLowEdge(i), mc_samp->GetBinLowEdge(i+1)));
            mc_samp->SetBinError(i , sqrt(mc_samp->GetBinContent(i)));
        }
        mc_samp->Scale(ngen/norm);
        if (iter==0) outfile->WriteTObject(mc_samp,"mc_samp");
        int entries1 = tree1->GetEntries();
        int entries2 = tree2->GetEntries();
        int entries3 = tree3->GetEntries();
        int entries4 = tree4->GetEntries();
//------Reweighing sample
        double perc=0;
        for(int event = 0; event < entries1; ++event) {
          tree1->GetEntry(event);
          Double_t percc = 100.0*( (event+1.0)/entries1);
          if (percc>=perc) {
            printf("\r[INFO] Reweighing mc acceptance: %d%%", (int)(100.0*( (event+1.0)/entries1)));
            std::cout << std::flush;
            perc++;
          }
          double weight = flux_model->Eval(gen)*gen*TMath::Log(prod_max/prod_min) / norm;
          if( (den.passed & 1) == 1) {
            mc_pass->Fill(den.reco_il1, weight);
            mc_pass_gen->Fill(gen, weight);
          }
        } std::cout << std::endl;
  //------L1,L1Unb,Track,Trigger
          perc=0;
          for(int event = 0; event < entries1; event++ ) {
            tree1->GetEntry(event);
            Double_t percc = 100.0*( (event+1.0)/entries1);
            if (percc>=perc) {
              printf("\r[INFO] Reweighing L1, L1Unb, Track, Trigger with <den> sample: %d%%", (int)(100.0*( (event+1.0)/entries1)));
              std::cout << std::flush;
              perc++;
            }
            if (mcType=="g") {
              // prendo il modello corretto in base alla specie
              if (den.species > 26 || den.species==0) continue;
              TF1* flux_model_event   = flux_models[den.species];
              double norm_event       = norms[den.species];
              // calcolo peso
              double weight = flux_model_event->Eval(gen) * gen
                            * TMath::Log(prod_max/prod_min) / norm_event;
              w = weight;
            }
            else { 
              double weight = flux_model->Eval(gen)*gen*TMath::Log(prod_max/prod_min) / norm;
              w = weight;
            }
            if( den.passed & (1<<1))  samp_l1_mc->Fill(den.reco_inn , w);
            if( den.passed & (1<<2))  pass_l1_mc->Fill(den.reco_inn , w);
            if( den.passed & (1<<5))  samp_tr_mc->Fill(den.reco_il1 , w);
            if( den.passed & (1<<6))  pass_tr_mc->Fill(den.reco_il1 , w);
            if( den.passed & (1<<7))  samp_tk_mc->Fill(den.reco_beta, w);
            if( den.passed & (1<<8))  pass_tk_mc->Fill(den.reco_beta, w);
            if( den.passed & (1<<9))  samp_l1u_mc->Fill(den.reco_inn, w);  
            if( den.passed & (1<<10)) pass_l1u_mc->Fill(den.reco_inn, w);
          } std::cout << std::endl;
  //------UToF
          perc=0;
          for(int event = 0; event < entries2; event++ ) {
            tree2->GetEntry(event);
            Double_t percc = 100.0*( (event+1.0)/entries2);
            if (percc>=perc) {
              printf("\r[INFO] Reweighing UToF with <den_ltof> sample: %d%%", (int)(100.0*( (event+1.0)/entries2)));
              std::cout << std::flush;
              perc++;
            }
            if (mcType=="g") {
              // prendo il modello corretto in base alla specie
              if (den_ltof.species > 26 || den_ltof.species==0) continue;
              TF1* flux_model_event   = flux_models[den_ltof.species];
              double norm_event       = norms[den_ltof.species];
              // calcolo peso
              double weight = flux_model_event->Eval(gen) * gen
                            * TMath::Log(prod_max/prod_min) / norm_event;
              w = weight;
            }
            else {
              double weight = flux_model->Eval(gen)*gen*TMath::Log(prod_max/prod_min) / norm;
              w = weight;
            }
            if( den_ltof.passed & (1<<3))  samp_tf_mc->Fill(den_ltof.reco_il1 , w);
            if( den_ltof.passed & (1<<4))  pass_tf_mc->Fill(den_ltof.reco_il1 , w);
          } std::cout << std::endl;
  //------TrackCharge
          perc=0;
          for(int event = 0; event < entries4; event++ ) {
            tree4->GetEntry(event);
            Double_t percc = 100.0*( (event+1.0)/entries4);
            if (percc>=perc) {
              printf("\r[INFO] Reweighing TrackCharge with <den_ltof_l9_1> sample: %d%%", (int)(100.0*( (event+1.0)/entries4)));
              std::cout << std::flush;
              perc++;
            }
            if (mcType=="g") {
              // prendo il modello corretto in base alla specie
              if (den_ltof_l9Unb_05.species > 26 || den_ltof_l9Unb_05.species==0) continue;
              TF1* flux_model_event   = flux_models[den_ltof_l9Unb_05.species];
              double norm_event       = norms[den_ltof_l9Unb_05.species];
              // calcolo peso
              double weight = flux_model_event->Eval(gen) * gen
                            * TMath::Log(prod_max/prod_min) / norm_event;
              w = weight;
            }
            else {
              double weight = flux_model->Eval(gen)*gen*TMath::Log(prod_max/prod_min) / norm;
              w = weight;
            }
            if( den_ltof_l9Unb_1.passed & (1<<11)) samp_tkch_mc->Fill(den_ltof_l9Unb_1.reco_il1, w);
            if( den_ltof_l9Unb_1.passed & (1<<12)) pass_tkch_mc->Fill(den_ltof_l9Unb_1.reco_il1, w);   
          } std::cout << std::endl;
  //------Reweigh mc track using rigidity map
          if (rigmap=="rigmap" && iter==0) {
            printf("\n[INFO] Reweighing track mc with rigidity map.. \n");
            RigidityMap(charge,samp_tk_mc);
            RigidityMap(charge,pass_tk_mc);
          }
          std::cout << "[DEBUG] pass_l1_mc directory: " << pass_l1_mc->GetDirectory() << std::endl;
          outfile->WriteTObject(mc_pass, Form("mc_pass_it%i",iter) );
          outfile->WriteTObject(mc_pass_gen, Form("mc_pass_gen_it%i",iter) );
  //------Rebin mc efficiencies
          if (rebin=="rebin") {
            printf("[INFO] Rebinning trigger mc efficiency.. \n");
            RebinByEfficiencyOverall(samp_tr_mc,SplineUtility::Efficiency::TriggerEff,charge);
            RebinByEfficiencyOverall(pass_tr_mc,SplineUtility::Efficiency::TriggerEff,charge);
            printf("[INFO] Rebinning l1 detect mc efficiency.. \n");
            RebinByEfficiencyOverall(samp_l1u_mc,SplineUtility::Efficiency::L1UnbEff,charge);
            RebinByEfficiencyOverall(pass_l1u_mc,SplineUtility::Efficiency::L1UnbEff,charge);
            printf("[INFO] Rebinning track charge mc efficiency.. \n");
            RebinByEfficiencyOverall(samp_tkch_mc,SplineUtility::Efficiency::TrackChEff,charge);
            RebinByEfficiencyOverall(pass_tkch_mc,SplineUtility::Efficiency::TrackChEff,charge);
            printf("[INFO] Rebinning track mc efficiency.. \n");
            RebinByEfficiencyOverall(samp_tk_mc,SplineUtility::Efficiency::TrackEff,charge);
            RebinByEfficiencyOverall(pass_tk_mc,SplineUtility::Efficiency::TrackEff,charge);
            printf("[INFO] Rebinning tof mc efficiency.. \n");
            RebinByEfficiencyOverall(samp_tf_mc,SplineUtility::Efficiency::TofEff,charge);
            RebinByEfficiencyOverall(pass_tf_mc,SplineUtility::Efficiency::TofEff,charge);
            printf("[INFO] Rebinning l1 mc efficiency.. \n");
            RebinByEfficiencyOverall(samp_l1_mc,SplineUtility::Efficiency::L1Eff,charge);
            RebinByEfficiencyOverall(pass_l1_mc,SplineUtility::Efficiency::L1Eff,charge);
            /*RebinHistAbove(samp_l1_mc,100.);
            RebinHistAbove(pass_l1_mc,100.);*/
          }
  //------Build current iteration efficiencies and corrections;
          l1_mc = divide(pass_l1_mc, samp_l1_mc, "alberto"); 
          tf_mc = divide(pass_tf_mc, samp_tf_mc, "alberto");
          tr_mc = divide(pass_tr_mc, samp_tr_mc, "alberto");
          tk_mc = divide(pass_tk_mc, samp_tk_mc, "alberto");      
          l1u_mc= divide(pass_l1u_mc,samp_l1u_mc,"alberto");      
          tkch_mc=divide(pass_tkch_mc,samp_tkch_mc,"alberto");
          outfile->WriteTObject(l1_mc, Form("l1_mc_it%i",iter) );
          outfile->WriteTObject(tf_mc, Form("tf_mc_it%i",iter) );
          outfile->WriteTObject(tr_mc, Form("tr_mc_it%i",iter) );
          outfile->WriteTObject(tk_mc, Form("tk_mc_it%i",iter) );       
          outfile->WriteTObject(l1u_mc,Form("l1u_mc_it%i",iter));       
          outfile->WriteTObject(tkch_mc,Form("tkch_mc_it%i",iter));
  //------Reweighing efficiencies
          damc_l1 = divide(l1_da, l1_mc, "ratio");
          damc_tf = divide(tf_da, tf_mc, "ratio");
          damc_tr = divide(tr_da, tr_mc, "ratio");
          damc_tk = divide(tk_da, tk_mc, "ratio");      
          damc_l1u= divide(l1u_da,l1u_mc,"ratio");     
          damc_tkch=divide(tkch_da,tkch_mc,"ratio");
          outfile->WriteTObject(damc_l1, Form("damc_l1_it%i",iter) );
          outfile->WriteTObject(damc_tf, Form("damc_tf_it%i",iter) );
          outfile->WriteTObject(damc_tr, Form("damc_tr_it%i",iter) );
          outfile->WriteTObject(damc_tk, Form("damc_tk_it%i",iter) );       
          outfile->WriteTObject(damc_l1u,Form("damc_l1u_it%i",iter));
          outfile->WriteTObject(damc_tkch,Form("damc_tkch_it%i",iter));
          printf("[INFO] Range for L1 pickup fitting : %f - %f GV \n",l1_range.first, l1_range.second);
          printf("[INFO] Range for Tof fitting : %f - %f GV \n",tf_range.first, tf_range.second);
          printf("[INFO] Range for Trigger fitting : %f - %f GV \n",tr_range.first, tr_range.second);
          printf("[INFO] Range for Track fitting : %f - %f GV\n",tk_range.first,tk_range.second);
          printf("[INFO] Range for L1 detect fitting : %f - %f GV \n",l1u_range.first, l1u_range.second);
          printf("[INFO] Range for Track charge fitting : %f - %f GV\n",tkch_range.first,tkch_range.second);
          printf("[INFO] Range for Acceptance fitting : %f - %f GV\n",acc_range.first,acc_range.second);  
          printf("[DEBUG] damc_l1 entries: %.0f\n", damc_l1->GetEntries());
          printf("[DEBUG] damc_tf entries: %.0f\n", damc_tf->GetEntries());
          printf("[DEBUG] damc_tr entries: %.0f\n", damc_tr->GetEntries());
          printf("[DEBUG] damc_tk entries: %.0f\n", damc_tk->GetEntries());
          printf("[DEBUG] damc_l1u entries: %.0f\n", damc_l1u->GetEntries());
          printf("[DEBUG] damc_tkch entries: %.0f\n", damc_tkch->GetEntries());
          auto check_hist = [&](TH1D* h, const char* name, double xmin, double xmax) {
          int bmin = h->FindBin(xmin);
          int bmax = h->FindBin(xmax);
          int nonzero = 0;
          for (int b = bmin; b <= bmax; ++b) {
              double c = h->GetBinContent(b);
              double e = h->GetBinError(b);
              if (c != 0 && std::isfinite(c) && e >= 0) ++nonzero;
          }
          printf("[CHECK] %s: bins %d-%d (%.2f-%.2f GV), nonzero = %d / %d\n",
                name, bmin, bmax, xmin, xmax, nonzero, bmax - bmin + 1);
          };
          check_hist(l1_da, "l1_da",  l1_range.first, l1_range.second);
          check_hist(l1_mc, "l1_mc",  l1_range.first, l1_range.second);
          check_hist(tr_da, "tr_da",  tr_range.first, tr_range.second);
          check_hist(tr_mc, "tr_mc",  tr_range.first, tr_range.second);
          check_hist(damc_l1, "damc_l1", l1_range.first, l1_range.second);
          check_hist(damc_tf, "damc_tf", tf_range.first, tf_range.second);
          check_hist(damc_tr, "damc_tr", tr_range.first, tr_range.second);
          check_hist(damc_tk, "damc_tk", tk_range.first, tk_range.second);
          check_hist(damc_l1u, "damc_l1u", l1u_range.first, l1u_range.second);
          check_hist(damc_tkch, "damc_tkch", tkch_range.first, tkch_range.second);
          final_damc_l1 = autospline(damc_l1, l1_range.first, l1_range.second, l1_knots.first, l1_knots.second);    
          final_damc_tf = autospline(damc_tf, tf_range.first, tf_range.second, tf_knots.first, tf_knots.second);     
          final_damc_tr = autospline(damc_tr, tr_range.first, tr_range.second, tr_knots.first, tr_knots.second);      
          final_damc_tk = autospline(damc_tk, tk_range.first, tk_range.second, tk_knots.first, tk_knots.second);     
          final_damc_l1u= autospline(damc_l1u,l1u_range.first,l1u_range.second,l1u_knots.first,l1u_knots.second);
          final_damc_tkch=autospline(damc_tkch,tkch_range.first,tkch_range.second,tkch_knots.first,tkch_knots.second);  
          outfile->WriteTObject(final_damc_l1, Form("final_damc_l1_it%i", iter));
          outfile->WriteTObject(final_damc_tf, Form("final_damc_tf_it%i", iter));
          outfile->WriteTObject(final_damc_tr, Form("final_damc_tr_it%i", iter));
          outfile->WriteTObject(final_damc_tk, Form("final_damc_tk_it%i", iter));
          outfile->WriteTObject(final_damc_l1u, Form("final_damc_l1u_it%i", iter));
          outfile->WriteTObject(final_damc_tkch,Form("final_damc_tkch_it%i",iter));

          if (charge!=15) {
            final_damc_tot = (TH1D*)final_damc_l1->Clone();
            final_damc_tot->Multiply(final_damc_tf);
            final_damc_tot->Multiply(final_damc_tr);
            final_damc_tot->Multiply(final_damc_tk);     
            final_damc_tot->Multiply(final_damc_l1u);
            final_damc_tot->Multiply(final_damc_tkch);  
            final_damc_tot->Multiply(final_daq);
            final_damc_tot->Multiply(final_l1ch);
        } else if (charge==15) { // If charge == 15 retrieve the total correction from ../P/Interpolation
          printf("[INFO] Avoiding data/mc and retrieving total correction from interpolation\n");
          TString intPath = "../IonsSelected/"+ionPath+"/Interpolation/"+ionPath+
                       Form("_interpolation%s_%s.root",sec_track_out.Data(),timePeriod.Data());   
          TFile *ff = new TFile(intPath.Data(), "READ");
          final_damc_tot = (TH1D*)ff->Get("final_damc_tot");
          final_damc_tot->Multiply(final_l1ch);
        }
//------Build current iteration acceptance
        acc = divide(mc_pass, mc_samp, "efficiency");
        acc->Scale(TMath::Pi() * 3.9 * 3.9);
        spline_acc = autospline(acc, acc_range.first, acc_range.second, acc_knots.first,acc_knots.second);
        final_acc = (TH1D*)spline_acc->Clone();
        final_acc->Multiply(final_damc_tot);

//------Dumping corrections values fot 10 GV---------
        if (iter==0) {
          std::cout << "[ --- DUMP ---] Total correction value @ 10 GV: " << final_damc_tot->GetBinContent(final_damc_tot->FindBin(10)) << std::endl;
        }

//------Build current iteration flux----------
        auto oldflux = (TH1D*)flux->Clone();
        flux = (TH1D*)counts->Clone();
        flux->Divide(lvt_25);
        flux->Divide(hwidth);
        flux->Divide(final_acc);
        unf_factor = divide(mc_pass, mc_pass_gen, "simple");
//------Write current iteration quantities
        outfile->WriteTObject(final_damc_tot, Form("final_damc_tot_it%i", iter));
        outfile->WriteTObject(acc, Form("acc_it%i", iter));
        outfile->WriteTObject(spline_acc, Form("spline_acc_it%i", iter));
        outfile->WriteTObject(final_acc, Form("final_acc_it%i", iter));
        outfile->WriteTObject(unf_factor, Form("unf_factor_it%i", iter));
//------Unfolding convergence criterion---------
        if(iter>0) {
            double max_rel_diff = 0;
            for(int i=prod_min; i<=flux->FindBin(prod_max); ++i) {
              double old_value = oldflux->GetBinContent(i);
              if (old_value != 0) {
                  auto diff = round(10000 * fabs(flux->GetBinContent(i) - old_value) / old_value) / 100;
                  if (diff > max_rel_diff) max_rel_diff = diff;
              }
            }
            std::cout<<"[UNFOLDING] Iteration: "<<iter<<", Max Relative Diff: "<<max_rel_diff<<"%"<<std::endl;
            if(max_rel_diff<=0.01 && iter>0) break;
        }
//------Rebuild the flux multiplied for next fit iteration
        //MultiplyByXPower(flux,2.7);
        outfile->WriteTObject(pass_tk_mc,"pass_tk_mc");
        outfile->WriteTObject(samp_tk_mc,"samp_tk_mc");
    }//iter
//------Writing final acceptance to use in BuildFlux
  outfile->WriteTObject(unf_factor,"unf_factor");
  outfile->WriteTObject(final_acc,"final_acc");
//------Writing efficiencies and corrections to be printed (only for Si and S)
  outfile->WriteTObject(l1_mc,"l1_mc");
  outfile->WriteTObject(tf_mc,"tf_mc");
  outfile->WriteTObject(tr_mc,"tr_mc");
  outfile->WriteTObject(tk_mc,"tk_mc");
  outfile->WriteTObject(l1u_mc,"l1u_mc");
  outfile->WriteTObject(tkch_mc,"tkch_mc");
  outfile->WriteTObject(damc_l1,"damc_l1");
  outfile->WriteTObject(damc_tf,"damc_tf");
  outfile->WriteTObject(damc_tr,"damc_tr");
  outfile->WriteTObject(damc_tk,"damc_tk");
  outfile->WriteTObject(damc_l1u,"damc_l1u");
  outfile->WriteTObject(damc_tkch,"damc_tkch");
  outfile->WriteTObject(final_damc_l1,"final_damc_l1");
  outfile->WriteTObject(final_damc_tf,"final_damc_tf");
  outfile->WriteTObject(final_damc_tr,"final_damc_tr");
  outfile->WriteTObject(final_damc_tk,"final_damc_tk");
  outfile->WriteTObject(final_damc_l1u,"final_damc_l1u");
  outfile->WriteTObject(final_damc_tkch,"final_damc_tkch");
  outfile->WriteTObject(final_damc_tot,"final_damc_tot");
//------Writing acceptance to be printed
  outfile->WriteTObject(acc,"acc");
  outfile->WriteTObject(spline_acc,"spline_acc");
//------Closing the output file
  outfile->Close();
}

//----Definitions------
TH1D *BuildPublishedFlux(unsigned int charge) {
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> values;
    std::vector<double> syst;
    std::vector<double> errors;
    readFluxTable(lower_bounds,upper_bounds,values,errors,charge,syst);
    auto published = buildPublished(lower_bounds,upper_bounds,values,errors);
    return published;
}
TGraph *build_histo_model(TH1D *histo) {
  int prod_Max = 550;
  auto model = new TGraphErrors(autospline(histo, prod_min, prod_Max,7,12));
  model->SetBit(TGraph::kIsSortedX);
  auto flux_model = new TGraphErrors();
  double y;
  for(double x=prod_min; x<=prod_Max; x=pow(10,log10(x)+0.01)){
    if(x<prod_Max) y = model->Eval(x, 0, "S");
    flux_model->SetPoint(flux_model->GetN(), x, y);
  }
  flux_model->SetBit(TGraph::kIsSortedX);
  return flux_model;
}
double integrate(const TGraph* fun, double min, double max) {
  double integral = 0;
  const double step = 0.00001;
  for (double x = min; x < max; x += step)
    integral += (fun->Eval(x) + 4*fun->Eval(x+step/2) + fun->Eval(x+step)) * step/6;
  return integral;
}
double integrate(const TF1* fun, double min, double max) {
  double integral = 0;
  const double step = 0.00001;
  for (double x = min; x < max; x += step)
    integral += (fun->Eval(x) + 4*fun->Eval(x+step/2) + fun->Eval(x+step)) * step/6;
  return integral;
}
void RigidityMap(unsigned int charge, TH1D *h) {
  //apri file e prendi rigamp
  TFile *file = new TFile("../IonsSelected/"+getIonPath(charge)+"/RigidityMap/rigmap.root" );
  auto rigmap = (TH2D*)file->Get("map_normalized");
  if (!rigmap) std::cerr << "Map not found" <<std::endl;
  //itera sugli eventi di h
  for (int i=1; i<=h->GetNbinsX(); i++) {
    //per ogni evento prendi rigidità rig (indice i) e content
    double content = h->GetBinContent(i);
    double r_low = h->GetBinLowEdge(h->GetBinCenter(i));
    //per ogni rigidità (indice i), fai la proiezione del bin della mappa dove cade 
    if (content!=0) {
      auto proj = rigmap->ProjectionY("proj",r_low,i);
      // prendi la media
      double mean = proj->GetMean();
      // assegna lo stesso bin content alla nuova rigidità
      h->SetBinContent(h->GetBinCenter(mean),content);
      delete proj;
    }
  }
  delete rigmap;
  delete file;
}
int getEntriesBySample(int sampleNum, TTree *t1, TTree *t2, TTree *t3, TTree *t4) {
  int a =1;
  switch(sampleNum) {
    case 1:
      a = t1->GetEntries();
      break;
    case 2:
      a = t2->GetEntries();
      break;
    case 3:
      a = t3->GetEntries();
      break;
    case 4:
      a = t4->GetEntries();
      break;
  }
  return a;
}
TString SamplePrefix(int sampleNum) {
  TString a;
  switch(sampleNum) {
    case 1:
      a = "den";
      break;
    case 2:
      a = "den_ltof";
      break;
    case 3:
      a = "den_ltof_l9_05";
      break;
    case 4:
      a = "den_ltof_l9_1";
      break;
  }
  return a;
}
double maxDaqFit(unsigned int charge) {
  double xmax;
  switch(charge) {
    case 9:
      xmax = 1000;
      break;
    case 14:
      xmax = 50;
      break;
    case 15:
      xmax = 15;
      break;
    case 16:
      xmax = 10;
      break;
    case 17:
      xmax = 15;
      break;
    case 18:
      xmax = 20;
      break;
    case 19:
      xmax = 20;
      break;
    case 20:
      xmax = 10;
      break;
    case 21:
      xmax = 10;
      break;
    case 22:
      xmax = 10;
      break;
    case 23:
      xmax = 9;
      break;
    case 24:
      xmax = 9;
      break;
    case 25:
      xmax = 8;
      break;
    case 26:
      xmax = 5;
      break;
    default:
      xmax = 10;
      break;
  }
  return xmax;
}
TH1D *buildDummyFlux(unsigned int charge, TString timePeriod, TString inp_sec_track) {
    //Getting counts
    TString fileName= Form("../IonsSelected/"+getIonPath(charge)+"/Counts/%scounts.root",inp_sec_track.Data());
    TFile *file = TFile::Open(fileName.Data());
    if (!file || file->IsZombie()) {
        printf("Errore nell'aprire il file\n");
    }
    int start_bin=getTimePeriod(timePeriod).first, stop_bin=getTimePeriod(timePeriod).second;
    auto counts = (TH1D *)((TH2D *)file->Get("IL1/rigidity"))->ProjectionY("r", start_bin, stop_bin);
    if (!counts) printf("Errore \n");
    //Getting livetime
    TString fileName1= "../IonsSelected/"+getIonPath(charge)+"/Livetime/livetime.root";
    TFile *file1 = TFile::Open(fileName1.Data());
    TH1D *lvt = (TH1D *)((TH2D *)file1->Get("lvt_25"))->ProjectionY("q", start_bin, stop_bin);
    //Getting Final Raw acceptance
    TString fileName2= Form("../IonsSelected/"+getIonPath(charge)+"/RawAcceptance/rawacc.root");
    TFile *file2 = TFile::Open(fileName2.Data());
    auto mc_pass_gen = (TH1D*)file2->Get("mc_pass_gen");
    auto mc_samp = (TH1D*)file2->Get("mc_samp");
    mc_pass_gen->Rebin(9);
    mc_samp->Rebin(9);
    auto acc = divide(mc_pass_gen, mc_samp, "efficiency");
    acc->Scale(TMath::Pi() * 3.9 * 3.9);
    auto acc_range= SetFitLimits(SplineUtility::Efficiency::Acc,charge);
    auto acc_knots =SetKnots(SplineUtility::Efficiency::Acc,charge);
    auto spline_acc = autospline(acc, acc_range.first, acc_range.second,acc_knots.first,acc_knots.second);
    //Bin Width
    auto hwidth = (TH1D*)hist_rig_highZ->Clone();
    for(int i=1; i<=hist_rig_highZ->GetNbinsX(); ++i) {
      hwidth->SetBinContent(i, hwidth->GetBinWidth(i));
      hwidth->SetBinError(i, 0);
    }
    //Raw Flux
    for (int ibin = 1; ibin <= counts->GetNbinsX(); ibin++) {
        double N = counts->GetBinContent(ibin);
        double LT = lvt->GetBinContent(ibin);
        double A = spline_acc->GetBinContent(ibin);
        double dR = hwidth->GetBinContent(ibin);
        double sigma_N = counts->GetBinError(ibin);
        double sigma_LT = lvt->GetBinError(ibin);
        double sigma_A = spline_acc->GetBinError(ibin);
        double sigma_dR = hwidth->GetBinError(ibin);
        double flux = (N / (LT * A * dR));
        double sigma_flux = flux * sqrt(
            pow(sigma_N / N, 2) +
            pow(sigma_LT / LT, 2) +
            pow(sigma_A / A, 2) +
            pow(sigma_dR / dR, 2)
        );

        counts->SetBinContent(ibin, flux);
        counts->SetBinError(ibin, sigma_flux);
    }
    return counts;
}


