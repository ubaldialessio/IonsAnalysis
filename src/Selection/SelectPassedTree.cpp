#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

void SetSample(Sample &den, bool condition, int bit, float& reco_var, float rigidity, int Z, bool subcondition, int subbit);
double FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain,  NAIA::Event &event);

int main(int argc, char *argv[]) {
  TH1::SetDefaultSumw2();
	auto tree1 = new TTree();
  auto tree2 = new TTree();
  auto tree3 = new TTree();
  auto tree4 = new TTree();
  auto ngen         = new TH1D("", "", 1, 0, 1);
	auto sample_tr    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
  auto mc_samp      = (TH1D*)hist_log->Clone("mc_samp");
	auto mc_pass_gen  = (TH1D*)hist_log->Clone("mc_pass_gen");
	auto mc_pass      = (TH1D*)hist_log->Clone("mc_pass");

    // this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;
	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> <sec_track>\n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString exename=argv[0],
          infilename=argv[2],
		      outname=argv[3],
          sec_track=argv[4],
		      ionPath=getIonPath(charge);
  NAIA::NAIAChain chain;  
	if(infilename.Contains(".root") /* && filesystem::exists(infilename.Data())*/ ){
		chain.Add(infilename.Data());
	}else if (infilename.Contains(".txt") /* && filesystem::exists(infilename.Data())*/ ){
		ifstream infilelist(infilename.Data());
		TString bufname;
		while(infilelist >> bufname) 
		    if(bufname.Contains(".root" /* && filesystem::exists(bufname.Data()) */))
		chain.Add(bufname.Data());
	}
	chain.SetupBranches();
	bool isMC = chain.IsMC();
  if (sec_track=="y") outname = "sec_track_"+outname;
	TString out = StreamUtility::getOutputDir(charge,exename,outname);
  std::cout << out << std::endl;
  Sample den, den_ltof, den_ltof_l9Unb_05, den_ltof_l9Unb_1;
  float reco_il1, reco_inn, beta, lt, reco_beta;
  short species;
	double gen;
	double Amass = 1;
	double weight = 1;
	Amass = GetDataMass(charge); // Getting atomic mass number
  tree1->Branch("gen",&gen, "gen/d");
  tree1->Branch("den_species",&den.species, "species/i");
  tree1->Branch("den_p", &den.passed, "passed/S");
  tree1->Branch("den_inn", &den.reco_inn, "reco_inn/F");
  tree1->Branch("den_il1", &den.reco_il1, "reco_il1/F");
  tree1->Branch("den_b", &den.reco_beta, "reco_beta/F");

  tree2->Branch("den_ltof_species",&den_ltof.species,"species/i");
  tree2->Branch("den_ltof_p", &den_ltof.passed, "passed/S");
  tree2->Branch("den_ltof_inn", &den_ltof.reco_inn, "reco_inn/F");
  tree2->Branch("den_ltof_il1", &den_ltof.reco_il1, "reco_il1/F");
  tree2->Branch("den_ltof_b", &den_ltof.reco_beta, "reco_beta/F");

  tree3->Branch("den_ltof_l9Unb_05_species",&den_ltof_l9Unb_05.species, "species/i");
  tree3->Branch("den_ltof_l9Unb_05_p", &den_ltof_l9Unb_05.passed, "passed/S");
  tree3->Branch("den_ltof_l9Unb_05_inn", &den_ltof_l9Unb_05.reco_inn, "reco_inn/F");
  tree3->Branch("den_ltof_l9Unb_05_il1", &den_ltof_l9Unb_05.reco_il1, "reco_il1/F");
  tree3->Branch("den_ltof_l9Unb_05_b", &den_ltof_l9Unb_05.reco_beta, "reco_beta/F");

  tree4->Branch("den_ltof_l9Unb_1_species",&den_ltof_l9Unb_1. species, "species/i");
  tree4->Branch("den_ltof_l9Unb_1_p", &den_ltof_l9Unb_1.passed, "passed/S");
  tree4->Branch("den_ltof_l9Unb_1_inn", &den_ltof_l9Unb_1.reco_inn, "reco_inn/F");
  tree4->Branch("den_ltof_l9Unb_1_il1", &den_ltof_l9Unb_1.reco_il1, "reco_il1/F");
  tree4->Branch("den_ltof_l9Unb_1_b", &den_ltof_l9Unb_1.reco_beta, "reco_beta/F");
	auto mcchain = chain.GetFileInfoTree();
	auto mcinfo = new NAIA::MCFileInfo(); 
  mcchain->SetBranchAddress("MCFileInfo", &mcinfo);
  		for(unsigned long long i = 0; i < mcchain->GetEntries(); ++i) {
			mcchain->GetEntry(i);
		    double rminDC = mcinfo->GetRMin();
		    double rmaxDC = mcinfo->GetRMax();
		    for(int ibin = 1; ibin <= mc_samp->GetNbinsX(); ++ibin){
		        double rmin = mc_samp->GetBinLowEdge(ibin), rmax = mc_samp->GetBinLowEdge(ibin+1);
		        if (rmin < rminDC) rmin = rminDC;
		        if (rmax > rmaxDC) rmax = rmaxDC;
		        if (rmin > rmaxDC || rmax < rminDC) continue;
		        mc_samp->AddBinContent(ibin, mcinfo->GetNGen() * TMath::Log(rmax / rmin) / TMath::Log(rmaxDC / rminDC));
		    }
		}    
  ngen->SetBinContent(1, mcinfo->GetNGen());
  auto mysel  = Mysel(charge);
  auto myselNoL1 = MyselNoL1(charge);
  auto num_tof= Num_tof(charge);
  auto den_tof = Den_tof(charge);
  auto den_tof_ltof=Den_tof_ltof(charge);
  auto den_tof_ltof_l9_05=Den_tof_ltof_l9_05(charge);
  auto den_tof_ltof_l9_1=Den_tof_ltof_l9_1(charge);
  auto num_l1=Num_l1(charge);
  auto den_l1=Den_l1(charge);
  auto den_l1_ltof=Den_l1_ltof(charge);
  auto den_l1_ltof_l9_05=Den_l1_ltof_l9_05(charge);
  auto den_l1_ltof_l9_1=Den_l1_ltof_l9_1(charge);
  auto num_l1Unb=Num_l1Unb(charge);
  auto den_l1Unb=Den_l1Unb(charge);
  auto den_l1Unb_ltof=Den_l1Unb_ltof(charge);
  auto den_l1Unb_ltof_l9_05=Den_l1Unb_ltof_l9_05(charge);
  auto den_l1Unb_ltof_l9_1=Den_l1Unb_ltof_l9_1(charge);
  auto num_track=Num_track(charge);
  auto den_track=Den_track(charge);
  auto num_trackCh = Num_trackCh(charge);
  auto den_trackCh = Den_trackCh(charge);
  auto den_trackCh_ltof = Den_trackCh_ltof(charge);
  auto den_trackCh_ltof_l9_05 = Den_trackCh_ltof_l9_05(charge);
  auto den_trackCh_ltof_l9_1 = Den_trackCh_ltof_l9_1(charge);
  auto den_trig = Den_trig(charge);
  auto hasGoodSecondTrack = GoodSecTrack();
  auto secTrackOnDiagonal = MySel::IsSecondTrackOnDiagonal(FIT,IL1,0.3);
/////////////////// LOOP ///////////////////////////////
	for(Event& event : chain) {
        int nucleus = event.mcTruthBase->Primary.Z;
        gen = event.mcTruthBase->Primary.GetGenMomentum()/event.mcTruthBase->Primary.Z;
		    den.passed = 0;
        den_ltof.passed = 0;
        den_ltof_l9Unb_05.passed = 0;
        den_ltof_l9Unb_1.passed = 0;
        if (myselNoL1(event)) {
          if (sec_track=="y") {
            if (hasGoodSecondTrack(event)) {
              if (secTrackOnDiagonal(event)) {
                continue;
              }
            }
			    }
            den.reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
            den.passed |= 1;
            den_ltof.passed |= 1;
            den_ltof_l9Unb_05.passed |= 1;
            den_ltof_l9Unb_1.passed |= 1;
            den.species = nucleus;
            den_ltof.species = nucleus;
            den_ltof_l9Unb_05.species = species;
            den_ltof_l9Unb_1.species = species;
            //counts
            reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
            int rbin = hist_rig_highZ->FindBin(reco_il1);
            double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
            mc_pass->Fill(reco_il1);
            mc_pass_gen->Fill(gen);
        }

        //All the den samples
        SetSample(den, den_l1(event), 1, den.reco_inn, event.trTrackBase->Rigidity[FIT][INN], nucleus, num_l1(event), 2);
        SetSample(den, den_tof(event), 3, den.reco_il1, event.trTrackBase->Rigidity[FIT][IL1], nucleus, num_tof(event), 4);
        SetSample(den, den_trig(event) && (event.evSummary->PhysBPatt & 65), 5, den.reco_il1, 
                  event.trTrackBase->Rigidity[FIT][IL1], nucleus, event.evSummary->IsPhysicsTrigger(), 6);
        if (den_track(event)) {
            den.reco_beta = FillHistos(*sample_tr,charge,Amass,weight,chain,event);
            den.passed |= 1<<7;
            den.species = nucleus;
            if(num_track(event)) {
              den.passed |= 1<<8;
              den.species = nucleus;
            }
        }                  
        SetSample(den, den_l1Unb(event), 9, den.reco_inn, event.trTrackBase->Rigidity[FIT][INN], nucleus, num_l1Unb(event), 10);
        SetSample(den, den_trackCh(event), 11, den.reco_il1, event.trTrackBase->Rigidity[FIT][IL1], nucleus, num_trackCh(event), 12);
        //All the den_ltof samples
        SetSample(den_ltof, den_l1_ltof(event), 1, den_ltof.reco_inn, event.trTrackBase->Rigidity[FIT][INN], nucleus, num_l1(event),2);
        SetSample(den_ltof, den_tof_ltof(event), 3, den_ltof.reco_il1, event.trTrackBase->Rigidity[FIT][IL1],nucleus,  num_tof(event), 4);
        SetSample(den_ltof, den_trig(event) && (event.evSummary->PhysBPatt & 65), 5, den_ltof.reco_il1, 
                  event.trTrackBase->Rigidity[FIT][IL1], nucleus,  event.evSummary->IsPhysicsTrigger(), 6);
        if (den_track(event)) {
            den_ltof.reco_beta = FillHistos(*sample_tr,charge,Amass,weight,chain,event);
            den_ltof.passed |= 1<<7;
            den_ltof.species = nucleus;
            if(num_track(event)) {
              den_ltof.passed |= 1<<8;
              den_ltof.species = nucleus;
            }
        }
        SetSample(den_ltof, den_l1Unb_ltof(event), 9, den_ltof.reco_inn, event.trTrackBase->Rigidity[FIT][INN], nucleus, num_l1Unb(event), 10);
        SetSample(den_ltof, den_trackCh_ltof(event), 11, den_ltof.reco_il1, event.trTrackBase->Rigidity[FIT][IL1], nucleus, num_trackCh(event), 12);
        //Al the den_ltof_l9Unb_05 samples
        SetSample(den_ltof_l9Unb_05, den_l1_ltof_l9_05(event), 1, den_ltof_l9Unb_05.reco_inn, event.trTrackBase->Rigidity[FIT][INN], nucleus, num_l1(event),2);
        SetSample(den_ltof_l9Unb_05, den_tof_ltof_l9_05(event), 3, den_ltof_l9Unb_05.reco_il1, event.trTrackBase->Rigidity[FIT][IL1], nucleus, num_tof(event), 4);
        SetSample(den_ltof_l9Unb_05, den_trig(event) && (event.evSummary->PhysBPatt & 65), 5, den_ltof_l9Unb_05.reco_il1, 
                  event.trTrackBase->Rigidity[FIT][IL1], nucleus, event.evSummary->IsPhysicsTrigger(), 6);
        if (den_track(event)) {
            den_ltof_l9Unb_05.reco_beta = FillHistos(*sample_tr,charge,Amass,weight,chain,event);
            den_ltof_l9Unb_05.passed |= 1<<7;
            den_ltof_l9Unb_05.species = nucleus;
            if(num_track(event)) {
              den_ltof_l9Unb_05.passed |= 1<<8;
              den_ltof_l9Unb_05.species = nucleus;
            }
        }
        SetSample(den_ltof_l9Unb_05, den_l1Unb_ltof_l9_05(event), 9, den_ltof_l9Unb_05.reco_inn, event.trTrackBase->Rigidity[FIT][INN],nucleus, num_l1Unb(event), 10);
        SetSample(den_ltof_l9Unb_05, den_trackCh_ltof_l9_05(event), 11, den_ltof_l9Unb_05.reco_il1, event.trTrackBase->Rigidity[FIT][IL1], nucleus,num_trackCh(event), 12);
        //All the den_ltof_l9Unb_1 samples
        SetSample(den_ltof_l9Unb_1, den_l1_ltof_l9_1(event), 1, den_ltof_l9Unb_1.reco_inn, event.trTrackBase->Rigidity[FIT][INN], nucleus,num_l1(event),2);
        SetSample(den_ltof_l9Unb_1, den_tof_ltof_l9_1(event), 3, den_ltof_l9Unb_1.reco_il1, event.trTrackBase->Rigidity[FIT][IL1],nucleus, num_tof(event), 4);
        SetSample(den_ltof_l9Unb_1, den_trig(event) && (event.evSummary->PhysBPatt & 65), 5, den_ltof_l9Unb_1.reco_il1, 
                  event.trTrackBase->Rigidity[FIT][IL1], nucleus, event.evSummary->IsPhysicsTrigger(), 6);
        if (den_track(event)) {
            den_ltof_l9Unb_1.reco_beta = FillHistos(*sample_tr,charge,Amass,weight,chain,event);
            den_ltof_l9Unb_1.passed |= 1<<7;
            den_ltof_l9Unb_1.species = nucleus;
            if(num_track(event)) {
              den_ltof_l9Unb_1.passed |= 1<<8;
              den_ltof_l9Unb_1.species = nucleus;
            }
        }
        SetSample(den_ltof_l9Unb_1, den_l1Unb_ltof_l9_1(event), 9, den_ltof_l9Unb_1.reco_inn, event.trTrackBase->Rigidity[FIT][INN], nucleus,num_l1Unb(event), 10);
        SetSample(den_ltof_l9Unb_1, den_trackCh_ltof_l9_1(event), 11, den_ltof_l9Unb_1.reco_il1, event.trTrackBase->Rigidity[FIT][IL1],nucleus, num_trackCh(event), 12);
        /*if (den_l1(event)) {
            den.reco_inn_l1pick = event.trTrackBase->Rigidity[FIT][INN];
            den.passed |= 1<<1;
            if(num_l1(event)) den.passed |= 1<<2;
        }
        if (den_tof(event)) {
            den.reco_il1_tof = event.trTrackBase->Rigidity[FIT][IL1];
            den.passed |= 1<<3;
            if(num_tof(event)) den.passed |= 1<<4;
        }
        if (den_trig(event)) {
            den.reco_il1_trig = event.trTrackBase->Rigidity[FIT][IL1];
            if (event.evSummary->PhysBPatt & 65 ) {
                den.passed |= 1<<5;
				        if (event.evSummary->IsPhysicsTrigger() ) den.passed |= 1<<6;
			      }
        }
        if (den_track(event)) {
            den.reco_beta = FillHistos(*sample_tr,charge,Amass,weight,chain,event);
            den.passed |= 1<<7;
            if(num_track(event)) den.passed |= 1<<8;
        }
        if (den_l1Unb(event)) {
            den.reco_inn_l1Unb = event.trTrackBase->Rigidity[FIT][INN];
            den.passed |= 1<<9;
            if (num_l1Unb(event)) den.passed |= 1<<10;
        }
        if(den_trackCh(event)) {
            den.reco_il1_trackCh = event.trTrackBase->Rigidity[FIT][IL1];
            den.passed |= 1<<11;
            if(num_trackCh(event)) den.passed |= 1<<12;
        }*/

        if(den.passed!=0 ) tree1->Fill();
        if(den_ltof.passed!=0 ) tree2->Fill();
        if(den_ltof_l9Unb_05.passed!=0 ) tree3->Fill();
        if(den_ltof_l9Unb_1.passed!=0 ) tree4->Fill();
    }
    auto outFile = new TFile(out, "recreate");
    outFile->WriteTObject(tree1, "reco_gen1");
    outFile->WriteTObject(tree2, "reco_gen2");
    outFile->WriteTObject(tree3, "reco_gen3");
    outFile->WriteTObject(tree4, "reco_gen4");
    outFile->WriteTObject(mc_samp, "mc_samp");
    outFile->WriteTObject(mc_pass_gen, "mc_pass_gen");
    outFile->WriteTObject(mc_pass, "mc_pass");
    delete outFile;
}
void SetSample(Sample &den, bool condition, int bit, float &reco_var, float rigidity, int Z ) {
    if (condition) {
        reco_var = rigidity;
        den.passed |= (1 << bit);
        den.species = Z;
    }
}
void SetSample(Sample &den, bool condition, int bit, float &reco_var, float rigidity, int Z, bool subcondition, int subbit) {
    if (condition) {
        reco_var = rigidity;
        den.passed |= (1 << bit);
        den.species = Z;
        if (subcondition) {
            den.passed |= (1 << subbit);
            den.species = Z;
        }
    }
}
double FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain,  NAIA::Event &event) {
  bool IsMC = chain.IsMC();
  float IGRFCutoff;
  if (!IsMC) {
    const NAIA::RTIInfo &rti_info = chain.GetEventRTIInfo();
    IGRFCutoff = rti_info.MaxIGRFCutoff[2][1];
  } else
    IGRFCutoff = event.mcTruthBase->Primary.GetGenMomentum() / event.mcTruthBase->Primary.Z;
  auto utime = event.header->UTCTime;
  float Beta = event.tofBaseSt->Beta[NAIA::Tof::BetaType::BetaH];
  const float amu = 0.938272075;
  double mass = Amass * amu;
  double BetaRig = mass * Beta / sqrt(1 - Beta * Beta) / charge;

  float EnergyD = 0;
  if (NAIA::ContainsKeys(event.ecalBase->Energy, NAIA::Ecal::EnergyRecoType::EnergyD)) {
    EnergyD = event.ecalBase->Energy[NAIA::Ecal::EnergyRecoType::EnergyD];
  }
  if (BetaRig <= 5) {
    histo.Fill(utime, BetaRig, weight);
	  return BetaRig;
  } else if (IGRFCutoff > 5 && IGRFCutoff <= 30) {
    histo.Fill(utime, IGRFCutoff, weight);
	  return IGRFCutoff;
  } else if (EnergyD > 30 && Efficiency::TrTrackEffSel::IsInsideEcal()(event)) {
    histo.Fill(utime, EnergyD, weight);
	  return EnergyD;
  }
  return EnergyD;
}