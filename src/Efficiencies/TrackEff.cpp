#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

double FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain,  NAIA::Event &event);

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto sample_tr    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_tr      = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);    
    // this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString exename=argv[0],
        infilename=argv[2],
		    outname=argv[3],
		    out, ionPath=getIonPath(charge);
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
	out = StreamUtility::getOutputDir(charge,exename,outname);
  unsigned int utime, icut_track = 0;
	float reco_beta;
	double Amass = 1;
	double weight = 1;
  double cutoff;
	Amass = GetDataMass(charge); // Getting atomic mass number
  auto den_track = Den_track(charge);
  auto num_track = Num_track(charge);

	for(Event& event : chain) { ////LOOP
    	utime = event.header->UTCTime;
		if (!isMC) {
  		if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
				continue;
			if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
				continue;
      if (isRunBad(utime)==true) //skipping trigger study and wrong configuration
					continue;
			if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
				continue;
		}
      ////TRACK EFF////
      if(den_track(event)) {
      		reco_beta = FillHistos(*sample_tr,charge,Amass,weight,chain,event);
      		if(num_track(event)) {
				    reco_beta = FillHistos(*pass_tr,charge,Amass,weight,chain,event);
			}
      	}
    } ////LOOP
    auto outfile = new TFile(out, "recreate");
    outfile->WriteTObject(sample_tr, "sample_tr");
    outfile->WriteTObject(pass_tr, "pass_tr");
    outfile->Close();
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