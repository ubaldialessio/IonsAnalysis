#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {
	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root>  \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString exename = argv[0], infilename=argv[2],
		      outname=argv[3],
		      out, ionName=getIonName(charge);

out = StreamUtility::getOutputDir(charge,exename,outname);
std::cout << out << std::endl;
auto map       = new TH2D("", ";Cutoff (GV); RecoRigidity (GV)", nRbins_HighZ-1, Rbins_HighZ, nRbins_HighZ-1, Rbins_HighZ);
auto sample_tr = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
auto pass_tr   = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);

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
double utime,cutoff,Amass,weight,reco_il1,reco_beta=0;
Amass = GetDataMass(charge); // Getting atomic mass number
NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;
auto mysel =
        ns::Trigger::HasPhysicsTrigger() && 
        ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),UTC) && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) && 
		    ns::InnerTracker::ChargeRMSLessThan(0.55,CRT) &&
        ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT) && 
        InnerTracker::L1NormResidualLessThan(10.0,FIT) &&
		    ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT) && 
        ns::Common::TrackerLayer::ChargeInRange(1,static_cast<float>(charge),CRT) &&
		    ns::TrackerLayer::ChargeStatus(1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) && 
        ns::Track::HitCut(1) && 
        ns::Track::L1FiducialVolume(FIT, IL1) && 
		    ns::Track::InnerFiducialVolume(FIT,IL1); 
    for(Event& event : chain) {
        utime = event.header->UTCTime;
        if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
          continue;
        if (utime < 1305756000 || utime >1620252000) //getting data only from 19th May 2011 to 6th May 2021
          continue;
        if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
          continue;
        if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
          continue;
        if (mysel(event)) {
            reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
            cutoff = chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
            map->Fill(cutoff,reco_il1);
        }
    }
  auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(map, "map");
	outfile->Close();
}
