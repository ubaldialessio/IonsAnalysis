#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto lvt_25	  	  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);

	// this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);

	TString exename = argv[0],
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
	out = StreamUtility::getOutputDir(charge,exename,outname);
	unsigned int utime, oldtime;
	double cutoff;
	float lt;
    for(Event& event : chain) {
			utime = event.header->UTCTime;
			if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
				continue;
			if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
				continue;
			if (isRunBad(utime)==true) //skipping trigger study and wrong configuration
				continue;
			if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ))  //RTI check
				continue;
			cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[1][1];
			lt    		= chain.GetEventRTIInfo().LivetimeFraction;
			if(utime != oldtime) {
			    for (int rbin = 1; rbin <=hist_rig_highZ->GetNbinsX(); rbin++) {
			        auto rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			        auto rcenter  = hist_rig_highZ->GetBinCenter(rbin);
			        if(rlowedge > 1.2*cutoff) lvt_25->Fill(utime, rcenter, lt);
			    } 
			oldtime=utime;
			}
    }
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(lvt_25, "lvt_25");
	outfile->Close();
}