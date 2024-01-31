#include "selection.h"

auto lvt_25 	  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
auto rigidity     = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
auto sample_tof   = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
auto pass_tof     = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);

int main(int argc, char *argv[]) {
	if (argc < 2) {
	    printf("Usage: \n");
	    printf("%s <input root file> <output> \n", argv[0]);
	    return 1;
	}

	//check if is a valid input
    bool validInput = true;
	//check(); 
	TString out = argv[2];

	//process
	if (validInput) {
		NAIA::NAIAChain chain;
		chain.Add(argv[1]);
		chain.SetupBranches();

		unsigned int utime, oldtime;
		float reco_il1;


/////////////////// LOOP ///////////////////////////////
		for(Event &event : chain) {
			// skipping photon-polarization period
			auto UTime = event.header->UTCTime;
			if (UTime>=1620025528 && UTime<1635856717) 
				continue;


			//RTI check
			if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )
				continue;
				

			utime  		= event.header->UTCTime;
			auto cutoff = chain.GetEventRTIInfo().MaxIGRFCutoff[1][1];
			auto lt     = chain.GetEventRTIInfo().LivetimeFraction;


			////MY SELECTION////
			if (mysel(event)) {
				reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
				int rbin = hist_rig->FindBin(reco_il1);
				double rlowedge = hist_rig->GetBinLowEdge(rbin);
				if(rlowedge > 1.2*cutoff) rigidity->Fill(utime, reco_il1);
			}


      		if(utime != oldtime) {
        		for (int rbin = 1; rbin <=hist_rig->GetNbinsX(); rbin++) {
          			auto rlowedge = hist_rig->GetBinLowEdge(rbin);
          			auto rcenter  = hist_rig->GetBinCenter(rbin);
          			if(rlowedge > 1.2*cutoff) lvt_25->Fill(utime, rcenter, lt);
        		} //lvt fill
      		} // if utime != oldtime


			////TOF EFF////
      		if(den_tof(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      		    int rbin = hist_rig->FindBin(reco_il1);
      		    double rlowedge = hist_rig->GetBinLowEdge(rbin);
      		    if(rlowedge > 1.2*cutoff){ 
      				sample_tof->Fill(utime, reco_il1);
      		        if(num_tof(event)) pass_tof->Fill(utime, reco_il1);
      		    }
      		}
    

        	
		} //event loop
	}
	
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(rigidity, "rigidity");
	outfile->WriteTObject(lvt_25, "lvt_25");
	outfile->WriteTObject(sample_tof, "sample_tof");
	outfile->WriteTObject(pass_tof, "pass_tof");
}
