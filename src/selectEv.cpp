#include "definition.h"
#include "selection.h"

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


		
		for(Event &event : chain) {
			//starting the selection
			if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )
				continue;
		
    		utime = event.header->UTCTime;
      		cutoff = chain.GetEventRTIInfo().MaxIGRFCutoff[0][1];
      		lt     = chain.GetEventRTIInfo().LivetimeFraction;
      		if(lt < 0.05
      			|| chain.GetEventRTIInfo().Zenith > 25
      			|| chain.GetEventRTIInfo().MeanAlignDiffExtLayer[0][1] > 35
      			|| chain.GetEventRTIInfo().IsInSAA()) continue;
      		if(utime != oldtime) {
        		for (int rbin = 1; rbin <=hist_rig->GetNbinsX(); rbin++) {
          		auto rlowedge = hist_rig->GetBinLowEdge(rbin);
          		auto rcenter  = hist_rig->GetBinCenter(rbin);
          		if(rlowedge > safety_factor[rbin]*cutoff) lvt_25->Fill(utime, rcenter, lt);
        		} //lvt fill
      		} // if utime != oldtime
		} //event loop
	}
