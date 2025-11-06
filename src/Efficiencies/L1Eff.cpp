#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto sample_l1_den   		 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1_den     		 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto sample_l1_den_ltof      = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1_den_ltof        = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);   
	auto sample_l1_den_ltof_l9_05= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1_den_ltof_l9_05  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);   
	auto sample_l1_den_ltof_l9_1 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1_den_ltof_l9_1   = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ); 
	auto qL1_den		         = new TH1D("", ";Q_L1;",550,0,28); 
	auto qL1_den_ltof		   	 = new TH1D("", ";Q_L1;",550,0,28);
	auto qL1_den_ltof_l9_05	   = new TH1D("", ";Q_L1;",550,0,28); 
	auto qL1_den_ltof_l9_1	   = new TH1D("", ";Q_L1;",550,0,28); 
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
	//check if is a valid input
    bool validInput = true;
	//process
	if (validInput) {
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
		unsigned int utime, icut_l1 = 0;
		float reco_inn;
		double cutoff, weight = 1, qL1;			
		auto num_l1=Num_l1(charge);
		auto den_l1=Den_l1(charge);

/////////////////// LOOP ///////////////////////////////
		for(Event& event : chain) {
			icut_l1 = 0;
			utime = event.header->UTCTime;
			if (!isMC) {
				if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
					continue;
				if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
				if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
					continue;
				cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
			}
			////L1 EFF/////
      		if(den_l1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[FIT][IL9];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qL1 = event.trTrackBase->LayerChargeXY[1][CRT];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
      					sample_l1_den->Fill(utime, reco_inn);
						qL1_den->Fill(qL1);
      		        	if(num_l1(event)) pass_l1_den->Fill(utime, reco_inn);
      		   	 	}
      		   	} else {
      		   		sample_l1_den->Fill(utime, reco_inn);
					qL1_den->Fill(qL1);
      		   		if (num_l1(event)) { 
						pass_l1_den->Fill(utime, reco_inn); 
					}
      		   	}
      		}
		} //event loop	
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(sample_l1_den, "sample_l1_den");
	outfile->WriteTObject(pass_l1_den, "pass_l1_den");
	outfile->WriteTObject(qL1_den,"qL1_den");
	outfile->WriteTObject(qL1_den_ltof,"qL1_den_ltof");
	outfile->WriteTObject(qL1_den_ltof_l9_05,"qL1_den_ltof_l9_05");
	outfile->WriteTObject(qL1_den_ltof_l9_1,"qL1_den_ltof_l9_1");
	outfile->Close();
	} //valid input
} //main