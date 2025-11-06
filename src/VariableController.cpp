#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"
#include "definition.h"

int main(int argc, char *argv[]) {
	if (argc < 2) {
	    printf("Usage: \n");
	    printf("%s <input root file> <output root file>  \n", argv[0]);
	    return 1;
	}

	//check if is a valid input
    bool validInput = true;
	//check();
	TString outname = argv[2],infilename=argv[1];

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

		//create the variables and branch that I want to store and use
		float lt, zenith, q_utof, q_ltof, q_tot_tof, beta, betaSt, chargel1, chargeInner, rigInnerL1,
			  chi2_spatial, chi2_time, unbl1_charge, reco_il1;
		double cutoff;
		unsigned int utime;
		TTree *tree=new TTree();
		tree->Branch("q_utof", &q_utof);
		tree->Branch("q_ltof", &q_ltof);
		tree->Branch("q_tot_tof", &q_tot_tof);
		tree->Branch("chargel1", &chargel1);
		tree->Branch("inner_charge", &chargeInner);
		tree->Branch("beta", &beta);
		tree->Branch("betaSt", &betaSt);
		tree->Branch("chi2_spatial", &chi2_spatial);
		tree->Branch("chi2_time", &chi2_time);
		tree->Branch("unbl1_charge", &unbl1_charge);

		//loop event
		for (Event& event : chain) {
			if(!event.CheckMask(NAIA::Category::HasTrack) || !event.CheckMask(NAIA::Category::HasTof) || 
			   !event.CheckMask(NAIA::Category::HasTofStandalone) ) continue;
			utime = event.header->UTCTime;
			if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
				continue;
			if (utime < 1305756000 || utime >1620252000) //getting data only from 19th May 2011 to 6th May 2021
				continue;
			if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
				continue;
			cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
			lt    		= chain.GetEventRTIInfo().LivetimeFraction;

			reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
			int rbin = hist_rig_highZ->FindBin(reco_il1);
			double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			if(rlowedge > 1.2*cutoff) {
				q_utof 	   = event.tofBase->Charge[UTC];
				q_ltof     = event.tofBase->Charge[LTC];
				q_tot_tof  = event.tofBase->Charge[TOT];
				beta       = event.tofBase->Beta[BTH];
				betaSt	   = event.tofBaseSt->Beta[BTH];
				chargel1   = event.trTrackBase->LayerChargeXY[0][CRT]; //0= layer 1
				chargeInner= event.trTrackBase->InnerCharge[CRT];
				rigInnerL1 = event.trTrackBase->Rigidity[FIT][IL1];
				chi2_time  = event.tofPlusSt->Chi2Time;
				chi2_spatial=event.tofPlusSt->Chi2Coo;
				unbl1_charge= event.extHitBase->Charge(EL1,CRT);
				tree->Fill();
			}
		}
		auto outfile = new TFile("/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/ControlVariables/"+outname, "recreate");
		outfile->WriteTObject(tree, "variables");
		outfile->Close();
	}
}
