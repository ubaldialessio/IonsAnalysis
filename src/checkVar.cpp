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
	TString out = argv[2];

	//process
	if (validInput) {
		NAIA::NAIAChain chain;
		chain.Add(argv[1]);
		chain.SetupBranches();

		//create the variables and branch that I want to store and use
		float lt, zenith, utof_edep, ltof_edep, charge_tof, charge_ltof, beta, chargel1, chargeInner, rigInnerL1;
		unsigned int utime;
		unsigned short int cutmask, npart;
		TTree *tree=new TTree();
		tree->Branch("utime", &utime);
		tree->Branch("charge_tof", &charge_tof);
		tree->Branch("charge_ltof", &charge_ltof);
		tree->Branch("chargel1", &chargel1);
		tree->Branch("inner_charge", &chargeInner);
		tree->Branch("beta", &beta);
		tree->Branch("npart", &npart);
		tree->Branch("cutmask", &cutmask);
		tree->Branch("utof_edep", &utof_edep);
		tree->Branch("ltof_edep", &ltof_edep);

		//loop event
		for (Event& event : chain) {
			if(!event.CheckMask(NAIA::Category::HasTrack) || !event.CheckMask(NAIA::Category::HasTof) ) continue;
			//the analysis
			utime      = event.header->UTCTime;
		 	npart      = event.evSummary->NParticle;
			charge_tof = event.tofBase->Charge[UTC];
			charge_ltof = event.tofBase->Charge[LTC];
		    beta       = event.tofBase->Beta[BTH];
			utof_edep  = event.tofPlus->LayerEdep[0]+event.tofPlus->LayerEdep[1];
			ltof_edep  = event.tofPlus->LayerEdep[2]+event.tofPlus->LayerEdep[3];
			chargel1 = event.trTrackBase->LayerChargeXY[0][CRT]; //0= layer 1
			chargeInner = event.trTrackBase->InnerCharge[CRT];
			rigInnerL1 = event.trTrackBase->Rigidity[FIT][IL1];

			tree->Fill();
		}
		auto outfile = new TFile(out, "recreate");
		outfile->WriteTObject(tree, "variables");
		outfile->Close();
	}
}
