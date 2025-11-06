#include "definition.h"
#include "binning.h"
#include "utils.h"

int main(int argc, char *argv[]) {
	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root>  \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString infilename=argv[2],
		    outname=argv[3],
		    out, ionPath=getIonName(charge);
            out="../Fragmentation/BelowL1/"+ionPath+"/l1_"+outname;
    ////////List definitions////////
    TList *l1TemplateList   = new TList();
	////////The binning template for the L1 distribution (to be fitted) need to be the same of the event selection (73 bins, high Z)
													  
	auto rigBinTemplate   = new TH1D("RigTemplateLightIons", "R (GV)", nRbins_HighZ - 1, Rbins_HighZ);
	
	for (Int_t ibin = 0; ibin < rigBinTemplate->GetNbinsX(); ibin++) {
	    l1TemplateList->Add(new TH1D(Form("L1Template_%03i", ibin),
	                               Form("%5.3f < R (GV) < %5.3f;Q_{L1X};Counts", rigBinTemplate->GetBinLowEdge(ibin + 1),
	                                    rigBinTemplate->GetBinLowEdge(ibin + 2)),
	                               600, 0, 20));
	}
    unsigned int utime;
    ////////L1 template selection////////
	TH1::SetDefaultSumw2();
    auto l1temp =
        BuildTemplatesSel::HitPattern() &&
        BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge), UTC) &&       
        BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge), LTC) &&
        ns::TrackerLayer::ChargeStatus(1) &&
        BuildTemplatesSel::InnerTrackerChargeInRange(static_cast<float>(charge), CRT) &&
        BuildTemplatesSel::InnerTrackerChargeRMSLessThan(0.07, CRT);
        if (charge < 9)
        l1temp &= BuildTemplatesSel::InnerTrackerNtrackLessThan(3) &&
        BuildTemplatesSel::HasGoodSecondTrTrack(FIT,INN, 0.2) &&
        BuildTemplatesSel::HasGoodTRDHits(charge);
	////////this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;
    ////////Building chain////////
	NAIA::NAIAChain chain;
	if(infilename.Contains(".root") /*&& filesystem::exists(infilename.Data())*/ ){
		chain.Add(infilename.Data());
	}else if (infilename.Contains(".txt") /*&& filesystem::exists(infilename.Data())*/ ){
		ifstream infilelist(infilename.Data());
		TString bufname;
		while(infilelist >> bufname) 
		    if(bufname.Contains(".root") /*&& filesystem::exists(bufname.Data())*/ )
		    chain.Add(bufname.Data());
	}
	chain.SetupBranches();
    /////////////////// LOOP ///////////////////////////////
    for(Event& event : chain) {
	    utime = event.header->UTCTime;
		if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
			continue;
		if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
			continue;
		if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
			continue;
		if (l1temp(event)) {
			auto Rigidity_IL1 = event.trTrackBase->Rigidity[FIT][IL1];
			int bin = rigBinTemplate->FindBin(Rigidity_IL1);
			if (bin < 1 || bin > rigBinTemplate->GetNbinsX())
				continue;
			auto L1charge = event.trTrackBase->LayerChargeXY[0][CRT];
			static_cast<TH1D *>(l1TemplateList->At(bin - 1))->Fill(L1charge);
		}
    }
	////////WRITE////////
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(l1TemplateList, Form("L1TemplateList_%i", charge), "Overwrite");
	outfile->Close();
}