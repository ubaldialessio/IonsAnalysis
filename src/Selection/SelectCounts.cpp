#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto rigidity     = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);

	// this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> -sec_track <y/n> \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString exename = argv[0],
		    infilename=argv[2],
		    outname=argv[3],
		    out, ionPath=getIonPath(charge),
			sec_track = argv[4];

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
		if (sec_track=="y") outname = "sec_track_"+outname;
		out = StreamUtility::getOutputDir(charge,exename,outname);
		unsigned int utime;
		float reco_il1;
		double  cutoff;
		unsigned int icut_sel = 0;
		std::vector<std::string> labels_Sel;
					labels_Sel = {"Total",
				              "Trigger",
				              "Beta>0.4",
				              "UTofChargeInRange",
				              "InnerHitPattern",
				              "N Hits y > 4",
							  "Charge RMS < 0.55",
				              "Inn tracker charge",
							  "L1NormRes<10",
							  "L1 charge asymmetry < 0.2",
				              "L1Charge In Range",
							  "L1 good status",
				              "IL1 #chi^{2}_{Y}",
				              "Inner #chi^{2}_{Y}",
				              "Hit on L1",
							  "L1 fid volume",
							  "Inner fid volume"};				
			auto ncuts_Sel = labels_Sel.size();
			auto counters_Sel = new TH1D("counters_Sel", "", ncuts_Sel, -0.5, ncuts_Sel - 0.5);
			auto fill_counters_Sel = [&counters_Sel, &icut_sel, &labels_Sel](Event &event) {
				 counters_Sel->GetXaxis()->SetBinLabel(icut_sel + 1, labels_Sel.at(icut_sel).c_str());
				 counters_Sel->Fill(icut_sel);
				 icut_sel++;
			};	
auto mysel = Mysel(charge);
auto hasGoodSecondTrack = GoodSecTrack();
auto secTrackOnDiagonal = MySel::IsSecondTrackOnDiagonal(FIT,IL1,0.3);

    for(Event& event : chain) {
        icut_sel = 0;
		utime = event.header->UTCTime;
		if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
			continue;
		if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
			continue;
		if (isRunBad(utime)==true) //skipping trigger study and wrong configuration
			continue;
		if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
			continue;
		cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
		////MY SELECTION////
		if (mysel(event)) {
			if (sec_track=="y") {
				if (hasGoodSecondTrack(event)) {
					if (secTrackOnDiagonal(event)) {
						continue;
					}
				}
			}
			reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
			int rbin = hist_rig_highZ->FindBin(reco_il1);
			double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			if(rlowedge > 1.2*cutoff) rigidity->Fill(utime, reco_il1);
		}
    }
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(counters_Sel, "counters_Sel");
	outfile->WriteTObject(rigidity, "rigidity");
	outfile->Close();
}