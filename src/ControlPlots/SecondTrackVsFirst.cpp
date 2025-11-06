#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

	const int nRbins_Track = 54;
	const double Rbins_Track[nRbins_Track] = {0.1,0.3,0.7,1.,1.5,2.15,2.40,2.67,2.97,3.29,3.64,
	4.02,4.43,4.88,5.37,5.90, 6.47,7.09,7.76,8.48,9.26, 10.1,11.0,12.0,13.0,
	14.1,15.3,16.6,18.0,19.5,21.1,22.8,24.7,26.7,28.8,33.5,38.9,45.1,
	52.2,60.3,69.7,80.5,93.0,108.,125.,147.,175.,211,259.,330.,441.,
	660.,1200.,3000.};
	auto hist_rig_Track= new TH1D("","", nRbins_Track-1,Rbins_Track);


	TH1::SetDefaultSumw2();
	auto firstTrack_Zinn	= new TH1D("",";Z_{inn};Counts",550,0,28);
	auto secondTrack_Zinn   = new TH1D("",";Z_{inn};Counts",550,0,28);
	auto firstVsSecond_Zinn = new TH2D("",";Z1_{inn};Z2_{inn}",550,0,28,550,0,28);
    auto firstTrackVsSecond                = new TH2D("", ";R1_{Inn} (GV); R2_{INN} (GV)", nRbins_Track-1, Rbins_Track, nRbins_Track-1, Rbins_Track);
    auto firstTrackVsSecond_2ndGood        = new TH2D("", ";R1_{Inn} (GV); R2_{INN} (GV)", nRbins_Track-1, Rbins_Track, nRbins_Track-1, Rbins_Track);
    auto firstTrackVsSecond_2ndGood_ZGT1_8 = new TH2D("", ";R1_{Inn} (GV); R2_{INN} (GV)", nRbins_Track-1, Rbins_Track, nRbins_Track-1, Rbins_Track);

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
		bool isMC = chain.IsMC();
		out = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/ControlVariables/SecondTrackVsFirst/"+getIonPath(charge)+"/"+outname;
		unsigned int utime;
		float first_track_reco_il1, sec_track_reco_il1, z1_inn, z2_inn;
		double  cutoff;
		auto myselStoP = MyselStoP(charge);
		auto goodSecTrack = GoodSecTrack();
    for(Event& event : chain) {
		utime = event.header->UTCTime;
		if (!isMC) {
			if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
				continue;
			if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
				continue;
			cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
		}
		////MY SELECTION////
		if (myselStoP(event)) {
			first_track_reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
			int rbin = hist_rig_Track->FindBin(first_track_reco_il1);
			double rlowedge = hist_rig_Track->GetBinLowEdge(rbin);
			if (!isMC) {
				if(rlowedge > 1.2*cutoff) {
					sec_track_reco_il1 = event.secondTrTrackBase->Rigidity[FIT][IL1];
					firstTrackVsSecond->Fill(first_track_reco_il1,sec_track_reco_il1);
					if (goodSecTrack(event)) {
						firstTrackVsSecond_2ndGood->Fill(first_track_reco_il1,sec_track_reco_il1);
						z1_inn = event.trTrackBase->InnerCharge[CRT];
						z2_inn = event.secondTrTrackBase->InnerCharge[CRT];
						firstTrack_Zinn->Fill(z1_inn);
						secondTrack_Zinn->Fill(z2_inn);
						firstVsSecond_Zinn->Fill(z1_inn,z2_inn);
						if (event.secondTrTrackBase->InnerCharge[CRT]>1.8) {
							firstTrackVsSecond_2ndGood_ZGT1_8->Fill(first_track_reco_il1,sec_track_reco_il1);
						}
					}
				}
			} else {
				sec_track_reco_il1 = event.secondTrTrackBase->Rigidity[FIT][IL1];
				firstTrackVsSecond->Fill(first_track_reco_il1,sec_track_reco_il1);
				if (goodSecTrack(event)) {
					firstTrackVsSecond_2ndGood->Fill(first_track_reco_il1,sec_track_reco_il1);
					z1_inn = event.trTrackBase->InnerCharge[CRT];
					z2_inn = event.secondTrTrackBase->InnerCharge[CRT];
					firstTrack_Zinn->Fill(z1_inn);
					secondTrack_Zinn->Fill(z2_inn);
					firstVsSecond_Zinn->Fill(z1_inn,z2_inn);
					if (event.secondTrTrackBase->InnerCharge[CRT]>1.8) {
						firstTrackVsSecond_2ndGood_ZGT1_8->Fill(first_track_reco_il1,sec_track_reco_il1);
					}
				}
			}
		}
    }
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(firstTrack_Zinn,"firstTrack_Zinn");
	outfile->WriteTObject(secondTrack_Zinn,"secondTrack_Zinn");
	outfile->WriteTObject(firstVsSecond_Zinn,"firstVsSecond_Zinn");
	outfile->WriteTObject(firstTrackVsSecond, "firstTrackVsSecond");
	outfile->WriteTObject(firstTrackVsSecond_2ndGood, "firstTrackVsSecond_2ndGood");
	outfile->WriteTObject(firstTrackVsSecond_2ndGood_ZGT1_8, "firstTrackVsSecond_2ndGood_ZGT1_8");
	outfile->Close();
}