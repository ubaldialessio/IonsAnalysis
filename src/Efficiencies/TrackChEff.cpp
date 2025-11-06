#include "utils.h"
#include "binning.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

    std::vector<unsigned int> UnderStudyCharge = {8,9,10,14,15,16,17,18,19,20,
                                                  21,22,23,24,25,26};

    TH1::SetDefaultSumw2();

    // Time markers
    unsigned int t0 = 1305853512, t1 = 1385487767, t2 = 1447346927, t3 = 1456503197,
                 t4 = 1582034309, t5 = 1582037855, t6 = 1620025528, t7 = 1635856717,
                 t8 = 1675341999;
    unsigned int t11= 1447346927 , t22 = 1456503197;	

    if (argc < 3) {
        printf("Usage: \n");
        printf("%s <path/to/input.root/.txt> <output.root> \n", argv[0]);
        return 1;
    }

    TString executable = argv[0],
            infilename = argv[1],
            outname = argv[2],
			out;
    TString exename = executable(executable.Last('/') + 1, executable.Length());
    // Setup chain
    NAIA::NAIAChain chain;
    if (infilename.Contains(".root")) {
        chain.Add(infilename.Data());
    } else if (infilename.Contains(".txt")) {
        std::ifstream infilelist(infilename.Data());
        TString bufname;
        while (infilelist >> bufname)
            if (bufname.Contains(".root"))
                chain.Add(bufname.Data());
    }
    chain.SetupBranches();
    bool isMC = chain.IsMC();

    unsigned int utime;
    float reco_il1, lt;
	double gen, cutoff;

    // Prepare ChargeAnalysisManager
    ChargeAnalysisManager manager(
        UnderStudyCharge,
        nTbins, Tbins,
        nRbins_HighZ, Rbins_HighZ,isMC,exename
    );

    // Map Z → ionPath
    std::map<unsigned int, std::string> ionPaths;
	for (auto charge : UnderStudyCharge) {
		ionPaths[charge] = std::string(getIonPath(charge).Data());
	}

    /////////////////// LOOP ///////////////////////////////
    for (NAIA::Event &event : chain) {

		utime = event.header->UTCTime;

        // Skip periods if real data
        if (!isMC) {
            if (utime >= t6 && utime < t7) continue; // photon-polarization
            if (isRunBad(utime)) continue;           // bad runs
            if (!NAIA::MatchAnyBit(event.header->Mask(),
                 NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof))
                continue;
            if (! NSL::Selections::DefaultRTISelection(chain.GetEventRTIInfo()) ||
                chain.GetEventRTIInfo().IsInSAA())
                continue;
            cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[1][1];
        }

        // Loop charges under study
        for (auto charge : UnderStudyCharge) {
            auto& ctx = manager.GetContext(charge);
            double Amass = ctx.Amass;
            double weight = 1.0;
            
            if (ctx.den_trackCh(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
				
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.Fill(charge, "IL1", "sample_trch_den", utime, reco_il1, weight);
      		        	if(ctx.num_trackCh(event))
                            manager.Fill(charge, "IL1", "pass_trch_den", utime, reco_il1, weight);
      		   	 	}
      		   	} else {
                    manager.Fill(charge, "IL1", "sample_trch_den", utime, reco_il1, weight);
      		    	if(ctx.num_trackCh(event))
                        manager.Fill(charge, "IL1", "pass_trch_den", utime, reco_il1, weight);
      		   	}
            }
            if (ctx.den_trackCh_ltof(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
				
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.Fill(charge, "IL1", "sample_trch_den_ltof", utime, reco_il1, weight);
      		        	if(ctx.num_trackCh(event))
                            manager.Fill(charge, "IL1", "pass_trch_den_ltof", utime, reco_il1, weight);
      		   	 	}
      		   	} else {
                    manager.Fill(charge, "IL1", "sample_trch_den_ltof", utime, reco_il1, weight);
      		        if(ctx.num_trackCh(event))
                        manager.Fill(charge, "IL1", "pass_trch_den_ltof", utime, reco_il1, weight);
      		   	}
            }
            if (ctx.den_trackCh_ltof_l9_05(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
				
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.Fill(charge, "IL1", "sample_trch_den_ltof_l9_05", utime, reco_il1, weight);
      		        	if(ctx.num_trackCh(event))
                            manager.Fill(charge, "IL1", "pass_trch_den_ltof_l9_05", utime, reco_il1, weight);
      		   	 	}
      		   	} else {
                    manager.Fill(charge, "IL1", "sample_trch_den_ltof_l9_05", utime, reco_il1, weight);
      		    	if(ctx.num_trackCh(event))
                        manager.Fill(charge, "IL1", "pass_trch_den_ltof_l9_05", utime, reco_il1, weight);               
      		   	}
            }
            if (ctx.den_trackCh_ltof_l9_1(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
				
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.Fill(charge, "IL1", "sample_trch_den_ltof_l9_1", utime, reco_il1, weight);
      		        	if(ctx.num_trackCh(event))
                            manager.Fill(charge, "IL1", "pass_trch_den_ltof_l9_1", utime, reco_il1, weight);
      		   	 	}
      		   	} else {
                    manager.Fill(charge, "IL1", "sample_trch_den_ltof_l9_1", utime, reco_il1, weight);
      		    	if(ctx.num_trackCh(event))
                        manager.Fill(charge, "IL1", "pass_trch_den_ltof_l9_1", utime, reco_il1, weight);
      		   	}
            }            
        }
    }
    // WRITE
	for (auto charge : UnderStudyCharge) {
		out = StreamUtility::getOutputDir(charge,exename,outname);
    	manager.Write(ionPaths, out.Data(), isMC, charge);
	}

    return 0;
}

