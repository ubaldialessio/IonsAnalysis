#include "binning.h"
#include "utils.h"

int main(int argc, char *argv[]) {

    std::vector<unsigned int> UnderStudyCharge = {8,  9,  10, 
                                                  14, 15, 16,
                                                  17, 18, 19, 20, 21, 22, 23, 26};

    TH1::SetDefaultSumw2();

    // Time markers
    unsigned int t0 = 1305853512, t1 = 1385487767, t2 = 1447346927, t3 = 1456503197,
                 t4 = 1582034309, t5 = 1582037855, t6 = 1620025528, t7 = 1635856717,
                 t8 = 1675341999;
    unsigned int t11= 1447346927 , t22 = 1456503197;

    auto noTriggers  = Efficiency::TriggerEff::HasTrigger(0);            // 00000000
	auto unbiased    = Efficiency::TriggerEff::HasTrigger(1)||           // 00000001 (only unbiased ToF)
	                   Efficiency::TriggerEff::HasTrigger(64)||          // 01000000 (only unbiased ECAL)
	                   Efficiency::TriggerEff::HasTrigger(65);           // 01000001 (unbiased ToF AND unbiased ECAL)
	auto physics     = NSL::Selections::Trigger::HasPhysicsTrigger();    // 00111110 (Any non-prescaled trigger)	

    if (argc < 3) {
        printf("Usage: \n");
        printf("%s <path/to/input.root/.txt> <output.root> \n", argv[0]);
        return 1;
    }

    TString executable = argv[0],
            infilename = argv[1],
            outname = argv[2];
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

    unsigned int utime, oldtime=0;
    float reco_il1, reco_inn, beta, lt, reco_beta;
	double gen, cutoff,qinn, ql9unb, qL1, qUnbL1, qtof;
    std::array<unsigned int, 24> daq;
	std::array<double, 8> daq_event_size;
    unsigned short nAcc_counters = 0, nAcc_clusters=0;

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

		daq = event.daq->JLength;
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
        }
        
        cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[1][1];

        // Loop charges under study
        for (auto charge : UnderStudyCharge) {
            auto& ctx = manager.GetContext(charge);
            double Amass = ctx.Amass;
            double weight = 1.0;

            if (ctx.den_l1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qL1 = event.trTrackBase->LayerChargeXY[1][NAIA::TrTrack::ChargeRecoType::YJ];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qL1_den", qL1);
                        manager.Fill(charge, "IL1", "sample_l1_den", utime, reco_inn, weight);
      		        	if(ctx.num_l1(event)) 
                            manager.Fill(charge, "IL1", "pass_l1_den", utime, reco_inn, weight);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qL1_den", qL1);
      		   	}
            }
            if (ctx.den_l1_ltof(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qL1 = event.trTrackBase->LayerChargeXY[1][NAIA::TrTrack::ChargeRecoType::YJ];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qL1_den_ltof", qL1);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qL1_den_ltof", qL1);
      		   	}
            }
            if (ctx.den_l1_ltof_l9_05(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qL1 = event.trTrackBase->LayerChargeXY[1][NAIA::TrTrack::ChargeRecoType::YJ];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qL1_den_ltof_l9_05", qL1);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qL1_den_ltof_l9_05", qL1);
      		   	}
            }
            if (ctx.den_l1_ltof_l9_1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qL1 = event.trTrackBase->LayerChargeXY[1][NAIA::TrTrack::ChargeRecoType::YJ];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qL1_den_ltof_l9_1", qL1);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qL1_den_ltof_l9_1", qL1);
      		   	}
            }

            if (ctx.den_l1Unb(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qUnbL1 = event.extHitBase->Charge(NAIA::UnbExtHitBaseData::ExtHit::L1,NAIA::TrTrack::ChargeRecoType::YJ);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qL1Unb_den", qUnbL1);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qL1Unb_den", qUnbL1);
      		   	}
            }
            if (ctx.den_l1Unb_ltof(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qUnbL1 = event.extHitBase->Charge(NAIA::UnbExtHitBaseData::ExtHit::L1,NAIA::TrTrack::ChargeRecoType::YJ);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qL1Unb_den_ltof", qUnbL1);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qL1Unb_den_ltof", qUnbL1);
      		   	}
            }
            if (ctx.den_l1Unb_ltof_l9_05(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qUnbL1 = event.extHitBase->Charge(NAIA::UnbExtHitBaseData::ExtHit::L1,NAIA::TrTrack::ChargeRecoType::YJ);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qL1Unb_den_ltof_l9_05", qUnbL1);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qL1Unb_den_ltof_l9_05", qUnbL1);
      		   	}
            }
            if (ctx.den_l1Unb_ltof_l9_1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qUnbL1 = event.extHitBase->Charge(NAIA::UnbExtHitBaseData::ExtHit::L1,NAIA::TrTrack::ChargeRecoType::YJ);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qL1Unb_den_ltof_l9_1", qUnbL1);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qL1Unb_den_ltof_l9_1", qUnbL1);
      		   	}
            }

            if (ctx.den_tof(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qtof = event.tofBase->Charge[NAIA::Tof::ChargeType::Upper];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qtof_den", qtof);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qtof_den", qtof);
      		   	}
            }
            if (ctx.den_tof_ltof(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qtof = event.tofBase->Charge[NAIA::Tof::ChargeType::Upper];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qtof_den_ltof", qtof);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qtof_den_ltof", qtof);
      		   	}
            }
            if (ctx.den_tof_ltof_l9_05(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qtof = event.tofBase->Charge[NAIA::Tof::ChargeType::Upper];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qtof_den_ltof_l9_05", qtof);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qtof_den_ltof_l9_05", qtof);
      		   	}
            }
            if (ctx.den_tof_ltof_l9_1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qtof = event.tofBase->Charge[NAIA::Tof::ChargeType::Upper];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qtof_den_ltof_l9_1", qtof);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qtof_den_ltof_l9_1", qtof);
      		   	}
            }

            if (ctx.den_trackCh(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                ql9unb = event.extHitBase->Charge(NAIA::UnbExtHitBaseData::ExtHit::L9, NAIA::TrTrack::ChargeRecoType::YJ);
				qinn   = event.trTrackBase->InnerCharge[NAIA::TrTrack::ChargeRecoType::YJ];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qInn_R_den", qinn);
                        manager.FillDistributions(charge, "IL1", "qUnbL9_R_den", ql9unb);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qInn_R_den", qinn);
                    manager.FillDistributions(charge, "IL1", "qUnbL9_R_den", ql9unb);
      		   	}
            }
            if (ctx.den_trackCh_ltof(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                ql9unb = event.extHitBase->Charge(NAIA::UnbExtHitBaseData::ExtHit::L9, NAIA::TrTrack::ChargeRecoType::YJ);
				qinn   = event.trTrackBase->InnerCharge[NAIA::TrTrack::ChargeRecoType::YJ];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qInn_R_den_ltof", qinn);
                        manager.FillDistributions(charge, "IL1", "qUnbL9_R_den_ltof", ql9unb);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qInn_R_den_ltof", qinn);
                    manager.FillDistributions(charge, "IL1", "qUnbL9_R_den_ltof", ql9unb);
      		   	}
            }
            if (ctx.den_trackCh_ltof_l9_05(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                ql9unb = event.extHitBase->Charge(NAIA::UnbExtHitBaseData::ExtHit::L9, NAIA::TrTrack::ChargeRecoType::YJ);
				qinn   = event.trTrackBase->InnerCharge[NAIA::TrTrack::ChargeRecoType::YJ];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qInn_R_den_ltof_l9_05", qinn);
                        manager.FillDistributions(charge, "IL1", "qUnbL9_R_den_ltof_l9_05", ql9unb);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qInn_R_den_ltof_l9_05", qinn);
                    manager.FillDistributions(charge, "IL1", "qUnbL9_R_den_ltof_l9_05", ql9unb);
      		   	}
            }
            if (ctx.den_trackCh_ltof_l9_1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                ql9unb = event.extHitBase->Charge(NAIA::UnbExtHitBaseData::ExtHit::L9, NAIA::TrTrack::ChargeRecoType::YJ);
				qinn   = event.trTrackBase->InnerCharge[NAIA::TrTrack::ChargeRecoType::YJ];
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.FillDistributions(charge, "IL1", "qInn_R_den_ltof_l9_1", qinn);
                        manager.FillDistributions(charge, "IL1", "qUnbL9_R_den_ltof_l9_1", ql9unb);
      		   	 	}
      		   	} else {
                    manager.FillDistributions(charge, "IL1", "qInn_R_den_ltof_l9_1", qinn);
                    manager.FillDistributions(charge, "IL1", "qUnbL9_R_den_ltof_l9_1", ql9unb);
      		   	}
            }            
        }
    }
    // WRITE
    manager.WriteDistributions(ionPaths, outname.Data(), isMC);
    return 0;
}
