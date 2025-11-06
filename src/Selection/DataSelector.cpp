#include "utils.h"
#include "binning.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

    std::vector<unsigned int> UnderStudyCharge = {14,15,16};

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

    unsigned int utime, oldtime=0;
    float reco_il1, reco_inn, beta, lt, reco_beta;
	double gen, cutoff;
    std::array<unsigned int, 24> daq;
	std::array<double, 8> daq_event_size;
    unsigned short nAcc_counters = 0, nAcc_clusters=0;

    // Prepare ChargeAnalysisManager
    ChargeAnalysisManager manager(
        UnderStudyCharge,
        nTbins, Tbins,
        nRbins_HighZ, Rbins_HighZ, isMC, exename
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
            if (!NSL::Selections::DefaultRTISelection(chain.GetEventRTIInfo()))
                continue;
            cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[1][1];
        }

        // Loop charges under study
        for (auto charge : UnderStudyCharge) {
            auto& ctx = manager.GetContext(charge);
            double Amass = ctx.Amass;
            double weight = 1.0;

            if (ctx.mysel(event)) {
                reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
                int rbin = hist_rig_highZ->FindBin(reco_il1);
                double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                if (!isMC) {
                    if (rlowedge > 1.2 * cutoff)
                        manager.Fill(charge, "IL1", "rigidity", utime, reco_il1);
                } else {
                    manager.Fill(charge, "IL1", "rigidity", utime, reco_il1);
                }
            }
            if (ctx.den_tof(event)) {
                reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_tof_den", utime, reco_il1);
      		    	if(ctx.num_tof(event)) 
                        manager.Fill(charge, "IL1", "pass_tof_den", utime, reco_il1);
      		    }
            }
            if (ctx.den_tof_ltof(event)) {
                reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_tof_den_ltof", utime, reco_il1);
      		    	if(ctx.num_tof(event)) 
                        manager.Fill(charge, "IL1", "pass_tof_den_ltof", utime, reco_il1);
      		    }
            }
            if (ctx.den_tof_ltof_l9_05(event)) {
                reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_tof_den_ltof_l9_05", utime, reco_il1);
      		    	if(ctx.num_tof(event)) 
                        manager.Fill(charge, "IL1", "pass_tof_den_ltof_l9_05", utime, reco_il1);
      		    }
            }
            if (ctx.den_tof_ltof_l9_05(event)) {
                reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_tof_den_ltof_l9_1", utime, reco_il1);
      		    	if(ctx.num_tof(event)) 
                        manager.Fill(charge, "IL1", "pass_tof_den_ltof_l9_1", utime, reco_il1);
      		    }
            }

            if (ctx.den_l1(event)) {
                reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_l1_den", utime, reco_inn);
      		    	if(ctx.num_l1(event)) 
                        manager.Fill(charge, "IL1", "pass_l1_den", utime, reco_inn);
      		    }
            }
            if (ctx.den_l1_ltof(event)) {
                reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_l1_den_ltof", utime, reco_inn);
      		    	if(ctx.num_l1(event)) 
                        manager.Fill(charge, "IL1", "pass_l1_den_ltof", utime, reco_inn);
      		    }
            }
            if (ctx.den_l1_ltof_l9_05(event)) {
                reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_l1_den_ltof_l9_05", utime, reco_inn);
      		    	if(ctx.num_l1(event)) 
                        manager.Fill(charge, "IL1", "pass_l1_den_ltof_l9_05", utime, reco_inn);
      		    }
            }
            if (ctx.den_l1_ltof_l9_1(event)) {
                reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_l1_den_ltof_l9_1", utime, reco_inn);
      		    	if(ctx.num_l1(event)) 
                        manager.Fill(charge, "IL1", "pass_l1_den_ltof_l9_1", utime, reco_inn);
      		    }
            }

            if (ctx.den_l1Unb(event)) {
                reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_l1Unb_den", utime, reco_inn);
      		    	if(ctx.num_l1Unb(event))
                        manager.Fill(charge, "IL1", "pass_l1Unb_den", utime, reco_inn);
      		    }
            }
            if (ctx.den_l1Unb_ltof(event)) {
                reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                   
                    manager.Fill(charge, "IL1", "sample_l1Unb_den_ltof", utime, reco_inn);
      		    	if(ctx.num_l1Unb(event))
                        manager.Fill(charge, "IL1", "pass_l1Unb_den_ltof", utime, reco_inn);
      		    }
            }
            if (ctx.den_l1Unb_ltof_l9_05(event)) {
                reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_l1Unb_den_ltof_l9_05", utime, reco_inn);
      		    	if(ctx.num_l1Unb(event))
                        manager.Fill(charge, "IL1", "pass_l1Unb_den_ltof_l9_05", utime, reco_inn);
      		    }
            }
            if (ctx.den_l1Unb_ltof_l9_1(event)) {
                reco_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
                if(rlowedge > 1.2*cutoff) { 
                    
                    manager.Fill(charge, "IL1", "sample_l1Unb_den_ltof_l9_1", utime, reco_inn);
      		    	if(ctx.num_l1Unb(event))
                        manager.Fill(charge, "IL1", "pass_l1Unb_den_ltof_l9_1", utime, reco_inn);
      		    }
            }

            if(ctx.den_track(event)) {
                weight=1;
      			reco_beta = FillHistos(*manager.GetHist(charge,"IL1","sample_tr"),
                                        charge, Amass, weight, chain, event);
      			if(ctx.num_track(event)) 
                    reco_beta = FillHistos(*manager.GetHist(charge,"IL1","pass_tr"),
                             charge, Amass, weight, chain, event);
      		}
            
            if (ctx.den_trackCh(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
				
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.Fill(charge, "IL1", "sample_trch_den", utime, reco_il1);
      		        	if(ctx.num_trackCh(event))
                            manager.Fill(charge, "IL1", "pass_trch_den", utime, reco_il1);
      		   	 	}
      		   	} else {
                    manager.Fill(charge, "IL1", "sample_trch_den", utime, reco_il1);
      		    	if(ctx.num_trackCh(event))
                        manager.Fill(charge, "IL1", "pass_trch_den", utime, reco_il1);
      		   	}
            }
            if (ctx.den_trackCh_ltof(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
				
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.Fill(charge, "IL1", "sample_trch_den_ltof", utime, reco_il1);
      		        	if(ctx.num_trackCh(event))
                            manager.Fill(charge, "IL1", "pass_trch_den_ltof", utime, reco_il1);
      		   	 	}
      		   	} else {
                    manager.Fill(charge, "IL1", "sample_trch_den_ltof", utime, reco_il1);
      		        if(ctx.num_trackCh(event))
                        manager.Fill(charge, "IL1", "pass_trch_den_ltof", utime, reco_il1);
      		   	}
            }
            if (ctx.den_trackCh_ltof_l9_05(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
				
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.Fill(charge, "IL1", "sample_trch_den_ltof_l9_05", utime, reco_il1);
      		        	if(ctx.num_trackCh(event))
                            manager.Fill(charge, "IL1", "pass_trch_den_ltof_l9_05", utime, reco_il1);
      		   	 	}
      		   	} else {
                    manager.Fill(charge, "IL1", "sample_trch_den_ltof_l9_05", utime, reco_il1);
      		    	if(ctx.num_trackCh(event))
                        manager.Fill(charge, "IL1", "pass_trch_den_ltof_l9_05", utime, reco_il1);               
      		   	}
            }
            if (ctx.den_trackCh_ltof_l9_1(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                
				
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
                        manager.Fill(charge, "IL1", "sample_trch_den_ltof_l9_1", utime, reco_il1);
      		        	if(ctx.num_trackCh(event))
                            manager.Fill(charge, "IL1", "pass_trch_den_ltof_l9_1", utime, reco_il1);
      		   	 	}
      		   	} else {
                    manager.Fill(charge, "IL1", "sample_trch_den_ltof_l9_1", utime, reco_il1);
      		    	if(ctx.num_trackCh(event))
                        manager.Fill(charge, "IL1", "pass_trch_den_ltof_l9_1", utime, reco_il1);
      		   	}
            }            

			////TRIGGER EFF WITH DIFFERENT PERIODS////	
            nAcc_clusters = event.evSummary->NAntiCluster;
            nAcc_counters = event.evSummary->NAcc;
			if (ctx.den_trig(event)) {
				weight = 1.0;
				reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
				int rbin = hist_rig_highZ->FindBin(reco_il1);
				double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                event.evSummary->RestorePhysBPatt(true);
				weight = event.evSummary->TriggerWeight();	
				if (!isMC) {
					if(rlowedge > 1.2*cutoff) {
						auto PhysBPatt = event.evSummary->PhysBPatt;
                        if (PhysBPatt & 65 && !event.evSummary->IsPhysicsTrigger()) 
							manager.Fill(charge, "IL1", "unbiased_trig", utime, reco_il1, weight);
                        if (event.evSummary->IsPhysicsTrigger()) {
                            manager.Fill(charge, "IL1", "phys_trig", utime, reco_il1, weight);
                            if (nAcc_clusters<5) manager.Fill(charge, "IL1", "phys_trig_acc4", utime, reco_il1, weight);
                            if (nAcc_clusters<7) manager.Fill(charge, "IL1", "phys_trig_acc7", utime, reco_il1, weight);
                            //Periods with nACC<5 are < t3
                            if (utime<t3) {
                                manager.Fill(charge, "IL1", "nPhys_1", utime, reco_il1);
                                if (nAcc_clusters<5) manager.Fill(charge, "IL1", "acc5_clusters_pre", utime, reco_il1); 
                                if (nAcc_counters<5) manager.Fill(charge, "IL1", "acc5_counters_pre", utime, reco_il1); 
                            }
                            //Periods with nAcc<=8 are t>=t3
                            if (utime >= t3) {
                                manager.Fill(charge, "IL1", "nPhys_2", utime, reco_il1);
								//for BZ trigger (bit3==1, bit 2,4,5,6==0), select NACC<5:
								//counters
                                if ( ((PhysBPatt & 0x4L)!=0) && ((PhysBPatt & 0x3AL)==0) && nAcc_counters<5 )
									manager.Fill(charge, "IL1", "acc5_counters", utime, reco_il1); 
								if ( ((PhysBPatt&0x4L) == 0) || ((PhysBPatt & 0x3AL) != 0))
									manager.Fill(charge, "IL1", "acc5_counters", utime, reco_il1); 
								//cluster
                                if ( ((PhysBPatt & 0x4L)!=0) && ((PhysBPatt & 0x3AL)==0) && nAcc_clusters<5 )
									manager.Fill(charge, "IL1", "acc5_clusters", utime, reco_il1); 
								if ( ((PhysBPatt&0x4L) == 0) || ((PhysBPatt & 0x3AL) != 0))
									manager.Fill(charge, "IL1", "acc5_clusters", utime, reco_il1); 

                                if (nAcc_counters<8) manager.Fill(charge, "IL1", "acc8_counters", utime, reco_il1); 							
                                if (nAcc_clusters<8) manager.Fill(charge, "IL1", "acc8_clusters", utime, reco_il1); 
                            }
                        }
				    }
			    }
     		}  
            //DAQ EFF WITH EMULATION
                for (int i=0; i<7; i++) daq_event_size[i] = 0;
                double daq_sdr_1st_half_size = daq[ 4]+daq[ 5]+daq[ 6]+daq[ 7]; // SDR0, SDR1
                double daq_sdr_2nd_half_size = daq[18]+daq[19]+daq[20]+daq[21]; // SDR2, SDR3
                double daq_trigger_size = daq[14]+daq[15];
                double daq_tracker_1st_half_size = daq[ 0]+daq[ 1]+daq[ 3]+daq[ 9]; // T0, ..., T3 
                double daq_tracker_2nd_half_size = daq[16]+daq[17]+daq[22]+daq[23]; // T4, ..., T7
                double daq_tracker_T3 = daq[1];
                double daq_trd_U0_size = daq[8];
                double daq_trd_U1_size = daq[2];
                double daq_rich_size = daq[10]+daq[11];
                double daq_ecal_size = daq[12]+daq[13];
                // 4 JINJ configuration (since Feb. 2020)
                daq_event_size[3] = daq_tracker_2nd_half_size; // JINJ-0
                daq_event_size[4] = daq_rich_size+daq_sdr_2nd_half_size+daq_trd_U0_size; // JINJ-1 
                daq_event_size[5] = daq_ecal_size+daq_sdr_1st_half_size+daq_trd_U1_size+daq_trigger_size+daq_tracker_T3; // JINJ-2
                daq_event_size[6] = daq_tracker_1st_half_size-daq_tracker_T3; // JINJ-3;
                // 2 JINJ configuration (since Feb. 2016)
                daq_event_size[1] = daq_event_size[3]+daq_event_size[4]; // JINJ-1 
                daq_event_size[2] = daq_event_size[5]+daq_event_size[6]; // JINJ-2
                // 1 JINJ (earlier data)
                daq_event_size[0] = daq_event_size[1]+daq_event_size[2];
			if (ctx.mysel(event)) {	
				reco_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
				int rbin = hist_rig_highZ->FindBin(reco_il1);
				double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				if (!isMC) {
					if(rlowedge > 1.2*cutoff) {
						//Period 1 -- 1JINJ
						// I exclude t1 <= t < t2. Period 1, with 1 JINJ, ends on t2 exlcuding t1 <= t < t2
						if (utime>=t1 && utime<=t2) continue;
						if (utime<t2) {
                            manager.Fill(charge, "IL1", "jinj1_DAQevent_length", utime, daq_event_size[0]);
							if (daq_event_size[0] <= 24500)
                                manager.Fill(charge, "IL1", "jinj1_belowBuffer", utime, reco_il1);
							if (daq_event_size[0] > 24500)
								manager.Fill(charge, "IL1", "jinj1_aboveBuffer", utime, reco_il1);
						}
						//Period 2 -- 2JINJ
						// Period 2, with 2 or 4 JINJ, starts on t2, excluding photon trigger	
						if (utime >=t2) {
                            manager.Fill(charge, "IL1", "jinj2_DAQevent_length", utime, daq_event_size[0]);
							if (daq_event_size[0] <= 24500)
                                manager.Fill(charge, "IL1", "jinj2_belowBuffer", utime, reco_il1);
							if (daq_event_size[0] > 24500)
                                manager.Fill(charge, "IL1", "jinj2_aboveBuffer", utime, reco_il1);
						}   
					}
				}
			}
        }
    }
	for (auto charge : UnderStudyCharge) {
		out = StreamUtility::getOutputDir(charge,exename,outname);
    	manager.Write(ionPaths, out.Data(), isMC, charge);
	}

    return 0;
}

