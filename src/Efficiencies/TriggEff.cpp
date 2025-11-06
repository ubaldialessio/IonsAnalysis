#include "utils.h"
#include "binning.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

    std::vector<unsigned int> UnderStudyCharge = {14};

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

    unsigned int utime;
    float reco_il1, lt;
	double cutoff;
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
 
 			////TRIGGER EFF WITH DIFFERENT PERIODS////
			///Trying to select nAcc < 5 only for BigZ trigger
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
     
        }
    }
    // WRITE
	for (auto charge : UnderStudyCharge) {
		out = StreamUtility::getOutputDir(charge,exename,outname);
    	manager.Write(ionPaths, out.Data(), isMC, charge);
	}

    return 0;
}

