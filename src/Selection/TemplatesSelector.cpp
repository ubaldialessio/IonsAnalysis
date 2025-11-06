#include "utils.h"
#include "binning.h"

int main(int argc, char *argv[]) {

    std::vector<unsigned int> UnderStudyCharge = {13,14,15,16,17,18};

    TH1::SetDefaultSumw2();
    TH1::AddDirectory(kFALSE);
    gROOT->SetBatch(kTRUE);


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
    float reco_il1, L1charge, l1charge, l2charge, cutoff, lt;
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
		ionPaths[charge] = std::string(getIonName(charge).Data());
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
		    lt    		= chain.GetEventRTIInfo().LivetimeFraction;
        }

        // Loop charges under study
        for (auto charge : UnderStudyCharge) {
            auto& ctx = manager.GetContext(charge);
            double Amass = ctx.Amass;
            double weight = 1.0;

            auto Rigidity_IL1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
			int bin = hist_rig_highZ->FindBin(Rigidity_IL1);
			double rlowedge = hist_rig_highZ->GetBinLowEdge(bin);
            if (rlowedge < 1.2*cutoff) continue;
			if (bin < 1 || bin > hist_rig_highZ->GetNbinsX())	continue;
            
            if (ctx.myselNoL1(event)) {
                L1charge = event.trTrackBase->LayerChargeXY[0][NAIA::TrTrack::ChargeRecoType::YJ];
                manager.FillTemplateL1Distribution(charge, "IL1", bin-1, L1charge, 1);
            }
            if (ctx.l1_temp(event)) {
				l1charge = event.trTrackBase->LayerChargeXY[0][NAIA::TrTrack::ChargeRecoType::YJ];
                manager.FillTemplateL1(charge, "IL1", bin-1, l1charge, 1);
            }
            if (ctx.l2_temp(event)) {
				l2charge = event.trTrackBase->LayerChargeXY[1][NAIA::TrTrack::ChargeRecoType::YJ];
                manager.FillTemplateL2(charge, "IL1", bin-1, l2charge, 1);
            }   
        }
    }
    // WRITE
    manager.WriteTemplateDistributions(ionPaths, outname.Data(), isMC);

    return 0;
}

