#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto jinj1_DAQevent_length      = new TH1D("","",2500, 0, 50000);
	auto jinj2_DAQevent_length      = new TH1D("","",2500, 0, 50000);
	auto jinj1_aboveBuffer    		= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto jinj1_belowBuffer    		= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto jinj2_aboveBuffer    		= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto jinj2_belowBuffer    		= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);

    // period 1:  t<T0,    nACC<5 and 1 JINJ  (NoROOM error = JINJ-1<24500 B) 					   t < 1447346927
    // period 2a: T0<t<T1, nACC<5 and 2 JINJs (NoROOM error = JINJ-1<24500 B || JINJ-2<24500 B)    1447346927 <= t < 1456503197
    // period 2b: t>T1,    nACC<8 and 2 JINJs (NoROOM error = JINJ-1<24500 B || JINJ-2<24500 B)    t > = 1456503197
	// I exclude t1 <= t < t2. Period 1, with 1 JINJ, ends on t2 exlcuding t1 <= t < t2
	// Period 2, with 2 or 4 JINJ, starts on t2, excluding photon trigger	
	unsigned int t1= 1447346927 , t2 = 1456503197;

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
		bool isMC = chain.IsMC();
		chain.SetupBranches();
		out = StreamUtility::getOutputDir(charge,exename,outname);
		unsigned int utime;
		int JLength;
		float reco_il1;
		double  cutoff;
		unsigned short nAcc_counters = 0,nAcc_clusters=0;
		std::array<unsigned int, 24> daq;
		std::array<double, 8> daq_event_size;
		auto mysel = Mysel(charge);
		for(Event& event : chain) {
			daq = event.daq->JLength;
			utime = event.header->UTCTime;
			if (!isMC) {
				if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
					continue;
				if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
					continue;
				if (isRunBad(utime)==true) //skipping trigger study and wrong configuration
					continue;
				if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
					continue;
				cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
			}

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

			nAcc_clusters = event.evSummary->NAntiCluster;
			nAcc_counters = event.evSummary->NAcc;
			
			if (mysel(event)) {	
				reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
				int rbin = hist_rig_highZ->FindBin(reco_il1);
				double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				if (!isMC) {
					if(rlowedge > 1.2*cutoff) {
						//Period 1 -- 1JINJ
						// I exclude t1 <= t < t2. Period 1, with 1 JINJ, ends on t2 exlcuding t1 <= t < t2
						if (utime>=t1 && utime<=t2) continue;
						if (utime<t2) {
							jinj1_DAQevent_length->Fill(daq_event_size[0]);
							if (daq_event_size[0] <= 24500)
								jinj1_belowBuffer->Fill(utime,reco_il1);
							if (daq_event_size[0] > 24500)
								jinj1_aboveBuffer->Fill(utime,reco_il1);
						}
						//Period 2 -- 2JINJ
						// Period 2, with 2 or 4 JINJ, starts on t2, excluding photon trigger	
						if (utime >=t2) {
							jinj2_DAQevent_length->Fill(daq_event_size[0]);
							if (daq_event_size[0] <= 24500)
								jinj2_belowBuffer->Fill(utime,reco_il1);
							if (daq_event_size[0] > 24500)
								jinj2_aboveBuffer->Fill(utime,reco_il1);
						}
					}
				}
			}
    }
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(jinj1_DAQevent_length, "jinj1_DAQevent_length");
	outfile->WriteTObject(jinj2_DAQevent_length, "jinj2_DAQevent_length");
	outfile->WriteTObject(jinj1_aboveBuffer, "jinj1_aboveBuffer");
	outfile->WriteTObject(jinj1_belowBuffer, "jinj1_belowBuffer");
	outfile->WriteTObject(jinj2_aboveBuffer, "jinj2_aboveBuffer");
	outfile->WriteTObject(jinj2_belowBuffer, "jinj2_belowBuffer");
	outfile->Close();
}