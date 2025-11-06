#include "definition.h"
#include "binning.h"
#include "utils.h"

double FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain,  NAIA::Event &event);

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	unsigned int t0 = 1305853512, t1 = 1385487767, t2 = 1447346927, t3 = 1456503197, t4 = 1582034309,
				 t5 = 1582037855, t6 = 1620025528, t7 = 1635856717, t8 = 1675341999;
				 unsigned int t11= 1447346927 , t22 = 1456503197;

	auto acc5_counters_pre 		    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto acc5_counters 		    	= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto acc8_counters  			= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto acc5_clusters_pre			= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto acc5_clusters     		    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto acc8_clusters  			= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto nPhys_1					= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto nPhys_2					= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);

	auto jinj1_DAQevent_length      = new TH1D("","",2500, 0, 50000);
	auto jinj2_DAQevent_length      = new TH1D("","",2500, 0, 50000);
	auto jinj1_aboveBuffer    		= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto jinj1_belowBuffer    		= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto jinj2_aboveBuffer    		= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto jinj2_belowBuffer    		= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);

	auto qInn_R		   = new TH2D("", ";R(GV);Q_inn", nRbins_HighZ-1, Rbins_HighZ, 550,0,28); 
	auto qUnbL9_R	   = new TH2D("", ";R(GV);Q_L9;" ,nRbins_HighZ-1, Rbins_HighZ, 550,0,28);

	auto lvt_25	  	  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto rigidity     = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto sample_tof   = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_tof     = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto sample_l1    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1      = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto unbiased_trig= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto phys_trig    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto sample_tr    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_tr      = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto sample_l1Unb = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1Unb   = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto sample_trch  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_trch    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ); 

	auto noTriggers  = trig::HasTrigger(0);                	// 00000000
	auto unbiased    = trig::HasTrigger(1)||             	// 00000001 (only unbiased ToF)
	                   trig::HasTrigger(64)||            	// 01000000 (only unbiased ECAL)
	                   trig::HasTrigger(65);              	// 01000001 (unbiased ToF AND unbiased ECAL)
	auto physics     = ns::Trigger::HasPhysicsTrigger();    // 00111110 (Any non-prescaled trigger)	


	// this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString exename=argv[0],
			infilename=argv[2],
		    outname=argv[3],
		    out, ionPath=getIonPath(charge);
	//check if is a valid input
    bool validInput = true;
	//process
	if (validInput) {
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
		if (isMC) out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/mc/"+outname;
		else out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/dat/"+outname;
		unsigned int utime, oldtime;
		float reco_il1, reco_inn, beta, lt, reco_beta;
		double gen, cutoff,qinn, ql9unb;;
		double Amass = 1;
		double weight = 1;
		unsigned short nAcc_counters = 0, nAcc_clusters=0;
		Amass = GetDataMass(charge); // Getting atomic mass number
		std::array<unsigned int, 24> daq;
		std::array<double, 8> daq_event_size;
		auto mysel  = Mysel(charge);
		auto num_tof= Num_tof(charge);
		auto den_tof_ltof=Den_tof_ltof(charge);
		auto num_l1=Num_l1(charge);
		auto den_l1=Den_l1(charge);
		auto num_l1Unb=Num_l1Unb(charge);
		auto den_l1Unb=Den_l1Unb(charge);
		auto num_track=Num_track(charge);
		auto den_track=Den_track(charge);
		auto num_trackCh = Num_trackCh(charge);
		auto den_trackCh_ltof_l9unb_05 = Den_trackCh_ltof_l9_05(charge);
		auto den_trig = Den_trig(charge);

/////////////////// LOOP ///////////////////////////////
		for(Event& event : chain) {
			daq = event.daq->JLength;
			utime = event.header->UTCTime;
			if (!isMC) {
				if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
					continue;
				if (isRunBad(utime)==true) //skipping trigger study and wrong configuration
					continue;
				if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
					continue;
				if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
					continue;
				cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
				lt    		= chain.GetEventRTIInfo().LivetimeFraction;
				if(utime != oldtime) {
			    	for (int rbin = 1; rbin <=hist_rig_highZ->GetNbinsX(); rbin++) {
			        	auto rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			         	auto rcenter  = hist_rig_highZ->GetBinCenter(rbin);
			          	if(rlowedge > 1.2*cutoff) lvt_25->Fill(utime, rcenter, lt);
			       	} 
			    oldtime=utime;
			    }
			} else {
				gen = event.mcTruthBase->Primary.GetGenMomentum()/event.mcTruthBase->Primary.Z;
			}	


			////MY SELECTION////
			if (mysel(event)) {
				reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
				int rbin = hist_rig_highZ->FindBin(reco_il1);
				double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				if (!isMC) {
					if(rlowedge > 1.2*cutoff) rigidity->Fill(utime, reco_il1);
				}
			}
			////TOF EFF////
      		if(den_tof_ltof(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_il1);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) { 
      					sample_tof->Fill(utime, reco_il1);
      		        	if(num_tof(event)) pass_tof->Fill(utime, reco_il1);
      		    	}
      		    } else {
      		    	sample_tof->Fill(utime, reco_il1);
      		    	if (num_tof(event) ) {
						pass_tof->Fill(utime, reco_il1);
					}
      		    }
      		}
			////L1 EFF/////
      		if(den_l1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[FIT][IL1];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
      					sample_l1->Fill(utime, reco_inn);
      		        	if(num_l1(event)) pass_l1->Fill(utime, reco_inn);
      		   	 	}
      		   	} else {
      		   		sample_l1->Fill(utime, reco_inn);
      		   		if (num_l1(event)) { 
						pass_l1->Fill(utime, reco_inn); 
					}
      		   	}
      		}

			////L1 UNB EFF/////
      		if(den_l1Unb(event)) {
      			reco_inn = event.trTrackBase->Rigidity[FIT][INN];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
      					sample_l1Unb->Fill(utime, reco_inn);
      		        	if(num_l1Unb(event)) pass_l1Unb->Fill(utime, reco_inn);
      		   	 	}
      		   	} else {
      		   		sample_l1Unb->Fill(utime, reco_inn);
      		   		if (num_l1Unb(event)) { 
						pass_l1Unb->Fill(utime, reco_inn); 
					}
      		   	}
      		}
			////TRACK EFF////
      		if(den_track(event)) {
      			reco_beta = FillHistos(*sample_tr,charge,Amass,weight,chain,event);
      			if(num_track(event)) {
					reco_beta = FillHistos(*pass_tr,charge,Amass,weight,chain,event);
				}
      		}
			////TRACK CHARGE EFF////
			if(den_trackCh_ltof_l9unb_05(event)) {
					reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
					int rbin = hist_rig_highZ->FindBin(reco_il1);
					double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
					ql9unb = event.extHitBase->Charge(EL9, CRT);
					qinn   = event.trTrackBase->InnerCharge[CRT];
				if (!isMC) {
						if(rlowedge > 1.2*cutoff) { 
							qInn_R->Fill(qinn,reco_il1);
							qUnbL9_R->Fill(ql9unb,reco_il1);
							sample_trch->Fill(utime, reco_il1);
							if(num_trackCh(event)) 
								pass_trch->Fill(utime, reco_il1);
						}
					} else {
						qInn_R->Fill(qinn,reco_il1);
						qUnbL9_R->Fill(ql9unb,reco_il1);
						sample_trch->Fill(utime, reco_il1);
						if (num_trackCh(event) ) 
							pass_trch->Fill(utime, reco_il1);
				}
			}
			nAcc_clusters = event.evSummary->NAntiCluster;
			nAcc_counters = event.evSummary->NAcc;
			////TRIGGER EFF WITH DIFFERENT PERIODS////	
			if (den_trig(event) &&nAcc_clusters<8 && nAcc_counters<8) {
				reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
				int rbin = hist_rig_highZ->FindBin(reco_il1);
				double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				if (!isMC) {
					if(rlowedge > 1.2*cutoff) {
						//Periods with nACC<5 are < t3
						if (utime<t3) {
							if (physics(event)) {
								nPhys_1->Fill(utime,reco_il1);
							}
								if (nAcc_clusters<5) acc5_clusters_pre->Fill(utime,reco_il1);
								if (nAcc_counters<5) acc5_counters_pre->Fill(utime,reco_il1);
						}
						//Periods with nAcc<=8 are t>=t3
						if (utime >= t3) {
							if (physics(event)) {
								nPhys_2->Fill(utime,reco_il1);
							}
								if (nAcc_counters<5) acc5_counters->Fill(utime,reco_il1);
								if (nAcc_clusters<5) acc5_clusters->Fill(utime,reco_il1);
								if (nAcc_counters<8) acc8_counters->Fill(utime,reco_il1);
								if (nAcc_clusters<8) acc8_clusters->Fill(utime,reco_il1);
						}
					}
					event.evSummary->RestorePhysBPatt(true);
					weight = event.evSummary->TriggerWeight();	
					if ( rlowedge > 1.2*cutoff && event.evSummary->PhysBPatt & 65 )
						unbiased_trig->Fill(utime, reco_il1, weight);
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

		} //event loop	
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(rigidity, "rigidity");
	outfile->WriteTObject(lvt_25, "lvt_25");
	outfile->WriteTObject(sample_tof, "sample_tof_den_ltof");
	outfile->WriteTObject(pass_tof, "pass_tof_den_ltof");
	outfile->WriteTObject(sample_l1, "sample_l1_den");
	outfile->WriteTObject(pass_l1, "pass_l1_den");
	outfile->WriteTObject(sample_tr, "sample_tr");
	outfile->WriteTObject(pass_tr, "pass_tr");
	outfile->WriteTObject(unbiased_trig, "unbiased_trig");
	outfile->WriteTObject(phys_trig, "phys_trig");
	outfile->WriteTObject(sample_l1Unb, "sample_l1Unb_den");
	outfile->WriteTObject(pass_l1Unb, "pass_l1Unb_den");
    outfile->WriteTObject(sample_trch, "sample_trch_den_ltof_l9_05");
    outfile->WriteTObject(pass_trch, "pass_trch_den_ltof_l9_05");
	outfile->WriteTObject(acc5_clusters_pre, "acc5_clusters_pre");
	outfile->WriteTObject(acc5_counters_pre,"acc5_counters_pre");
	outfile->WriteTObject(acc8_clusters, "acc8_clusters");
	outfile->WriteTObject(acc8_counters, "acc8_counters");
	outfile->WriteTObject(acc5_clusters, "acc5_clusters");
	outfile->WriteTObject(acc5_counters, "acc5_counters");
	outfile->WriteTObject(nPhys_1,"nPhys_1");
	outfile->WriteTObject(nPhys_2,"nPhys_2");
	outfile->WriteTObject(jinj1_DAQevent_length, "jinj1_DAQevent_length");
	outfile->WriteTObject(jinj2_DAQevent_length, "jinj2_DAQevent_length");
	outfile->WriteTObject(jinj1_aboveBuffer, "jinj1_aboveBuffer");
	outfile->WriteTObject(jinj1_belowBuffer, "jinj1_belowBuffer");
	outfile->WriteTObject(jinj2_aboveBuffer, "jinj2_aboveBuffer");
	outfile->WriteTObject(jinj2_belowBuffer, "jinj2_belowBuffer");
	outfile->Close();
	} //valid input
} //main

double FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain,  NAIA::Event &event) {
  bool IsMC = chain.IsMC();
  float IGRFCutoff;
  if (!IsMC) {
    const NAIA::RTIInfo &rti_info = chain.GetEventRTIInfo();
    IGRFCutoff = rti_info.MaxIGRFCutoff[2][1];
  } else
    IGRFCutoff = event.mcTruthBase->Primary.GetGenMomentum() / event.mcTruthBase->Primary.Z;
  auto utime = event.header->UTCTime;
  float Beta = event.tofBaseSt->Beta[NAIA::Tof::BetaType::BetaH];
  const float amu = 0.938272075;
  double mass = Amass * amu;
  double BetaRig = mass * Beta / sqrt(1 - Beta * Beta) / charge;

  float EnergyD = 0;
  if (NAIA::ContainsKeys(event.ecalBase->Energy, NAIA::Ecal::EnergyRecoType::EnergyD)) {
    EnergyD = event.ecalBase->Energy[NAIA::Ecal::EnergyRecoType::EnergyD];
  }
  if (BetaRig <= 5) {
    histo.Fill(utime, BetaRig, weight);
	  return BetaRig;
  } else if (IGRFCutoff > 5 && IGRFCutoff <= 30) {
    histo.Fill(utime, IGRFCutoff, weight);
	  return IGRFCutoff;
  } else if (EnergyD > 30 && Efficiency::TrTrackEffSel::IsInsideEcal()(event)) {
    histo.Fill(utime, EnergyD, weight);
	  return EnergyD;
  }
  return EnergyD;
}