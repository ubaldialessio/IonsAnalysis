#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto sample_l1Unb_den   		 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1Unb_den     		 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto sample_l1Unb_den_ltof      = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1Unb_den_ltof        = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);   
	auto sample_l1Unb_den_ltof_l9_05= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1Unb_den_ltof_l9_05  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);   
	auto sample_l1Unb_den_ltof_l9_1 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_l1Unb_den_ltof_l9_1   = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ); 
	auto qL1Unb_den		            = new TH1D("", ";Q_L1Unb;",550,0,28); 
	auto qL1Unb_den_ltof		   	= new TH1D("", ";Q_L1Unb;",550,0,28);
	auto qL1Unb_den_ltof_l9_05	    = new TH1D("", ";Q_L1Unb;",550,0,28); 
	auto qL1Unb_den_ltof_l9_1	    = new TH1D("", ";Q_L1Unb;",550,0,28); 
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
		out = StreamUtility::getOutputDir(charge,exename,outname);
		unsigned int utime, icut_l1u = 0;
		float reco_inn;
		double cutoff, weight = 1, qUnbL1;			
		std::vector<std::string> labels_l1UnbEff;
					labels_l1UnbEff = {"Total",
							              "Trigger",
							              "Beta",
							              "Upper Tof Charge",
							              "Inner Hit Pattern",
							              "Inner Hit y>4",
							              "Inner Charge In Range",									  
							              "Inner #chi^{2}_{Y}<10",
							              "L1Fid",
							              "InnerFid"};						
			auto ncuts_l1UnbEff = labels_l1UnbEff.size();
			auto counters_l1UnbEff = new TH1D("counters_l1UnbEff", "", ncuts_l1UnbEff, -0.5, ncuts_l1UnbEff - 0.5);
			auto fill_counters_l1UnbEff = [&counters_l1UnbEff, &icut_l1u, &labels_l1UnbEff](Event &event) {
				 counters_l1UnbEff->GetXaxis()->SetBinLabel(icut_l1u + 1, labels_l1UnbEff.at(icut_l1u).c_str());
				 counters_l1UnbEff->Fill(icut_l1u);
				 icut_l1u++;
			};
auto num_l1Unb =
	track::IsHitPresent(EL1) &&
    track::UnbExtHitChargeStatus(EL1) &&
    track::UnbExtHitChargeInRange(EL1,static_cast<float>(charge),CRT,"UNB");
auto den_l1Unb = 
	ns::Trigger::HasPhysicsTrigger() && 
	ns::Tof::BetaInRange(0.4,inf,BTH) && 
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.2,static_cast<float>(charge)+0.2,UTC) && 
	ns::InnerTracker::HitPattern() && 
	ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
	ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT) && 
	ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN) && 
	ns::Track::L1FiducialVolume(FIT, INN) && 
	ns::Track::InnerFiducialVolume(FIT, INN);
auto den_l1Unb_ltof = 
	den_l1Unb &= 
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC);
auto den_l1Unb_ltof_l9_05 = 
	den_l1Unb_ltof &=
	track::ExtChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+100.,EL9,CRT);
auto den_l1Unb_ltof_l9_1 = 
	den_l1Unb_ltof &=
	track::ExtChargeInRange(static_cast<float>(charge)-0.,static_cast<float>(charge)+100.,EL9,CRT);

/////////////////// LOOP ///////////////////////////////
		for(Event& event : chain) {
			icut_l1u = 0;
			utime = event.header->UTCTime;
			if (!isMC) {
				if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
					continue;
				if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
					continue;
				if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
					continue;
				cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
			}
			////L1 UNB EFF/////
			if(den_l1Unb_ltof_l9_05(event)) {
      			reco_inn = event.trTrackBase->Rigidity[FIT][INN];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qUnbL1 = event.extHitBase->Charge(EL1,CRT);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
      					sample_l1Unb_den_ltof_l9_05->Fill(utime, reco_inn);
						qL1Unb_den_ltof_l9_05->Fill(qUnbL1);
      		        	if(num_l1Unb(event)) pass_l1Unb_den_ltof_l9_05->Fill(utime, reco_inn);
      		   	 	}
      		   	} else {
      		   		sample_l1Unb_den_ltof_l9_05->Fill(utime, reco_inn);
					qL1Unb_den_ltof_l9_05->Fill(qUnbL1);
      		   		if (num_l1Unb(event)) { 
						pass_l1Unb_den_ltof_l9_05->Fill(utime, reco_inn); 
					}
      		   	}
      		}
      		if(den_l1Unb(event)) {
      			reco_inn = event.trTrackBase->Rigidity[FIT][INN];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qUnbL1 = event.extHitBase->Charge(EL1,CRT);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
      					sample_l1Unb_den->Fill(utime, reco_inn);
						qL1Unb_den->Fill(qUnbL1);
      		        	if(num_l1Unb(event)) pass_l1Unb_den->Fill(utime, reco_inn);
      		   	 	}
      		   	} else {
      		   		sample_l1Unb_den->Fill(utime, reco_inn);
					qL1Unb_den->Fill(qUnbL1);
      		   		if (num_l1Unb(event)) pass_l1Unb_den->Fill(utime, reco_inn); 
      		   	}
      		}
			if(den_l1Unb_ltof(event)) {
      			reco_inn = event.trTrackBase->Rigidity[FIT][INN];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qUnbL1 = event.extHitBase->Charge(EL1,CRT);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
      					sample_l1Unb_den_ltof->Fill(utime, reco_inn);
						qL1Unb_den_ltof->Fill(qUnbL1);
      		        	if(num_l1Unb(event)) pass_l1Unb_den_ltof->Fill(utime, reco_inn);
      		   	 	}
      		   	} else {
      		   		sample_l1Unb_den_ltof->Fill(utime, reco_inn);
					qL1Unb_den_ltof->Fill(qUnbL1);
      		   		if (num_l1Unb(event)) { 
						pass_l1Unb_den_ltof->Fill(utime, reco_inn); 
					}
      		   	}
      		}
			if(den_l1Unb_ltof_l9_1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[FIT][INN];
      		    int rbin = hist_rig_highZ->FindBin(reco_inn);
      		    double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
				qUnbL1 = event.extHitBase->Charge(EL1,CRT);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
      					sample_l1Unb_den_ltof_l9_1->Fill(utime, reco_inn);
						qL1Unb_den_ltof_l9_1->Fill(qUnbL1);
      		        	if(num_l1Unb(event)) pass_l1Unb_den_ltof_l9_1->Fill(utime, reco_inn);
      		   	 	}
      		   	} else {
      		   		sample_l1Unb_den_ltof_l9_1->Fill(utime, reco_inn);
					qL1Unb_den_ltof_l9_1->Fill(qUnbL1);
      		   		if (num_l1Unb(event)) { 
						pass_l1Unb_den_ltof_l9_1->Fill(utime, reco_inn); 
					}
      		   	}
      		}
		} //event loop	
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(counters_l1UnbEff, "counters_l1UnbEff");
	outfile->WriteTObject(sample_l1Unb_den, "sample_l1Unb_den");
	outfile->WriteTObject(pass_l1Unb_den, "pass_l1Unb_den");
	outfile->WriteTObject(sample_l1Unb_den_ltof, "sample_l1Unb_den_ltof");
	outfile->WriteTObject(pass_l1Unb_den_ltof, "pass_l1Unb_den_ltof");
	outfile->WriteTObject(sample_l1Unb_den_ltof_l9_05, "sample_l1Unb_den_ltof_l9_05");
	outfile->WriteTObject(pass_l1Unb_den_ltof_l9_05, "pass_l1Unb_den_ltof_l9_05");
	outfile->WriteTObject(sample_l1Unb_den_ltof_l9_1, "sample_l1Unb_den_ltof_l9_1");
	outfile->WriteTObject(pass_l1Unb_den_ltof_l9_1, "pass_l1Unb_den_ltof_l9_1");
	outfile->WriteTObject(qL1Unb_den,"qL1Unb_den");
	outfile->WriteTObject(qL1Unb_den_ltof,"qL1Unb_den_ltof");
	outfile->WriteTObject(qL1Unb_den_ltof_l9_05,"qL1Unb_den_ltof_l9_05");
	outfile->WriteTObject(qL1Unb_den_ltof_l9_1,"qL1Unb_den_ltof_l9_1");
	outfile->Close();
	} //valid input
} //main