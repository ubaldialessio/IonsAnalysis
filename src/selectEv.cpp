#include "definition.h"
#include "binning.h"
#include "utils.h"

void FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain, NAIA::Event &event);
double GetDataMass(unsigned int charge);

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	
	auto tree		  = new TTree();
	auto mc_samp      = (TH1D*)hist_log->Clone("mc_samp");
	auto mc_pass_gen  = (TH1D*)hist_log->Clone("mc_pass_gen");
	auto mc_pass      = (TH1D*)hist_log->Clone("mc_pass");

	auto lvt_25	  	  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto rigidity     = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto sample_tof   = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto pass_tof     = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto sample_l1    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto pass_l1      = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto sample_trig  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto phys_trig    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto unbiased_trig= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto sample_tr    = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);
	auto pass_tr      = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins-1, Rbins);

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
	TString infilename=argv[2],
		   outname=argv[3],
		   out, ionPath=getIonPath(charge);
	//check if is a valid input
    bool validInput = true;
	//process
	if (validInput) {
		NAIA::NAIAChain chain;
		if(infilename.Contains(".root") && filesystem::exists(infilename.Data()) ){
		    chain.Add(infilename.Data());
		}else if (infilename.Contains(".txt") && filesystem::exists(infilename.Data()) ){
		    ifstream infilelist(infilename.Data());
		    TString bufname;
		    while(infilelist >> bufname) 
		      if(bufname.Contains(".root") && filesystem::exists(bufname.Data()))
		        chain.Add(bufname.Data());
		}
		chain.SetupBranches();
		bool isMC = chain.IsMC();
		if (isMC) out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/mc/"+outname;
		else out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/dat/"+outname;
		unsigned int utime, oldtime;
		float reco_il1, reco_inn, beta, lt;
		double gen, cutoff;
		double Amass = 1;
		double weight = 1;
		Amass = GetDataMass(charge); // Getting atomic mass number
		if(isMC){
		    tree->Branch("utime", &utime);
		    tree->Branch("reco_inn", &reco_inn);
		    tree->Branch("reco_il1", &reco_il1);
		    //tree->Branch("reco_beta", &reco_beta);
		    tree->Branch("gen", &gen);
		    tree->Branch("beta", &beta);
		    //tree->Branch("passed", &passed);
		    ///////////////////////////////////////
		    auto mcchain = chain.GetFileInfoTree();
		    auto mcinfo = new NAIA::MCFileInfo();
		    mcchain->SetBranchAddress("MCFileInfo", &mcinfo);
		    for(unsigned long long i = 0; i < mcchain->GetEntries(); ++i){
		      mcchain->GetEntry(i);
		      double rminDC = mcinfo->GetRMin();
		      double rmaxDC = mcinfo->GetRMax();
		      unsigned int ngen = mcinfo->EventNo.second - mcinfo->EventNo.first + 1;
		      for(int ibin = 1; ibin <= mc_samp->GetNbinsX(); ++ibin){
		        double rmin = mc_samp->GetBinLowEdge(ibin), rmax = mc_samp->GetBinLowEdge(ibin+1);
		        if (rmin < rminDC) rmin = rminDC;
		        if (rmax > rmaxDC) rmax = rmaxDC;
		        if (rmin > rmaxDC || rmax < rminDC) continue;
		        mc_samp->AddBinContent(ibin, ngen * TMath::Log(rmax / rmin) / TMath::Log(rmaxDC / rminDC));
		      }
		    }
		}
		unsigned int icut_sel = 0, icut_track = 0, icut_tof = 0, icut_trig = 0, icut_l1 = 0;
		std::vector<std::string> labels_Sel;
					labels_Sel = {"Total",
				              "Trigger",
				              "BetaInRange",
				              "UTofChargeInRange",
				              "TofGoodPath",
				              "InnerHitPattern",
				              "NGoodCl",
				              "InnerChReso<0.55",
				              "InnerChInRange",
				              "Inner #chi^{2}_{Y}",
				              "InnerL1 #chi^{2}_{Y}",
				              "L1HitCut",
				              "L1Charge In Range",
				              "L1NormRes"};			
			auto ncuts_Sel = labels_Sel.size();
			auto counters_Sel = new TH1D("counters_Sel", "", ncuts_Sel, -0.5, ncuts_Sel - 0.5);
			auto fill_counters_Sel = [&counters_Sel, &icut_sel, &labels_Sel](Event &event) {
				 counters_Sel->GetXaxis()->SetBinLabel(icut_sel + 1, labels_Sel.at(icut_sel).c_str());
				 counters_Sel->Fill(icut_sel);
				 icut_sel++;
			};		
		std::vector<std::string> labels_trackEff;
					labels_trackEff = {"Total",
							              "Trigger",
							              "TofSt_Beta",
							              "L1 UnbHit",
							              "L1 UnbCharge",
							              "InsideINN",
							              "InsideL1"
							              "InnFidV",
							              "L1FidV",
							              "Tof #chi^{2} time",
							              "Tof #chi^{2} coo",
							              "Upper TofSt Charge",
							              "Lower TofSt Charge",
							              "Tot TofSt Charge"};
			auto ncuts_trackEff = labels_trackEff.size();
			auto counters_trackEff = new TH1D("counters_trackEff", "", ncuts_trackEff, -0.5, ncuts_trackEff - 0.5);
			auto fill_counters_trackEff = [&counters_trackEff, &icut_track, &labels_trackEff](Event &event) {
				 counters_trackEff->GetXaxis()->SetBinLabel(icut_track + 1, labels_trackEff.at(icut_track).c_str());
				 counters_trackEff->Fill(icut_track);
				 icut_track++;
			};		
		std::vector<std::string> labels_tofEff;
					labels_tofEff = {"Total",
										  "Trigger",
							              "InnerHitPattern",
							              "InnerNHits>4",
							              "NGoodCl",
							              "InnerChInRange",
							              "Inner #chi^{2}_{Y}",
							              "InnerL1 #chi^{2}_{Y}",
							              "L1HitCut",
							              "IsInL1",
							              "GoodSecTrack",
							              "L1 Charge In Range",
							              "Inner Charge Resolution <0.55",
							              "L1NormRes"};
				
			auto ncuts_tofEff = labels_tofEff.size();
			auto counters_tofEff = new TH1D("counters_tofEff", "", ncuts_tofEff, -0.5, ncuts_tofEff - 0.5);
			auto fill_counters_tofEff = [&counters_tofEff, &icut_tof, &labels_tofEff](Event &event) {
				 counters_tofEff->GetXaxis()->SetBinLabel(icut_tof + 1, labels_tofEff.at(icut_tof).c_str());
				 counters_tofEff->Fill(icut_tof);
				 icut_tof++;
			};		
		std::vector<std::string> labels_trigEff;
					labels_trigEff = {"Total",
										  "TofGoodPath",
							              "Beta",
							              "Upper Tof Charge",
							              "Inner Charge Pattern",
							              "NGoodCl>2",
							              "Inner Charge In Range",
							              "Inner #chi^{2}_{Y}",
							              "InnerL1 #chi^{2}_{Y}",
							              "L1 Hit",
							              "L1NormRes",
							              "L1 Charge In Range",
							              "Inner Charge Resolution <0.55",
							              "IsInsideL1"};				
			auto ncuts_trigEff = labels_trigEff.size();
			auto counters_trigEff = new TH1D("counters_trigEff", "", ncuts_trigEff, -0.5, ncuts_trigEff - 0.5);
			auto fill_counters_trigEff = [&counters_trigEff, &icut_trig, &labels_trigEff](Event &event) {
				 counters_trigEff->GetXaxis()->SetBinLabel(icut_trig + 1, labels_trigEff.at(icut_trig).c_str());
				 counters_trigEff->Fill(icut_trig);
				 icut_trig++;
			};		
		std::vector<std::string> labels_l1Eff;
					labels_l1Eff = {"Total",
							              "Trigger",
							              "TofGoodPath",
							              "Beta",
							              "Upper Tof Charge",
							              "Inner Charge Pattern",
							              "Inner Track NHits>4",
							              "NGoodCl",
							              "Inner Charge In Range",
							              "Inner #chi^{2}_{Y}",
							              "L1Fid",
							              "InnerFid",
							              "IsInL1",
							              "GoodTRDHits",
							              "InnerChRMS <0.55"};				
			auto ncuts_l1Eff = labels_l1Eff.size();
			auto counters_l1Eff = new TH1D("counters_l1Eff", "", ncuts_l1Eff, -0.5, ncuts_l1Eff - 0.5);
			auto fill_counters_l1Eff = [&counters_l1Eff, &icut_l1, &labels_l1Eff](Event &event) {
				 counters_l1Eff->GetXaxis()->SetBinLabel(icut_l1 + 1, labels_l1Eff.at(icut_l1).c_str());
				 counters_l1Eff->Fill(icut_l1);
				 icut_l1++;
			};
auto mysel =
	ns::Trigger::HasPhysicsTrigger().AddPostHook(fill_counters_Sel) && //OK
	ns::Tof::BetaInRange(0.0,inf,BTH).AddPostHook(fill_counters_Sel) && //OK
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.75,static_cast<float>(charge)+0.75,UTC).AddPostHook(fill_counters_Sel) && //OK
    (ns::Tof::GoodPathlength(0b0001) ||ns::Tof::GoodPathlength(0b0010) ).AddPostHook(fill_counters_Sel) && //OK
    ns::InnerTracker::HitPattern().AddPostHook(fill_counters_Sel) && //OK
    ns::InnerTracker::NGoodClustersGreaterThan(2, 0x10013D).AddPostHook(fill_counters_Sel) && //OK
    ns::InnerTracker::ChargeRMSLessThan(0.55, CRT).AddPostHook(fill_counters_Sel) && //OK
    ns::InnerTracker::ChargeInRange(static_cast<float>(charge)-0.3,static_cast<float>(charge)+0.7,CRT).AddPostHook(fill_counters_Sel) && //OK
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1).AddPostHook(fill_counters_Sel) && //OK
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN).AddPostHook(fill_counters_Sel) && //OK
    ns::Track::HitCut(1).AddPostHook(fill_counters_Sel) && //OK
	ns::TrackerLayer::ChargeInRange(1,0,static_cast<float>(charge)+0.8,CRT).AddPostHook(fill_counters_Sel) && //OK
	InnerTracker::L1NormResidualLessThan(10.0,FIT).AddPostHook(fill_counters_Sel); //OK
auto num_tof =
	ns::Tof::BetaInRange(0.0,inf,BTH) && //OK
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.75,static_cast<float>(charge)+0.75,UTC) && //OK
	(ns::Tof::GoodPathlength(0b0001) || ns::Tof::GoodPathlength(0b0010)); //OK
auto den_tof = 
	ns::Trigger::HasPhysicsTrigger().AddPostHook(fill_counters_tofEff) && //OK
	ns::InnerTracker::HitPattern().AddPostHook(fill_counters_tofEff) &&  //OK
	ns::InnerTracker::ChargeInRange(static_cast<float>(charge)-0.3,static_cast<float>(charge)+0.7,CRT).AddPostHook(fill_counters_tofEff) && //OK
	ns::InnerTracker::NGoodClustersGreaterThan(2, 0x10013D).AddPostHook(fill_counters_tofEff) && //OK
	ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1).AddPostHook(fill_counters_tofEff) && //OK
	ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN).AddPostHook(fill_counters_tofEff) && //OK
	ns::Track::HitCut(1).AddPostHook(fill_counters_tofEff) && //OK
	MySel::IsInsideL1(FIT, INN, 0).AddPostHook(fill_counters_tofEff) && //OK
	MySel::HasGoodSecondTrTrack(FIT, INN, 0.2).AddPostHook(fill_counters_tofEff) && //OK
	ns::TrackerLayer::ChargeInRange(1,0,static_cast<float>(charge)+0.8,CRT).AddPostHook(fill_counters_tofEff) && //OK
	ns::InnerTracker::ChargeRMSLessThan(0.55, CRT).AddPostHook(fill_counters_tofEff) && //OK
	InnerTracker::L1NormResidualLessThan(10.0,FIT).AddPostHook(fill_counters_tofEff); //OK	
auto num_l1 =
	ns::TrackerLayer::ChargeInRange(1,0,static_cast<float>(charge)+0.8,CRT) && 
	ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) && 
    ns::Track::L1FiducialVolume(FIT, IL1) && 
    ns::Track::HitCut(1) && 
    InnerTracker::L1NormResidualLessThan(10.0,FIT); 
auto den_l1 = 
	ns::Trigger::HasPhysicsTrigger().AddPostHook(fill_counters_l1Eff) && 
	ns::Tof::GoodPathlength(0b0011).AddPostHook(fill_counters_l1Eff) && 
	ns::Tof::BetaInRange(0.0,inf,BTH).AddPostHook(fill_counters_l1Eff) && 
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.2,static_cast<float>(charge)+0.2,UTC).AddPostHook(fill_counters_l1Eff) && 
	ns::InnerTracker::HitPattern().AddPostHook(fill_counters_l1Eff) && 
	ns::InnerTracker::NHitsGreaterThan(4, YSD).AddPostHook(fill_counters_l1Eff) &&
	ns::InnerTracker::NGoodClustersGreaterThan(2, 0x10013D).AddPostHook(fill_counters_l1Eff) &&
	ns::InnerTracker::ChargeInRange(static_cast<float>(charge)-0.2,static_cast<float>(charge)+0.2,CRT).AddPostHook(fill_counters_l1Eff) && 
	ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN).AddPostHook(fill_counters_l1Eff) && 
	ns::Track::L1FiducialVolume(FIT, INN).AddPostHook(fill_counters_l1Eff) && 
	ns::Track::InnerFiducialVolume(FIT, INN).AddPostHook(fill_counters_l1Eff) && 
	MySel::IsInsideL1(FIT, INN, 0).AddPostHook(fill_counters_l1Eff) && 
	trd::HasGoodTRDHits(charge).AddPostHook(fill_counters_l1Eff) &&
	ns::InnerTracker::ChargeRMSLessThan(0.55, CRT).AddPostHook(fill_counters_l1Eff); 
auto num_track =
	ns::InnerTracker::HitPattern() && //OK
	ns::InnerTracker::NGoodClustersGreaterThan(2, 0x10013D) && //OK
	ns::InnerTracker::ChargeRMSLessThan(0.55, CRT) && //OK
    ns::InnerTracker::ChargeInRange(static_cast<float>(charge)-0.3f,static_cast<float>(charge)+0.7,CRT) && //OK
    ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN); //OK
auto den_track =
	ns::Trigger::HasPhysicsTrigger().AddPostHook(fill_counters_trackEff) && //OK
	TofSt::tofBetaInRange(0.4,inf,BTH).AddPostHook(fill_counters_trackEff) && //OK
	track::UnbExtHitChargeStatus(EL1).AddPostHook(fill_counters_trackEff) && //OK
	track::UnbExtHitChargeInRange(EL1,static_cast<float>(charge),CRT,"").AddPostHook(fill_counters_trackEff) && //OK
	track::IsInsideInner(4.5).AddPostHook(fill_counters_trackEff) && //OK
	track::IsInsideL1(1).AddPostHook(fill_counters_trackEff) && //OK
	track::InnerFiducialVolume().AddPostHook(fill_counters_trackEff) &&
	track::L1FiducialVolume().AddPostHook(fill_counters_trackEff) &&
	TofSt::tofChi2TimeLessThan(2).AddPostHook(fill_counters_trackEff) && //OK
	TofSt::tofChi2CooLessThan(2).AddPostHook(fill_counters_trackEff) && //OK
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.3,static_cast<float>(charge)+0.3, UTC).AddPostHook(fill_counters_trackEff) && //OK
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.3,static_cast<float>(charge)+0.3, LTC).AddPostHook(fill_counters_trackEff) && //OK
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.3,static_cast<float>(charge)+0.3, TOT).AddPostHook(fill_counters_trackEff); //OK
auto den_trig =
	ns::Tof::GoodPathlength(0b0011).AddPostHook(fill_counters_trigEff) && //OK
	ns::Tof::BetaInRange(0.0,inf,BTH).AddPostHook(fill_counters_trigEff) && //OK
	ns::Tof::ChargeInRange(static_cast<float>(charge)-0.75,static_cast<float>(charge)+0.75,UTC).AddPostHook(fill_counters_trigEff) && //OK
    ns::InnerTracker::HitPattern().AddPostHook(fill_counters_trigEff) && //OK
	ns::InnerTracker::NGoodClustersGreaterThan(2, 0x10013D).AddPostHook(fill_counters_trigEff) && //OK
    ns::InnerTracker::ChargeInRange(static_cast<float>(charge)-0.3,static_cast<float>(charge)+0.7,CRT).AddPostHook(fill_counters_trigEff) && //OK
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1).AddPostHook(fill_counters_trigEff) && //OK
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN).AddPostHook(fill_counters_trigEff) && //OK
    ns::Track::HitCut(1).AddPostHook(fill_counters_trigEff) && //OK
    InnerTracker::L1NormResidualLessThan(10.0,FIT).AddPostHook(fill_counters_trigEff) && //OK
	ns::TrackerLayer::ChargeInRange(1,0,static_cast<float>(charge)+0.8,CRT).AddPostHook(fill_counters_trigEff) && //OKok
	ns::InnerTracker::ChargeRMSLessThan(0.55, CRT).AddPostHook(fill_counters_trigEff) && //OK
	trig::IsInsideL1(FIT,IL1,0).AddPostHook(fill_counters_trigEff); //OK

/////////////////// LOOP ///////////////////////////////
		for(Event& event : chain) {
			icut_sel = 0;
			icut_track = 0; 
			icut_tof = 0;
			icut_trig = 0;
			icut_l1 = 0;
			utime = event.header->UTCTime;
			if (!isMC) {
				if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
					continue;
				if (!NAIA::MatchAnyBit(event.header->Mask(), cat)) 	// check charge with TOF or Tracker
					continue;
				if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
					continue;

				cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[0][1];
				lt    		= chain.GetEventRTIInfo().LivetimeFraction;

				if(utime != oldtime) {
			    	for (int rbin = 1; rbin <=hist_rig->GetNbinsX(); rbin++) {
			        	auto rlowedge = hist_rig->GetBinLowEdge(rbin);
			         	auto rcenter  = hist_rig->GetBinCenter(rbin);
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
				int rbin = hist_rig->FindBin(reco_il1);
				double rlowedge = hist_rig->GetBinLowEdge(rbin);
				if (!isMC) {
					if(rlowedge > 1.2*cutoff) rigidity->Fill(utime, reco_il1);
				} else {
					mc_pass->Fill(reco_il1);
					mc_pass_gen->Fill(gen);
				}
			}
			////TOF EFF////
      		if(den_tof(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      		    int rbin = hist_rig->FindBin(reco_il1);
      		    double rlowedge = hist_rig->GetBinLowEdge(rbin);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) { 
      					sample_tof->Fill(utime, reco_il1);
      		        	if(num_tof(event)) pass_tof->Fill(utime, reco_il1);
      		    	}
      		    } else {
      		    	sample_tof->Fill(utime, reco_il1);
      		    	if (num_tof(event) ) pass_tof->Fill(utime, reco_il1);
      		    }
      		}
			////L1 EFF/////
      		if(den_l1(event)) {
      			reco_inn = event.trTrackBase->Rigidity[FIT][INN];
      		    int rbin = hist_rig->FindBin(reco_inn);
      		    double rlowedge = hist_rig->GetBinLowEdge(rbin);
      		    if (!isMC) {
      		    	if(rlowedge > 1.2*cutoff) {
      					sample_l1->Fill(utime, reco_inn);
      		        	if(num_l1(event)) pass_l1->Fill(utime, reco_inn);
      		   	 	}
      		   	} else {
      		   		sample_l1->Fill(utime, reco_inn);
      		   		if (num_l1(event)) pass_l1->Fill(utime, reco_inn);
      		   	}
      		}

			////TRACK EFF////
      		if(den_track(event)) {
      			FillHistos(*sample_tr,charge,Amass,weight,chain,event);
      			if(num_track(event)) FillHistos(*pass_tr,charge,Amass,weight,chain,event);
      		}
			////TRIGGER EFF/////
      		if(den_trig(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      			int rbin = hist_rig->FindBin(reco_il1);
        		double rlowedge = hist_rig->GetBinLowEdge(rbin);
        		if (!isMC) {
        			event.evSummary->RestorePhysBPatt(true);
        			double weight = event.evSummary->TriggerWeight();	
        		}
        		if (!isMC) { 
        			if(rlowedge > 1.2*cutoff) {
	      				if (physics(event)) {
	      					sample_trig->Fill(utime, reco_il1, weight);	
	      					phys_trig->Fill(utime, reco_il1, weight);
	      				}
                  		if (unbiased(event)) {
                  			sample_trig->Fill(utime, reco_il1, weight);
                  			unbiased_trig->Fill(utime, reco_il1);
          				}
          			}
          		} else {
          			if (physics(event)) {
          				sample_trig->Fill(utime, reco_il1, weight);
          				phys_trig->Fill(utime, reco_il1, weight);
          			}
          			if (unbiased(event)) {
          				unbiased_trig->Fill(utime,reco_il1,weight);
          				sample_trig->Fill(utime, reco_il1, weight);
          			}
          		}
     		}
      		if (isMC) tree->Fill();	
		} //event loop	
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(counters_Sel, "counters_Sel");
	outfile->WriteTObject(counters_tofEff, "counters_tofEff");
	outfile->WriteTObject(counters_trigEff, "counters_trigEff");
	outfile->WriteTObject(counters_l1Eff, "counters_l1Eff");
	outfile->WriteTObject(counters_trackEff, "counters_trackEff");
	outfile->WriteTObject(rigidity, "rigidity");
	outfile->WriteTObject(lvt_25, "lvt_25");
	outfile->WriteTObject(sample_tof, "sample_tof");
	outfile->WriteTObject(pass_tof, "pass_tof");
	outfile->WriteTObject(sample_l1, "sample_l1");
	outfile->WriteTObject(pass_l1, "pass_l1");
	outfile->WriteTObject(sample_tr, "sample_tr");
	outfile->WriteTObject(pass_tr, "pass_tr");
	outfile->WriteTObject(sample_trig, "sample_trig");
	outfile->WriteTObject(phys_trig, "phys_trig");
	outfile->WriteTObject(unbiased_trig, "unbiased_trig");
	if (isMC) {
	    //outfile->WriteTObject(tree, "reco_gen");
	    outfile->WriteTObject(mc_samp, "mc_samp");
	    outfile->WriteTObject(mc_pass_gen, "mc_pass_gen");
	    outfile->WriteTObject(mc_pass, "mc_pass");
	    //outfile->WriteTObject(corr_mc_pass, "corr_mc_pass");
	}
	outfile->Close();
	} //valid input
} //main


//FUNCTION DEFINITION
void FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain,  NAIA::Event &event) {
  bool IsMC = chain.IsMC();
  float IGRFCutoff;
  if (!IsMC) {
    const NAIA::RTIInfo &rti_info = chain.GetEventRTIInfo();
    IGRFCutoff = rti_info.MaxIGRFCutoff[3][1];
  } else
    IGRFCutoff = event.mcTruthBase->Primary.GetGenMomentum() / event.mcTruthBase->Primary.Z;
  auto utime = event.header->UTCTime;
  float Beta = event.tofBaseSt->Beta[NAIA::Tof::BetaType::BetaH];
  const float amu = 0.93146;
  double mass = Amass * amu;
  double BetaRig = mass * Beta / sqrt(1 - Beta * Beta) / charge;

  float EnergyD = 0;
  if (NAIA::ContainsKeys(event.ecalBase->Energy, NAIA::Ecal::EnergyRecoType::EnergyD)) {
    EnergyD = event.ecalBase->Energy[NAIA::Ecal::EnergyRecoType::EnergyD];
  }
  if (BetaRig <= 5) {
    histo.Fill(utime, BetaRig, weight);
  } else if (IGRFCutoff > 5 && IGRFCutoff <= 30) {
    histo.Fill(utime, IGRFCutoff, weight);
  } else if (EnergyD > 30 && Efficiency::TrTrackEffSel::IsInsideEcal()(event)) {
    histo.Fill(utime, EnergyD, weight);
  }
}

double GetDataMass(unsigned int charge) {
  double A = 0;
  switch (charge) {
  case 1:
    A = 1;
    break;
  case 2:
    A = 4;
    break;
  case 3:
    A = 6.5;
    break;
  case 4:
    A = 7.6;
    break;
  case 5:
    A = 10.65;
    break;
  case 6:
    A = 12;
  case 7:
    A = 14;
    break;
  case 8:
    A = 16;
    break;
  case 10:
    A = 20;
    break;
  case 12:
    A = 24;
    break;
  case 14:
    A = 28;
    break;
  case 16:
    A = 32;
    break;
  }
  return A;
}
