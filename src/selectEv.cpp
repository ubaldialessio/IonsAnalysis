#include "selection.h"

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

	auto noTriggers  = trg::HasTrigger(0);                	// 00000000
	auto unbiased    = trg::HasTrigger(1)  ||             	// 00000001 (only unbiased ToF)
	                   trg::HasTrigger(64) ||            	// 01000000 (only unbiased ECAL)
	                   trg::HasTrigger(65);              	// 01000001 (unbiased ToF AND unbiased ECAL)
	auto physics     = ns::Trigger::HasPhysicsTrigger();    // 00111110 (Any non-prescaled trigger)


	// this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	if (argc < 2) {
	    printf("Usage: n");
	    printf("%s <input root file> <output> n", argv[0]);
	    return 1;
	}
	TString infilename = argv[1], out = argv[2];

	//check if is a valid input
    bool validInput = true;
	//check(); 

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

		unsigned int utime, oldtime;
		float reco_il1, reco_inn, beta, lt;
		double gen, cutoff;
		double Amass = 1;
		double weight = 1;
		Amass = GetDataMass(6); // Getting atomic mass number

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
		

/////////////////// LOOP ///////////////////////////////
		for(Event& event : chain) {
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
				gen = event.mcTruthBase->Primary.GetGenMomentum();
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
      			FillHistos(*sample_tr,6,Amass,weight,chain,event);
      			if(num_track(event)) FillHistos(*pass_tr,6,Amass,weight,chain,event);
      		}


			////TRIGGER EFF/////
      		if(den_trig(event)) {
      			reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      			int rbin = hist_rig->FindBin(reco_il1);
        		double rlowedge = hist_rig->GetBinLowEdge(rbin);
        		event.evSummary->RestorePhysBPatt(true);
        		double weight = event.evSummary->TriggerWeight();
        		if (!isMC) {
        			if(rlowedge > 1.2*cutoff) {
	      				if (physics(event) || unbiased(event)) sample_trig->Fill(utime, reco_il1, weight);
          				if (physics(event))  phys_trig->Fill(utime, reco_il1, weight);
                  		if (unbiased(event)) unbiased_trig->Fill(utime, reco_il1);
          			}
          		} else {
          			if (physics(event) || unbiased(event) ) sample_trig->Fill(utime, reco_il1);
          			if (physics(event) ) phys_trig->Fill(utime, reco_il1);
          			if (unbiased(event) ) unbiased_trig->Fill(utime,reco_il1);
          		}
     		}
      		if (isMC) tree->Fill();	
		} //event loop
	
	//WRITE//
	auto outfile = new TFile(out, "recreate");
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
