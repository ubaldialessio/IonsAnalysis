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
	auto L1charge     = new TH1D("", "L1 Charge (Z)", 300, 0, 20);
	auto allRigL1	  = new TH1D("", ";R(GV)", nRbins-1, Rbins);
	auto genRig		  = (TH1D*)hist_log->Clone("genRig");

	// setup binning template useed for L1/L2 templates.
	const int nRigbins_LightIons = 24;
	const double Rigbins_LightIons[nRigbins_LightIons] = {0.8,  1.16, 1.51, 1.92, 2.40,  2.97, 3.64, 4.43,
	                                                  5.37, 6.47, 7.76, 9.26, 11.00, 13.0, 15.3, 18.0,
	                                                  21.1, 24.7, 28.8, 33.5, 38.9,  45.1, 52.2, 5000};
	
	const int nRigbins_HeavyIons = 21;
	const double Rigbins_HeavyIons[nRigbins_HeavyIons] = {0.8,  2.40, 2.97, 3.64, 4.43, 5.37, 6.47, 7.76, 9.26, 11.00, 13.0,
	                                                  15.3, 18.0, 21.1, 24.7, 28.8, 33.5, 38.9, 45.1, 52.2, 5000};
	TList *l1TemplateList = new TList();
	TList *l2TemplateList = new TList();
	auto rigBinTemplate = new TH1D("RigTemplateLightIons", "R (GV)", nRigbins_LightIons - 1, Rigbins_LightIons);
	
	for (Int_t ibin = 0; ibin < rigBinTemplate->GetNbinsX(); ibin++) 
	    l1TemplateList->Add(new TH1D(Form("L1Template_%03i", ibin),
	                               Form("%5.3f < R (GV) < %5.3f;Q_{L1X};Counts", rigBinTemplate->GetBinLowEdge(ibin + 1),
	                                    rigBinTemplate->GetBinLowEdge(ibin + 2)),
	                               600, 0, 20));
	                               
	for (Int_t ibin = 0; ibin < rigBinTemplate->GetNbinsX(); ibin++) 
    	l2TemplateList->Add(new TH1D(Form("L2Template_%03i", ibin),
                               Form("%5.3f < R (GV) < %5.3f;Q_{L2X};Counts", rigBinTemplate->GetBinLowEdge(ibin + 1),
                                    rigBinTemplate->GetBinLowEdge(ibin + 2)),
                               600, 0, 20));

	// this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	if (argc < 4) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> <AboveL1/BelowL1> \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString infilename=argv[2],
		    outname=argv[3],
		    out, ionPath=getIonName(charge), where=argv[4];
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
		if (where.Contains("AboveL1") && isMC==true) out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/AboveL1/"+ionPath+"/"+outname;
		if (where.Contains("BelowL1") ) out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"+ionPath+"/"+outname;
		unsigned int utime, oldtime;
		float reco_il1, reco_inn, beta, lt, l1charge;
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
		
auto mysel =
	ns::Trigger::HasPhysicsTrigger() && //OK
	ns::Tof::BetaInRange(0.0,inf,BTH) && //OK
	ns::PG::Tof::ChargeInRange(static_cast<float>(charge),UTC) && //OK
    (ns::Tof::GoodPathlength(0b0001) ||ns::Tof::GoodPathlength(0b0010) ) && //OK
    ns::InnerTracker::HitPattern() && //OK
    ns::InnerTracker::NGoodClustersGreaterThan(2, 0x10013D) && //OK
    ns::InnerTracker::ChargeRMSLessThan(0.55, CRT) && //OK
    ns::PG::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT) && //OK
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) && //OK
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) && //OK
    ns::Track::HitCut(1) && //OK
	ns::PG::TrackerLayer::ChargeInRange(1,static_cast<float>(charge),CRT) && //OK
	InnerTracker::L1NormResidualLessThan(10.0,FIT); //OK
auto l1temp =
    BuildTemplatesSel::HitPattern() &&
	BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge), UTC) &&       
	BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge), LTC) &&
	ns::TrackerLayer::ChargeStatus(1) &&
    BuildTemplatesSel::InnerTrackerChargeInRange(static_cast<float>(charge), CRT) &&
    BuildTemplatesSel::InnerTrackerChargeRMSLessThan(0.07, CRT);
	if (charge < 9)
    l1temp &= BuildTemplatesSel::InnerTrackerNtrackLessThan(3) &&
  	BuildTemplatesSel::HasGoodSecondTrTrack(FIT,INN, 0.2);
  	//BuildTemplatesSel::HasGoodTRDHits(charge);
auto l2temp = 
    BuildTemplatesSel::HitPattern() &&
    BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge),UTC) &&
    BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge),LTC) &&
    BuildTemplatesSel::InnerTrackerChargeInRange(static_cast<float>(charge), CRT) &&
    BuildTemplatesSel::InnerTrackerChargeRMSLessThan(0.07, CRT) && ns::TrackerLayer::ChargeStatus(1) &&
    BuildTemplatesSel::TrackerLayerChargeInRange(1, static_cast<float>(charge), CRT) && ns::TrackerLayer::ChargeStatus(2);
	if (charge < 9)
    l2temp &= BuildTemplatesSel::InnerTrackerNtrackLessThan(3) &&
	BuildTemplatesSel::HasGoodSecondTrTrack(FIT,INN, 0.2);



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
				gen = event.mcTruthBase->Primary.GetGenMomentum()/event.mcTruthBase->Primary.Z;
				allRigL1->Fill(event.trTrackBase->Rigidity[FIT][IL1]);
				genRig->Fill(gen);
			}	
			////MY SELECTION////
			if (where.Contains("AboveL1") ) {
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
			}

			if (where.Contains("BelowL1") ) {
				if (l1temp(event)) {
					auto Rigidity_IL1 = event.trTrackBase->Rigidity[FIT][IL1];
					int bin = rigBinTemplate->FindBin(Rigidity_IL1);
					if (bin < 1 || bin > rigBinTemplate->GetNbinsX())
						continue;
					auto L1charge = event.trTrackBase->LayerChargeXY[0][CRT];
					static_cast<TH1D *>(l1TemplateList->At(bin - 1))->Fill(L1charge);
				}
				if (l2temp(event)) {
					auto Rigidity_IL1 = event.trTrackBase->Rigidity[FIT][IL1];
					int bin = rigBinTemplate->FindBin(Rigidity_IL1);
					if (bin < 1 || bin > rigBinTemplate->GetNbinsX())
						continue;
					auto L2charge = event.trTrackBase->LayerChargeXY[1][CRT];
				    static_cast<TH1D *>(l2TemplateList->At(bin - 1))->Fill(L2charge);
				}
			}
      		if (isMC) tree->Fill();	
		} //event loop	
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(rigidity, "rigidity");
	outfile->WriteTObject(lvt_25, "lvt_25");
	outfile->WriteTObject(L1charge, "L1charge");
	outfile->WriteTObject(l1TemplateList, Form("L1TemplateList_%i", charge), "Overwrite");
	outfile->WriteTObject(l2TemplateList, Form("L2TemplateList_%i", charge), "Overwrite");
	if (isMC) {
	    outfile->WriteTObject(mc_samp, "mc_samp");
	    outfile->WriteTObject(mc_pass_gen, "mc_pass_gen");
	    outfile->WriteTObject(mc_pass, "mc_pass");
	    outfile->WriteTObject(allRigL1, "allRigL1");
	    outfile->WriteTObject(genRig, "genRig");
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
