#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto sample_tof_den   		 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_tof_den     		 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto sample_tof_den_ltof      = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_tof_den_ltof        = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);   
	auto sample_tof_den_ltof_l9_05= new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_tof_den_ltof_l9_05  = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);   
	auto sample_tof_den_ltof_l9_1 = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ);
	auto pass_tof_den_ltof_l9_1   = new TH2D("", ";Utime;R(GV)", nTbins-1, Tbins, nRbins_HighZ-1, Rbins_HighZ); 
	auto qtof_den		           = new TH1D("", ";Q_tof;",550,0,28); 
	auto qtof_den_ltof		   	   = new TH1D("", ";Q_tof;",550,0,28);
	auto qtof_den_ltof_l9_05	   = new TH1D("", ";Q_tof;",550,0,28); 
	auto qtof_den_ltof_l9_1	       = new TH1D("", ";Q_tof;",550,0,28); 
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
    unsigned int utime, icut_tof = 0;
	float reco_beta, reco_il1, qtof;
	double cutoff, weight = 1;
		std::vector<std::string> labels_tofEff;
					labels_tofEff = {"Total",
										  "Trigger",
							              "InnerHitPattern",
										  "N hits y > 4",
							              "Inn tracker charge",
							              "L1NormRes<10",
										  "L1 Charge asymmetry < 0.2",
										  "L1 Charge In Range",
										  "L1 hit",
							              "IL1 #chi^{2}_{Y}",
							              "Inner #chi^{2}_{Y}",
							              "L1HitCut",
							              "L1 fid volume",
										  "Inn fid volume",								
							              "NoGoodSecTrack",
										  "Ltof charge cut"};
			auto ncuts_tofEff = labels_tofEff.size();
			auto counters_tofEff = new TH1D("counters_tofEff", "", ncuts_tofEff, -0.5, ncuts_tofEff - 0.5);
			auto fill_counters_tofEff = [&counters_tofEff, &icut_tof, &labels_tofEff](Event &event) {
				 counters_tofEff->GetXaxis()->SetBinLabel(icut_tof + 1, labels_tofEff.at(icut_tof).c_str());
				 counters_tofEff->Fill(icut_tof);
				 icut_tof++;
			};	
	auto num_tof = Num_tof(charge);
	auto den_tof = Den_tof(charge);
	auto den_tof_ltof = Den_tof_ltof(charge);
	auto den_tof_ltof_l9_05 = Den_tof_ltof_l9_05(charge);
	auto den_tof_ltof_l9_1 = Den_tof_ltof_l9_1(charge);

	for(Event& event : chain) { ////LOOP
        icut_tof = 0; 
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
		////TOF EFF////
      	if(den_tof(event)) {
      		reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      		int rbin = hist_rig_highZ->FindBin(reco_il1);
      		double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			qtof = event.tofBase->Charge[UTC];
      		if (!isMC) {
      		    if(rlowedge > 1.2*cutoff) { 
      				sample_tof_den->Fill(utime, reco_il1);
					qtof_den->Fill(qtof);
      		        if(num_tof(event)) pass_tof_den->Fill(utime, reco_il1);
      		    }
      		} else {
      		    sample_tof_den->Fill(utime, reco_il1);
				qtof_den->Fill(qtof);
      		    if (num_tof(event) ) {
					pass_tof_den->Fill(utime, reco_il1);
				}
      		}
      	}
		if(den_tof_ltof(event)) {
      		reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      		int rbin = hist_rig_highZ->FindBin(reco_il1);
      		double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			qtof = event.tofBase->Charge[UTC];
      		if (!isMC) {
      		    if(rlowedge > 1.2*cutoff) { 
      				sample_tof_den_ltof->Fill(utime, reco_il1);
					qtof_den_ltof->Fill(qtof);
      		        if(num_tof(event)) pass_tof_den_ltof->Fill(utime, reco_il1);
      		    }
      		} else {
      		    sample_tof_den_ltof->Fill(utime, reco_il1);
				qtof_den_ltof->Fill(qtof);
      		    if (num_tof(event) ) {
					pass_tof_den_ltof->Fill(utime, reco_il1);
				}
      		}
      	}
		if(den_tof_ltof_l9_05(event)) {
      		reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      		int rbin = hist_rig_highZ->FindBin(reco_il1);
      		double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			qtof = event.tofBase->Charge[UTC];
      		if (!isMC) {
      		    if(rlowedge > 1.2*cutoff) { 
      				sample_tof_den_ltof_l9_05->Fill(utime, reco_il1);
					qtof_den_ltof_l9_05->Fill(qtof);
      		        if(num_tof(event)) pass_tof_den_ltof_l9_05->Fill(utime, reco_il1);
      		    }
      		} else {
      		    sample_tof_den_ltof_l9_05->Fill(utime, reco_il1);
				qtof_den_ltof_l9_05->Fill(qtof);
      		    if (num_tof(event) ) {
					pass_tof_den_ltof_l9_05->Fill(utime, reco_il1);
				}
      		}
      	}
		if(den_tof_ltof_l9_1(event)) {
      		reco_il1 = event.trTrackBase->Rigidity[FIT][IL1];
      		int rbin = hist_rig_highZ->FindBin(reco_il1);
      		double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			qtof = event.tofBase->Charge[UTC];
      		if (!isMC) {
      		    if(rlowedge > 1.2*cutoff) { 
      				sample_tof_den_ltof_l9_1->Fill(utime, reco_il1);
					qtof_den_ltof_l9_1->Fill(qtof);
      		        if(num_tof(event)) pass_tof_den_ltof_l9_1->Fill(utime, reco_il1);
      		    }
      		} else {
      		    sample_tof_den_ltof_l9_1->Fill(utime, reco_il1);
				qtof_den_ltof_l9_1->Fill(qtof);
      		    if (num_tof(event) ) {
					pass_tof_den_ltof_l9_1->Fill(utime, reco_il1);
				}
      		}
      	}
    } ////LOOP
    auto outfile = new TFile(out, "recreate");
    outfile->WriteTObject(counters_tofEff, "counters_tofEff");
    outfile->WriteTObject(sample_tof_den, "sample_tof_den");
    outfile->WriteTObject(pass_tof_den, "pass_tof_den");
	outfile->WriteTObject(sample_tof_den_ltof, "sample_tof_den_ltof");
    outfile->WriteTObject(pass_tof_den_ltof, "pass_tof_den_ltof");
	outfile->WriteTObject(sample_tof_den_ltof_l9_05, "sample_tof_den_ltof_l9_05");
    outfile->WriteTObject(pass_tof_den_ltof_l9_05, "pass_tof_den_ltof_l9_05");
	outfile->WriteTObject(sample_tof_den_ltof_l9_1, "sample_tof_den_ltof_l9_1");
    outfile->WriteTObject(pass_tof_den_ltof_l9_1, "pass_tof_den_ltof_l9_1");
	outfile->WriteTObject(qtof_den,"qtof_den");
	outfile->WriteTObject(qtof_den_ltof,"qtof_den_ltof");
	outfile->WriteTObject(qtof_den_ltof_l9_05,"qtof_den_ltof_l9_05");
	outfile->WriteTObject(qtof_den_ltof_l9_1,"qtof_den_ltof_l9_1");
    outfile->Close();
}