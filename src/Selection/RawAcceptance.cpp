#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto mc_samp      = (TH1D*)hist_log->Clone("mc_samp");
	auto mc_pass_gen  = (TH1D*)hist_log->Clone("mc_pass_gen");
	auto mc_pass      = (TH1D*)hist_log->Clone("mc_pass");

	// this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> -sec_track <y/n>\n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString exename = argv[0],infilename=argv[2], outname=argv[3], out, ionPath=getIonPath(charge), sec_track=argv[4];
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
		if (sec_track=="y") outname = "sec_track_"+outname;
		out = StreamUtility::getOutputDir(charge,exename,outname);
		unsigned int utime;
		float reco_il1;
		double gen;
		///////////////////////////////////////
		auto mcchain = chain.GetFileInfoTree();
		auto mcinfo = new NAIA::MCFileInfo();
		mcchain->SetBranchAddress("MCFileInfo", &mcinfo);
		for(unsigned long long i = 0; i < mcchain->GetEntries(); ++i) {
			mcchain->GetEntry(i);
		    double rminDC = mcinfo->GetRMin();
		    double rmaxDC = mcinfo->GetRMax();
		    for(int ibin = 1; ibin <= mc_samp->GetNbinsX(); ++ibin){
		        double rmin = mc_samp->GetBinLowEdge(ibin), rmax = mc_samp->GetBinLowEdge(ibin+1);
		        if (rmin < rminDC) rmin = rminDC;
		        if (rmax > rmaxDC) rmax = rmaxDC;
		        if (rmin > rmaxDC || rmax < rminDC) continue;
		        mc_samp->AddBinContent(ibin, mcinfo->GetNGen() * TMath::Log(rmax / rmin) / TMath::Log(rmaxDC / rminDC));
		    }
		}
auto mysel  = Mysel(charge);
auto hasGoodSecondTrack = GoodSecTrack();
auto secTrackOnDiagonal = MySel::IsSecondTrackOnDiagonal(FIT,IL1,0.3);

    for(Event& event : chain) {
		utime = event.header->UTCTime;
		gen = event.mcTruthBase->Primary.GetGenMomentum()/event.mcTruthBase->Primary.Z;
		////MY SELECTION////
		if (mysel(event)) {
			if (sec_track=="y") {
				if (hasGoodSecondTrack(event)) {
					if (secTrackOnDiagonal(event)) {
						continue;
					}
				}
			}
			reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
			int rbin = hist_rig_highZ->FindBin(reco_il1);
			double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
			mc_pass->Fill(reco_il1);
			mc_pass_gen->Fill(gen);
		}
	} //event loop	
    auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(mc_samp, "mc_samp");
	outfile->WriteTObject(mc_pass_gen, "mc_pass_gen");
	outfile->WriteTObject(mc_pass, "mc_pass");
	outfile->Close();
}
