#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

std::array<int, 2> PIDtoParticle(int pID);
int pidToNucleusCharge(int pID);

int main(int argc, char *argv[]) {

	TH1::SetDefaultSumw2();
	auto hRig_IL1 	  = new TH1D ("", ";R (GV);Counts", nRbins_HighZ-1, Rbins_HighZ);
	auto mc_samp      = (TH1D*)hist_log->Clone("mc_samp");
    auto mc_pass      = (TH1D*)hist_log->Clone("mc_pass");
    auto mc_pass_gen  = (TH1D*)hist_log->Clone("mc_pass_gen");
	auto mc_pass_gen_acc4  = (TH1D*)hist_log->Clone("mc_pass_gen_acc4");
	auto mc_pass_acc4      = (TH1D*)hist_log->Clone("mc_pass_acc4");
    auto mc_pass_gen_acc7  = (TH1D*)hist_log->Clone("mc_pass_gen_acc7");
	auto mc_pass_acc7      = (TH1D*)hist_log->Clone("mc_pass_acc7");
    auto events_acc4  = new TTree();
    auto events_acc7  = new TTree();
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	if (argc < 4) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> <secTrack>\n", argv[0]);
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
    out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/Fragments/"+outname;
	if (sec_track=="y") out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/Fragments/sec_track_"+outname;
	unsigned int utime;
	float reco_il1;
	double gen;
    unsigned short nAcc_counters = 0;
    unsigned int runNum=-999, evNum=-999;

    events_acc4->Branch("runNum",&runNum,"runNum/I");
    events_acc4->Branch("evNum",&evNum,"evNum/I");
    events_acc7->Branch("runNum",&runNum,"runNum/I");
    events_acc7->Branch("evNum",&evNum,"evNum/I");


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
	auto mysel = Mysel(charge);
    auto hasGoodSecondTrack = GoodSecTrack();
    auto secTrackOnDiagonal = MySel::IsSecondTrackOnDiagonal(FIT,IL1,0.3);
    auto l1ClusterCut = MySel::L1ClusterCutBetween(10,30,YSD,ALL_LAYER);


    for(Event& event : chain) {
        auto header = event.header;
		utime = header->UTCTime;
        runNum= header->Run;
        evNum= header->EventNo;

		gen = event.mcTruthBase->Primary.GetGenMomentum()/event.mcTruthBase->Primary.Z;
        
		//Check that the interaction happens only Above L1
	    bool frag_to_target = false;
        for (const auto &sec : event.mcTruthPlus->Secondaries) {
            auto Z = sec.Z;
            auto Position = sec.Position;
            auto GenPosition = sec.GetGenPosition();
            double Gen_x = GenPosition.X();
            double Gen_y = GenPosition.Y();
            double Gen_z = GenPosition.Z();
            //secondary Z match charge (nucleus under study) && secondary is generated above L1
            if(Z == charge && Gen_z >= 159.04) {
                frag_to_target = true;
            }
        }
        if(!frag_to_target)
            continue;

		////MY SELECTION////
		if (mysel(event)) {
            nAcc_counters = event.evSummary->NAcc;
			/*if (sec_track=="y") {
				std::cout << "Searching for good second track ... " << std::endl;
				if (hasGoodSecondTrack(event)) {
					std::cout << "Good second track found" << std::endl;
					if (secTrackOnDiagonal(event)) {
						std::cout << "Skipping event" << std::endl;
						continue;
					}
				}
			}*/
            //if (!l1ClusterCut(event)) continue;
			reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
			int rbin = hist_rig_highZ->FindBin(reco_il1);
			double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
            mc_pass->Fill(reco_il1);
            mc_pass_gen->Fill(gen);
            if (nAcc_counters<=4) {
                mc_pass_acc4->Fill(reco_il1);
			    mc_pass_gen_acc4->Fill(gen);
                events_acc4->Fill();
            }
            if (nAcc_counters<=7) {
                mc_pass_acc7->Fill(reco_il1);
			    mc_pass_gen_acc7->Fill(gen);
                events_acc7->Fill();
            }
			hRig_IL1->Fill(reco_il1);
		}
	} //event loop	
    auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(mc_samp, "mc_samp");
    outfile->WriteTObject(mc_pass_acc4, "mc_pass");
	outfile->WriteTObject(mc_pass_gen_acc4, "mc_pass_gen");
	outfile->WriteTObject(mc_pass_acc4, "mc_pass_acc4");
	outfile->WriteTObject(mc_pass_gen_acc4, "mc_pass_gen_acc4");
    outfile->WriteTObject(mc_pass_acc7, "mc_pass_acc7");
	outfile->WriteTObject(mc_pass_gen_acc7, "mc_pass_gen_acc7");
	outfile->WriteTObject(hRig_IL1,"hRig_IL1");
    outfile->WriteTObject(events_acc4,"events_acc4");
    outfile->WriteTObject(events_acc7,"events_acc7");
	outfile->Close();
}

//this works only for naia 1.2.0
int pidToNucleusCharge(int pID) {
    switch(pID) {
        case 11:   // e-
            return -1;
        case -2112: // anti-neutron
            return 0;
        case 2212: // p
            return 1;
        case +1000020030: // HE3C         Z:   2 A:   3 mass: 2.809230
            return 2;
        case +1000020040: // ALPHA
            return 2;
        case +1000030070: // LI7
            return 3;
        case +1000040070: // BE7
            return 4;
        case +1000040090: // BE9
            return 4;
        case +1000040100: // BE10
            return 4;
        case +1000050100: // B10
            return 5;
        case +1000050110: // B11
            return 5;
        case +1000060120: // C12
            return 6; 
        case +1000060130: // C13
            return 6;
        case +1000070140: // N14
            return 7;
        case +1000070150: // N15
            return 7;
        case +1000080160: // O16
            return 8;     
        case +1000080170: // O17
            return 8;
        case +1000080180: // O18
            return 8;   
        case +1000150310: // P31
            return 15;
        case +1000160320: // S32
            return 16;
        case +1000170350: // CL35
            return 17;
        case +1000180360: // AR36
            return 18;
        case +1000190390: // K39
            return 19;
        case +1000200400: // CA40
            return 20;
        case +1000250550: // MN55
            return 25;
        case +1000260560: // FE56
            return 26;
    }
    return 0;
}

std::array<int, 2> PIDtoParticle(int pID) {
  //lxplus -> /cvmfs/ams.cern.ch/Offline/vdev/F/Geant3_utils.F
  // maybe better with a map
  if (abs(pID)>10000) return {pID/abs(pID)*int(abs(pID)/100)%100,abs(pID)%100};
  switch (pID) { 
    case   2: return { 1,  0}; // positron 
    case   3: return {-1,  0}; // electron
    case  13: return { 0,  1}; // neutron
    case  14: return { 1,  1}; // proton
    case  15: return {-1,  1}; // antiproton
    case  45: return { 1,  2}; // H2
    case  46: return { 1,  3}; // H3
    case  47: return { 2,  4}; // He4 
    case  49: return { 2,  3}; // He3 
    case  61: return { 3,  6}; // Li6
    case  62: return { 3,  7}; // Li7
    case  63: return { 4,  7}; // Be7
    case  64: return { 4,  8}; // Be9
    case 114: return { 4, 10}; // Be10
    case  65: return { 5, 10}; // B10
    case  66: return { 5, 11}; // B11
    case  67: return { 6, 12}; // C12
    case  68: return { 7, 14}; // N14 
    case 118: return { 7, 15}; // N15  
    case  69: return { 8, 16}; // O16  
    case  70: return { 9, 19}; // F19  
    case  71: return {10, 20}; // Ne20 
    case  72: return {11, 23}; // Na23 
    case  73: return {12, 24}; // Mg24 
    case  74: return {13, 27}; // Al27 
    case  75: return {14, 28}; // Si28 
    case  76: return {15, 31}; // P31  
    case  77: return {16, 32}; // S32  
    case  78: return {17, 35}; // Cl35 
    case  79: return {18, 36}; // Ar36 
    case  80: return {19, 39}; // K39  
    case  81: return {20, 40}; // Ca40 
    case  82: return {21, 45}; // Sc45 
    case  83: return {22, 48}; // Ti48 
    case  84: return {23, 51}; // V51  
    case  85: return {24, 52}; // Cr52 
    case  86: return {25, 55}; // Mn55 
    case  87: return {26, 56}; // Fe56 
    case  88: return {27, 59}; // Co59 
    case  89: return {28, 58}; // Ni58 
    case  90: return {29, 63}; // Cu63 
    case  91: return {30, 64}; // Zn64 
    case  92: return {32, 74}; // Ge74 
    case  93: return {34, 80}; // Se80 
    case  94: return {36, 84}; // Kr84 
    case  95: return {38, 88}; // Sr88 
    case  96: return {40, 90}; // Zr90 
    case  97: return {42, 98}; // Mo98 
    case  98: return {46,106}; // Pd106
    case  99: return {48,114}; // Cd114
    case 100: return {50,120}; // Sn120
    case 101: return {54,132}; // Xe132
    case 102: return {56,138}; // Ba138
    case 103: return {58,140}; // Ce140
    case 104: return {62,152}; // Sm152
    case 105: return {66,164}; // Dy164
    case 106: return {70,174}; // Yb174
    case 107: return {74,184}; // W184 
    case 108: return {78,194}; // Pt194
    case 109: return {79,197}; // Au197
    case 110: return {80,202}; // Hg202
    case 111: return {82,208}; // Pb208
    case 112: return {92,238}; // U238
  }
  return {0,0};
}