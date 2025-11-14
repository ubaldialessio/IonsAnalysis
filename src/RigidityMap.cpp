#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int main(int argc, char *argv[]) {

      // Time markers
    unsigned int t0 = 1305853512, t1 = 1385487767, t2 = 1447346927, t3 = 1456503197,
                 t4 = 1582034309, t5 = 1582037855, t6 = 1620025528, t7 = 1635856717,
                 t8 = 1675341999;
	if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root>  \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString exename = argv[0], infilename=argv[2],
		      outname=argv[3],
		      out, ionName=getIonName(charge);

out = StreamUtility::getOutputDir(charge,exename,outname);
std::cout << out << std::endl;
auto map       = new TH2D("", ";Cutoff (GV); Inner R (GV)", nRbins_HighZ-1, Rbins_HighZ, nRbins_HighZ-1, Rbins_HighZ);
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
double utime,cutoff,Amass,weight,reco_il1,reco_beta=0;
Amass = GetDataMass(charge); // Getting atomic mass number
NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;
auto mysel = Mysel(charge);
    for(Event& event : chain) {
        utime = event.header->UTCTime;
        if (utime >= t6 && utime < t7) continue; // photon-polarization
        if (isRunBad(utime)) continue;           // bad runs
        if (!NAIA::MatchAnyBit(event.header->Mask(),
            NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof))
          continue;
        if (!NSL::Selections::DefaultRTISelection(chain.GetEventRTIInfo()))
          continue;
        if (mysel(event)) {
            cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[1][1];
            reco_il1  = event.trTrackBase->Rigidity[FIT][INN];
            map->Fill(cutoff,reco_il1);
        }
    }
  auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(map, "map");
	outfile->Close();
}
