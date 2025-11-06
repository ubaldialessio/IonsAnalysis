#include "binning.h"
#include "utils.h"

int main(int argc, char *argv[]) {

    std::vector<unsigned int> UnderStudyCharge = {8,  9,  10, 
                                                  14, 15, 16,
                                                  17, 18, 19, 20, 21, 22, 23, 24, 25, 26};

    TH1::SetDefaultSumw2();	

    if (argc < 3) {
        printf("Usage: \n");
        printf("%s <path/to/input.root/.txt> <output.root> \n", argv[0]);
        return 1;
    }

    TString executable = argv[0],
            infilename = argv[1],
            outname = argv[2];
    TString exename = executable(executable.Last('/') + 1, executable.Length());
    // Setup chain
    NAIA::NAIAChain chain;
    if (infilename.Contains(".root")) {
        chain.Add(infilename.Data());
    } else if (infilename.Contains(".txt")) {
        std::ifstream infilelist(infilename.Data());
        TString bufname;
        while (infilelist >> bufname)
            if (bufname.Contains(".root"))
                chain.Add(bufname.Data());
    }
    chain.SetupBranches();
    bool isMC = chain.IsMC();

    // Prepare ChargeAnalysisManager
    ChargeAnalysisManager manager(
        UnderStudyCharge,
        nTbins, Tbins,
        nRbins_HighZ, Rbins_HighZ,isMC,exename
    );

    // Map Z → ionPath
    std::map<unsigned int, std::string> ionPaths;
	for (auto charge : UnderStudyCharge) {
		ionPaths[charge] = std::string(getIonPath(charge).Data());
	}

    /////////////////// MC Sample //////////////////////////
    manager.FillMcSampleFromFileInfo(chain,"IL1");

    /////////////////// LOOP ///////////////////////////////
    for (NAIA::Event &event : chain) {

        manager.FillMcTrees(chain,event,"IL1");

    }
    // WRITE
    manager.WriteMcTrees(ionPaths,outname.Data(),isMC);

    return 0;
}