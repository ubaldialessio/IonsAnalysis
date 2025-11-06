#ifndef STREAMUTILITY_H
#define STREAMUTILITY_H
#include "includes.h"

//To be replaced based on where the repository is installed
TString repoDir = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis";

namespace StreamUtility {
    std::map<TString , TString> InpDirectories { {"SelectCounts","/Counts/"}, 
                                              {"Livetime","/Livetime/"}, 
                                              {"RawAcceptance","/RawAcceptance/"}, 
                                              {"SelectPassedTree","/Passed/"},
                                              {"UnfoldAcceptance","/Unfolding/"},
                                              {"L1Eff","/Efficiencies/L1/"},
											  {"L1EffUnb","/Efficiencies/L1Unb/"},
											  {"L1ChargeCut","Efficiencies/L1ChargeCut/"},
                                              {"L9Eff","/Efficiencies/L9/"},
                                              {"TofEff","/Efficiencies/Tof/"},
                                              {"TrackEff","/Efficiencies/Track/"},
											  {"TrackChEff","/Efficiencies/TrackCharge/"},
                                              {"TriggEff","/Efficiencies/Trigger/"},
											  {"DaqEff","/Efficiencies/Daq/"},
                                              {"BuildFlux","/Flux/"},
											  {"BuildFluxForJack","/Flux/"},
											  {"RigidityMap","/RigidityMap/"},
											  {"VariableController","/../../ControlVariables/"},
											  {"AboveL1","/Fragments/"},
											  {"EvaluateFragments","/Fragments/Result/"},
											  {"SecondTrackVsFirst","/ControlVariables/SecondTrackVsFirst/"},
											  {"DataSelector","/dat/"}
                                              };
    struct rootFiles{
        std::pair<std::vector<TFile*>,std::vector<TFile*>> effFiles;
        TFile *dal1;
        TFile *datf;
        TFile *datr;
        TFile *datk;
		TFile *dal1u;
		TFile *datkch;
		TFile *daq;
		TFile *l1ch;
        TFile *mcl1;
        TFile *mctf;
        TFile *mctr;
        TFile *mctk;
		TFile *mcl1u;
		TFile *mctkch;
        TFile *livetime;
        TFile *rawacc;
        TFile *counts;
        TFile *passed;
    };
	TString formatOutnameOptions(unsigned int charge, TString outname, TString sec_track, TString l1_cut,
								 TString mcType, TString timePeriod);
    TString getOutputDir(unsigned int charge, TString exename, TString outname);
    TString getInputDir (unsigned int charge, TString exename, TString outname);
    rootFiles loadFile  (unsigned int charge, TString opt2,    TString mcType, TString sec_track);
    std::pair<std::vector<TFile*>,std::vector<TFile*>> getEfficiencyFiles(TString option, unsigned int charge);
}

TString StreamUtility::formatOutnameOptions(unsigned int charge, TString outname, TString sec_track, TString l1_cut, 
											TString mcType, TString timePeriod) {
	TString result, mc;
	TString SecTrack = sec_track == "y" ?  "_sec_track" : "";
	TString L1Cut    = l1_cut == "y" ? "_l1_cut" : "";
	if (mcType=="s") mc = "s";
	if (mcType=="g") mc == "global";
	if (mcType=="")  mc == "";
	result = Form("%s%s%s%s%s.root",outname.Data(),SecTrack.Data(),L1Cut.Data(),mc.Data(),timePeriod.Data());
	std::cout << result << std::endl;
	return result;
}
TString StreamUtility::getOutputDir(unsigned int charge, TString exename, TString outname) {
    TString executable = exename(exename.Last('/') + 1, exename.Length());
    std::cout << "Executable name: " << executable.Data() << std::endl;
    TString p = repoDir+"/IonsSelected/"+getIonPath(charge)+StreamUtility::InpDirectories[executable]+outname;
	std::cout << "Output: " << p << std::endl;
    return p;
}
TString StreamUtility::getInputDir(unsigned int charge, TString exename, TString outname) {
    TString p = "";
    return p;
}
std::pair<std::vector<TFile*>,std::vector<TFile*>> StreamUtility::getEfficiencyFiles(TString option, unsigned int charge) {
	std::pair<std::vector<TFile*>,std::vector<TFile*>> p = {};
	std::vector<TString> effDir = {"L1","Tof","Trigger","Track","L1Unb","TrackCharge"};
	TString ionPath=getIonPath(charge);
	TString ionName=getIonName(charge);
	if (option == "single") {
		TString path = repoDir+"/IonsSelected/"+ionPath+"/Efficiencies/";
		//Build the path for the 5 eff files
		std::vector<TString> da_eff, mc_eff;
		for (const auto i : effDir) {
			da_eff.push_back( Form("%s%s/%s_%s.root",path.Data(), i.Data(),"data",i.Data()) );
			mc_eff.push_back( Form("%s%s/%s_%s.root",path.Data(), i.Data(),"mc",i.Data()) );
		}
		//Build the root file for data and mc and insert inside p
		for (const auto i : da_eff) {
			auto d = new TFile(i.Data());
			p.first.push_back(d);
		}
		for (const auto i : mc_eff) {
			auto m = new TFile(i.Data());
			p.second.push_back(nullptr);
		}
		return p;
	}
	if (option == "all") {
		TString datPath = repoDir+"/IonsSelected/"+ionPath+"/dat/dat.root";
		auto d = new TFile(datPath.Data() );
		p.first.push_back(d);
		TString mcPath = repoDir+"/IonsSelected/"+ionPath+"/mc/mc.root";
		auto m = new TFile(mcPath.Data() );
		p.second.push_back(m);
	}
	return p;
}
StreamUtility::rootFiles StreamUtility::loadFile(unsigned int charge, TString opt2, TString mcType, TString sec_track) {
	rootFiles r;
	auto ipath = getIonPath(charge);
	TString inp_sec_track;
	sec_track == "y" ? inp_sec_track = "_sec_track_" : "";
	r.effFiles = getEfficiencyFiles(opt2,charge); //need to finish the implementation. I missed the retrieving of mc_pass_gen and mc_samp and mc_samp 
	if (opt2=="single") {
		r.livetime = new TFile( Form(repoDir+"/IonsSelected/%s/Livetime/livetime.root",ipath.Data()) );
		r.rawacc   = new TFile( Form(repoDir+"/IonsSelected/%s/Passed/%spassed.root",ipath.Data(),inp_sec_track.Data()) );
		r.counts   = new TFile( Form(repoDir+"/IonsSelected/%s/Counts/%scounts.root",ipath.Data(),inp_sec_track.Data()) );
		r.daq	   = new TFile( Form(repoDir+"/IonsSelected/%s/Efficiencies/Daq/daq.root",ipath.Data()) );
		r.l1ch     = new TFile( Form(repoDir+"/IonsSelected/%s/Efficiencies/L1ChargeCut/l1ChargeCut.root",ipath.Data()) );
		if (mcType=="g")
        	r.passed   = new TFile( Form(repoDir+"/IonsSelected/%s/Passed/%spassed_global.root",ipath.Data(),inp_sec_track.Data()));
		else
			r.passed   = new TFile( Form(repoDir+"/IonsSelected/%s/Passed/%spassed.root",ipath.Data(),inp_sec_track.Data()));
		r.dal1	   = r.effFiles.first[0];
		r.datf	   = r.effFiles.first[1];
		r.datr	   = r.effFiles.first[2];
		r.datk	   = r.effFiles.first[3];
		r.dal1u    = r.effFiles.first[4];
		r.datkch   = r.effFiles.first[5];
		/*r.mcl1     = r.effFiles.second[0];
		r.mctf     = r.effFiles.second[1];
		r.mctr     = r.effFiles.second[2];
		r.mctk     = r.effFiles.second[3];
		r.mcl1u    = r.effFiles.second[4];
		r.mctkch   = r.effFiles.second[5];*/
	}
	if (opt2=="all") {
		r.livetime = r.effFiles.first[0];
		r.rawacc   = r.effFiles.second[0];
		r.counts   = r.effFiles.first[0];
        r.passed   = new TFile( Form(repoDir+"/IonsSelected/%s/Passed/passed.root",ipath.Data()));
		r.dal1	   = r.effFiles.first[0];
		r.datf	   = r.effFiles.first[0];
		r.datr	   = r.effFiles.first[0];
		r.datk	   = r.effFiles.first[0];
		r.mcl1     = r.effFiles.second[0];
		r.mctf     = r.effFiles.second[0];
		r.mctr     = r.effFiles.second[0];
		r.mctk     = r.effFiles.second[0];
	}
	return r;
}
/*StreamUtility::rootFiles StreamUtility::loadFile(unsigned int charge, TString opt2, TString mcType, TString sec_track) {
    rootFiles r;
    auto ipath = getIonPath(charge);

    // opzionale: suffisso per sec_track
    TString inp_sec_track;
    if (sec_track == "y")
        inp_sec_track = "_sec_track_";

    // recupera i file efficienza
    r.effFiles = getEfficiencyFiles(opt2, charge);

    if (opt2 == "single") {
        // File principali
        r.livetime = new TFile(Form(repoDir+"/IonsSelected/%s/Livetime/livetime.root",
                                    ipath.Data()), "READ");
        r.rawacc   = new TFile(Form(repoDir+"/IonsSelected/%s/RawAcceptance/%srawacc.root",
                                    ipath.Data(), inp_sec_track.Data()), "READ");
        r.counts   = new TFile(Form(repoDir+"/IonsSelected/%s/Counts/%scounts.root",
                                    ipath.Data(), inp_sec_track.Data()), "READ");
        r.daq      = new TFile(Form(repoDir+"/IonsSelected/%s/Efficiencies/Daq/daq.root",
                                    ipath.Data()), "READ");

        if (mcType == "g")
            r.passed = new TFile(Form(repoDir+"/IonsSelected/%s/Passed/%spassed_global.root",
                                      ipath.Data(), inp_sec_track.Data()), "READ");
        else
            r.passed = new TFile(Form(repoDir+"/IonsSelected/%s/Passed/%spassed.root",
                                      ipath.Data(), inp_sec_track.Data()), "READ");

        // File efficienze (unico file con dentro IL1 e FS)
        r.effData = r.effFiles.first[0];
        r.effMC   = r.effFiles.second[0];
    }

    if (opt2 == "all") {
        // In questo caso usi dat.root e mc.root
        r.livetime = r.effFiles.first[0];   // dentro c'è già IL1/FS
        r.rawacc   = r.effFiles.second[0];
        r.counts   = r.effFiles.first[0];
        r.passed   = new TFile(Form(repoDir+"/IonsSelected/%s/Passed/passed.root",
                                    ipath.Data()), "READ");

        r.effData = r.effFiles.first[0];
        r.effMC   = r.effFiles.second[0];
    }

    return r;
}*/


void printInfo(unsigned int charge, TString exename, TString outname) {

}

#endif