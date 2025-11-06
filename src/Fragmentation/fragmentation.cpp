#include "definition.h"
#include "binning.h"
#include "utils.h"

int main(int argc, char *argv[]) {
	if (argc < 4) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> <SecTrack> \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString infilename=argv[2],
		    outname=argv[3],
		    out, ionPath=getIonName(charge), sec_track = argv[4];

	TH1::SetDefaultSumw2();

	// setup binning template useed for L1/L2 templates.
	const int nRigbins_LightIons = 24;
	const double Rigbins_LightIons[nRigbins_LightIons] = {0.8,  1.16, 1.51, 1.92, 2.40,  2.97, 3.64, 4.43,
	                                                  5.37, 6.47, 7.76, 9.26, 11.00, 13.0, 15.3, 18.0,
	                                                  21.1, 24.7, 28.8, 33.5, 38.9,  45.1, 52.2, 5000};
	
	const int nRigbins_HeavyIons = 21;
	const double Rigbins_HeavyIons[nRigbins_HeavyIons] = {0.8,  2.40, 2.97, 3.64, 4.43, 5.37, 6.47, 7.76, 9.26, 11.00, 13.0,
	                                                  15.3, 18.0, 21.1, 24.7, 28.8, 33.5, 38.9, 45.1, 52.2, 5000};

	const int nP_rigBins = 13; 
	const double P_rigBins[nP_rigBins] = {0.8,3.64,5.37,7.76,11.00,15.3,21.1,38.9,70.00,150.00,300.00,500.00,1000.00};												  
	TList *l1TemplateList   = new TList();
	TList *l2TemplateList   = new TList();
	TList *l1Distribution   = new TList();
	TList *l1KernelTreeList = new TList();
	TList *l2KernelTreeList = new TList();
	TList *l1DistKernTreeList=new TList();

	float L1charge = 0;
	float L2charge = 0;
  	double L1weight = 1;
	double L2weight = 1;
	double L1Dweight = 1;

	auto l1rigBinTemplate = new TH1D("l1rigBinTemplate", "R (GV)", (charge==15 || charge==16 || charge==17)? nRbins_HighZ - 1 : nRbins_HighZ-1,
																   (charge==15|| charge==16 || charge==17)? Rbins_HighZ : Rbins_HighZ);
	auto rigBinTemplate   = new TH1D("RigTemplateLightIons", "R (GV)", (charge==15|| charge==16 || charge==17)? nRbins_HighZ - 1 : nRbins_HighZ-1,
																       (charge==15|| charge==16 || charge==17)? Rbins_HighZ : Rbins_HighZ);

	//L1 template
	for (Int_t ibin = 0; ibin < l1rigBinTemplate->GetNbinsX(); ibin++) {
	    l1TemplateList->Add(new TH1D(Form("L1Template_%03i", ibin),
	                               Form("%5.3f < R (GV) < %5.3f;Q_{L1X};Counts", l1rigBinTemplate->GetBinLowEdge(ibin + 1),
	                                    l1rigBinTemplate->GetBinLowEdge(ibin + 2)),
	                               600, 0, 32));
		l1KernelTreeList->Add(new TTree(Form("L1KTree_%03i", ibin), Form("%5.3f < R (GV) < %5.3f", l1rigBinTemplate->GetBinLowEdge(ibin + 1),
                                        l1rigBinTemplate->GetBinLowEdge(ibin + 2))));
		static_cast<TTree *>(l1KernelTreeList->At(ibin))->SetDirectory(0);
    	static_cast<TTree *>(l1KernelTreeList->At(ibin))->Branch("L1charge", &(L1charge), "L1charge/F");
    	static_cast<TTree *>(l1KernelTreeList->At(ibin))->Branch("L1weight", &L1weight, "L1weight/D");
	}
	
	//L2 template
	for (Int_t ibin = 0; ibin < rigBinTemplate->GetNbinsX(); ibin++) {
		l2TemplateList->Add(new TH1D(Form("L2Template_%03i", ibin),
                               Form("%5.3f < R (GV) < %5.3f;Q_{L2X};Counts", rigBinTemplate->GetBinLowEdge(ibin + 1),
                                        rigBinTemplate->GetBinLowEdge(ibin + 2)),
                               		600, 0, 32));
		l2KernelTreeList->Add(new TTree(Form("L2KTree_%03i", ibin), Form("%5.3f < R (GV) < %5.3f", rigBinTemplate->GetBinLowEdge(ibin + 1),
                                        rigBinTemplate->GetBinLowEdge(ibin + 2))));
		static_cast<TTree *>(l2KernelTreeList->At(ibin))->SetDirectory(0);
    	static_cast<TTree *>(l2KernelTreeList->At(ibin))->Branch("L2charge", &(L2charge), "L2charge/F");
    	static_cast<TTree *>(l2KernelTreeList->At(ibin))->Branch("L2weight", &L2weight, "L2weight/D");
	}   

	//L1 charge distribution
	for (Int_t ibin = 0; ibin < l1rigBinTemplate->GetNbinsX(); ibin++) {
		l1Distribution->Add(new TH1D(Form("L1Distribution_%03i", ibin),
                               Form("%5.3f < R (GV) < %5.3f;Q_{L1};Counts", l1rigBinTemplate->GetBinLowEdge(ibin + 1),
                                        l1rigBinTemplate->GetBinLowEdge(ibin + 2)),
                               		600, 0, 32));
		l1DistKernTreeList->Add(new TTree(Form("L1DKTree_%03i", ibin), Form("%5.3f < R (GV) < %5.3f", l1rigBinTemplate->GetBinLowEdge(ibin + 1),
                                        l1rigBinTemplate->GetBinLowEdge(ibin + 2))));
		static_cast<TTree *>(l1DistKernTreeList->At(ibin))->SetDirectory(0);
    	static_cast<TTree *>(l1DistKernTreeList->At(ibin))->Branch("L1Dcharge", &(L1charge), "L1Dcharge/F");
    	static_cast<TTree *>(l1DistKernTreeList->At(ibin))->Branch("L1Dweight", &L1Dweight, "L1Dweight/D");
	}               



	// this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

	//check if is a valid input
    bool validInput = true;
	//process
	if (validInput) {
		NAIA::NAIAChain chain;
		if(infilename.Contains(".root") /*&& filesystem::exists(infilename.Data())*/ ){
		    chain.Add(infilename.Data());
		}else if (infilename.Contains(".txt") /*&& filesystem::exists(infilename.Data())*/ ){
		    ifstream infilelist(infilename.Data());
		    TString bufname;
		    while(infilelist >> bufname) 
		      if(bufname.Contains(".root") /*&& filesystem::exists(bufname.Data())*/ )
		        chain.Add(bufname.Data());
		}
		chain.SetupBranches();
		out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"+ionPath+"/"+outname;
		if (sec_track=="y") out="/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"+ionPath+"/sec_track_"+outname;
		unsigned int utime, oldtime;
		float reco_il1, reco_inn, beta, lt, l1charge, l2charge;
		double gen, cutoff;
		double Amass = 1;
		double weight = 1;
		Amass = GetDataMass(charge); // Getting atomic mass number
		
auto mysel = Mysel(charge);
auto myselNoL1 = MyselNoL1(charge);
auto l1temp =
    BuildTemplatesSel::InnerTrackerChargeInRange(static_cast<float>(charge), CRT) &&
	BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge), UTC) &&       
	BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge), LTC) &&
	ns::TrackerLayer::ChargeStatus(1);
	//ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT);
auto l2temp = 
    BuildTemplatesSel::TrackerLayerChargeInRange(1, static_cast<float>(charge), CRT) &&
	BuildTemplatesSel::InnerTrackerChargeInRange(static_cast<float>(charge), CRT) &&
    BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge),UTC) &&
    BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge),LTC) &&
	ns::TrackerLayer::ChargeStatus(1) &&
	ns::TrackerLayer::ChargeStatus(2);
    //BuildTemplatesSel::L3L8Charge(static_cast<float>(charge), CRT) &&				    //L3-L8
	//ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT);
auto hasGoodSecondTrack = GoodSecTrack();
auto secTrackOnDiagonal = MySel::IsSecondTrackOnDiagonal(FIT,IL1,0.3);
auto l1ClusterCut = MySel::L1ClusterCutBetween(10,30,YSD,ALL_LAYER);

/////////////////// LOOP ///////////////////////////////
		for(Event& event : chain) {
			utime = event.header->UTCTime;
				if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
					continue;
				if (isRunBad(utime)==true) //skipping trigger study and wrong configuration
					continue;
				if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
					continue;
				cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
				if (l1temp(event)) {
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
					auto Rigidity_IL1 = event.trTrackBase->Rigidity[FIT][IL1];
					int bin = l1rigBinTemplate->FindBin(Rigidity_IL1);
					double rlowedge = hist_rig_highZ->GetBinLowEdge(bin);
					if (bin < 1 || bin > l1rigBinTemplate->GetNbinsX())
						continue;
					if (rlowedge < 1.2*cutoff) continue;
					l1charge = event.trTrackBase->LayerChargeXY[0][CRT];
					static_cast<TH1D *>(l1TemplateList->At(bin - 1))->Fill(l1charge,L1weight);
					static_cast<TTree *>(l1KernelTreeList->At(bin - 1))->Fill();
				}
				if (l2temp(event)) {
					auto Rigidity_IL1 = event.trTrackBase->Rigidity[FIT][IL1];
					int bin = rigBinTemplate->FindBin(Rigidity_IL1);
					double rlowedge = hist_rig_highZ->GetBinLowEdge(bin);
					if (bin < 1 || bin > rigBinTemplate->GetNbinsX())
						continue; 
					if (rlowedge < 1.2*cutoff) continue;
					l2charge = event.trTrackBase->LayerChargeXY[1][CRT];
				    static_cast<TH1D *>(l2TemplateList->At(bin - 1))->Fill(l2charge,L2weight);
					static_cast<TTree *>(l2KernelTreeList->At(bin - 1))->Fill();
				}
				if (myselNoL1(event)) {
					auto Rigidity_IL1 = event.trTrackBase->Rigidity[FIT][IL1];
					int bin = l1rigBinTemplate->FindBin(Rigidity_IL1);
					double rlowedge = hist_rig_highZ->GetBinLowEdge(bin);
					if (bin < 1 || bin > l1rigBinTemplate->GetNbinsX())
						continue;
					if (rlowedge < 1.2*cutoff) continue;
					auto L1charge = event.trTrackBase->LayerChargeXY[0][CRT];
					static_cast<TH1D *>(l1Distribution->At(bin - 1))->Fill(L1charge,L1weight);
					static_cast<TTree *>(l1DistKernTreeList->At(bin - 1))->Fill();					
				}
		} //event loop	
	//WRITE//
	auto outfile = new TFile(out, "recreate");
	outfile->WriteTObject(l1TemplateList, Form("L1TemplateList_%i", charge), "Overwrite");
	outfile->WriteTObject(l2TemplateList, Form("L2TemplateList_%i", charge), "Overwrite");
	outfile->WriteTObject(l1Distribution, Form("L1Distribution_%i", charge), "Overwrite");
	outfile->WriteTObject(l1KernelTreeList, Form("l1KernelTreeList_%i", charge), "Overwrite");
	outfile->WriteTObject(l2KernelTreeList, Form("l2KernelTreeList_%i", charge), "Overwrite");
	outfile->WriteTObject(l1DistKernTreeList, Form("l1DistKernTreeList_%i", charge), "Overwrite");

	outfile->Close();
	} //valid input
} //main

