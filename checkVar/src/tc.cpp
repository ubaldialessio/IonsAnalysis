#include "analysis.h"
using namespace std;
#define P1 1
#define P2 2
#define P3 4
#define P4 8
#define P5 16
#define P6 32
#define P7 64
#define P8 128
#define P9 256
#define P10 512
#define P11 1024
#define P12 2048
#define P13 4096
#define P14 8192
#define P15 16384
#define P16 32768

/*
bool IsInsideucial(NAIA::TrTrack::FitPositionHeight layer, NAIA::Event &event) {
  constexpr float fid_R[] = {62, 62, 46, 46, 46, 46, 46, 46, 43}; // Beware, for L9 the check is X only, not R
  constexpr float fid_Y[] = {47, 40, 44, 44, 36, 36, 44, 44, 29};
  if(!NAIA::ContainsKeys(event.tofBaseSt->InterpolateAtZ(NAIA::TrTrack::fitPositionHeightZ[layer]))) return false;
  double x = event.tofBaseSt->InterpolateAtZ(NAIA::TrTrack::fitPositionHeightZ[layer])[0];
  double y = event.tofBaseSt->InterpolateAtZ(NAIA::TrTrack::fitPositionHeightZ[layer])[1];
  if(layer == NAIA::TrTrack::FitPositionHeight::Layer9){
    if (fabs(x) > fid_R[layer] || fabs(y) > fid_Y[layer]) return false;
  }
  else{
	    if(sqrt(x * x + y * y) > fid_R[layer] || fabs(y) > fid_Y[layer]) return false;
  }
  return true;
}

class tInnerucialVolume : public NSL::Selection{
	public:
  tInnerucialVolume() {
	m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
	  int nhits = 0;
      if(IsInsideucial(L2, event)) nhits++;
      if(IsInsideucial(L3, event)) nhits++;
      if(IsInsideucial(L4, event)) nhits++;
      if(IsInsideucial(L5, event)) nhits++;
      if(IsInsideucial(L6, event)) nhits++;
      if(IsInsideucial(L7, event)) nhits++;
      if(IsInsideucial(L8, event)) nhits++;
      if(nhits < 5) return false;
      if(!IsInsideucial(L2, event))  return false;
      if(!(IsInsideucial(L3, event) || IsInsideucial(L4, event))) return false;
      if(!(IsInsideucial(L5, event) || IsInsideucial(L6, event))) return false;
      if(!(IsInsideucial(L7, event) || IsInsideucial(L8, event))) return false;
      return true;
    });
  }
};

class tL1ucialVolume : public NSL::Selection{
	public:
  tL1ucialVolume(){
	m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
	  return IsInsideucial(L1, event);
    });
  }
};
*/
int main(int argc, char **argv){

  NAIA::NAIAChain chain;
  TString outfilename = GetOutname_FillChain(chain, argc, argv);
  chain.SetupBranches();
  bool isMC = chain.IsMC();

//eff track

  float lt, zenith;
  float gen, utof_edep, ltof_edep;
  float charge_tof, charge_ltof, beta;
  unsigned int utime;
  unsigned short int cutmask, npart, pid;
  TTree *tree=new TTree();
  tree->Branch("utime", &utime);
  tree->Branch("charge_tof", &charge_tof);
  tree->Branch("charge_ltof", &charge_ltof);
  tree->Branch("beta", &beta);
  tree->Branch("npart", &npart);
  tree->Branch("cutmask", &cutmask);
  tree->Branch("utof_edep", &utof_edep);
  tree->Branch("ltof_edep", &ltof_edep);
  if(isMC){
    tree->Branch("gen", &gen);
    tree->Branch("pid", &pid);
  }

  auto den  =
    ms::Trig::PhysicsTrigger()                                             &&
    ms::TofSt::tBetaInRange(0.3, std::numeric_limits<float>::max(), BTH)   &&
    ms::TofSt::tChi2TimeLessThan(2)                                        &&
    ms::TofSt::tChi2CooLessThan(2)                                         &&
    ms::TofSt::tChargeInRange(0.7, 1.3, UTC)                               &&
    ms::TofSt::tChargeInRange(0.7, std::numeric_limits<float>::max(), LTC) &&
    ms::TofSt::tInnerFiducialVolume()                                      &&
    ms::TofSt::tL1FiducialVolume()                                         &&
    ms::ExtChargeInRange(0.6, 1.9, EL1, CRT);
  auto num =
    ns::Track::ChiSquareLessThan(10, YSD, FIT, INN)                        &&
    ns::InnerTracker::ChargeInRange(0.7, 1.5, CRT)                         &&
    ns::InnerTracker::ChargeRMSLessThan(0.4, CRT)                          &&
    ns::InnerTracker::HitPattern()                                         &&
    ns::Track::InnerFiducialVolume(FIT, INN)                               &&
    ns::InnerTracker::NHitsGreaterThan(4, YSD);

//eff_track
/*
  auto trig = ms::Trig::PhysicsTrigger();
  auto binr = ms::TofSt::tBetaInRange(0.3, std::numeric_limits<float>::max(), BTH);
  auto ch2t = ms::TofSt::tChi2TimeLessThan(2);
  auto ch2c = ms::TofSt::tChi2CooLessThan(2);
  auto zToF = ms::TofSt::tChargeInRange(0.7, 1.3, UTC);
  auto zLoF = ms::TofSt::tChargeInRange(0.7, std::numeric_limits<float>::max(), LTC);
  auto tinf = tInnerucialVolume();
  auto tl1f = tL1ucialVolume();
  auto exch = ms::ExtChargeInRange(0.6, 1.9, EL1, CRT);

  auto c2IN = ns::Track::ChiSquareLessThan(10, YSD, FIT, INN);
  auto zINN = ns::InnerTracker::ChargeInRange(0.7, 1.5, CRT);
  auto zRMS = ns::InnerTracker::ChargeRMSLessThan(0.4, CRT);
  auto hitp = ns::InnerTracker::HitPattern();
  auto hitn = ns::InnerTracker::NHitsGreaterThan(4, YSD);
  auto Ifid = ns::Track::InnerFiducialVolume(FIT, INN);
*/
// select_events
/*
  float lt, zenith;
  float gen, utof_edep, ltof_edep, reco_inn, reco, c2INN, c2IL1, cutoff;
  float charge_inn, rms, charge_l1, charge_tof, charge_ltof, beta;
  unsigned int utime;
  unsigned short int cutmask, npart, pid;
  TTree *tree=new TTree();
  tree->Branch("utime", &utime);
  tree->Branch("reco", &reco);
  tree->Branch("reco_inn", &reco_inn);
  tree->Branch("c2INN", &c2INN);
  tree->Branch("c2IL1", &c2IL1);
  tree->Branch("charge_inn", &charge_inn);
  tree->Branch("charge_l1", &charge_l1);
  tree->Branch("charge_tof", &charge_tof);
  tree->Branch("charge_ltof", &charge_ltof);
  tree->Branch("beta", &beta);
  tree->Branch("rms", &rms);
  tree->Branch("cutoff", &cutoff);
  tree->Branch("npart", &npart);
  tree->Branch("cutmask", &cutmask);
  tree->Branch("utof_edep", &utof_edep);
  tree->Branch("ltof_edep", &ltof_edep);
  if(isMC){
    tree->Branch("gen", &gen);
    tree->Branch("pid", &pid);
  }

  auto trig = ms::Trig::PhysicsTrigger();
  auto binr = ns::Tof::BetaInRange(0.3, std::numeric_limits<float>::max(), BTH);
  auto zToF = ns::Tof::ChargeInRange(0.5, 2.5, UTC);
  auto c2IN = ns::Track::ChiSquareLessThan(10, YSD, FIT, INN);
  auto hitp = ns::InnerTracker::HitPattern();
  auto hitn = ns::InnerTracker::NHitsGreaterThan(4, YSD);
  auto c2L1 = ns::Track::ChiSquareLessThan(10, YSD, FIT, IL1);
  auto zINN = ns::InnerTracker::ChargeInRange(0.7, 1.5, CRT);
  auto zRMS = ns::InnerTracker::ChargeRMSLessThan(0.4, CRT);
  auto zIL1 = ns::TrackerLayer::ChargeInRange(1, 0.6, 1.9, CRT);
  auto Ifid = ns::Track::InnerFiducialVolume(FIT, IL1);
  auto Lfid = ns::Track::L1FiducialVolume(FIT, IL1);
  auto Lres = ns::Track::L1NormResidualLessThan(10, FIT);
  auto Mass = ms::MassCut(BTH, FIT, INN);
  auto Part = ms::NPartCut();
*/

  for(NAIA::Event &event : chain){

    cutmask = 0;
/*
    if(!event.trTrackBase->FitIDExists(FIT, INN) || !event.trTrackBase->FitIDExists(FIT, IL1))
      continue;
    reco_inn   = event.trTrackBase->Rigidity[FIT][INN];
    reco       = event.trTrackBase->Rigidity[FIT][IL1];
    charge_inn = event.trTrackBase->InnerCharge[CRT];
    charge_l1  = event.trTrackBase->LayerChargeXY[0][CRT];
    c2INN      = event.trTrackBase->TrChiSq[FIT][INN][YSD];
    c2IL1      = event.trTrackBase->TrChiSq[FIT][IL1][YSD];
    rms        = event.trTrackBase->InnerChargeRMS[CRT];
*/


    utime      = event.header->UTCTime;
    npart      = event.evSummary->NParticle;
    charge_tof = event.tofBase->Charge[UTC];
    charge_ltof = event.tofBase->Charge[LTC];
    beta       = event.tofBase->Beta[BTH];
    utof_edep  = event.tofPlus->LayerEdep[0]+event.tofPlus->LayerEdep[1];
    ltof_edep  = event.tofPlus->LayerEdep[2]+event.tofPlus->LayerEdep[3];

    if(isMC){
      gen = event.mcTruthBase->Primary.GetGenMomentum();
      pid = 0;
      for(int i=0; i<9; i++){
        short int parpid = event.mcTruthPlus->TrackMCHits[i].pID;
        if(!parpid) continue;
        if(pid!=parpid) pid = parpid;
      }
    }
    else{
//      cutoff = chain.GetEventRTIInfo().MaxIGRFCutoff[0][1];
      lt = chain.GetEventRTIInfo().LivetimeFraction;
      zenith = chain.GetEventRTIInfo().Zenith;
      if (lt < 0.05 || zenith > 25
      || chain.GetEventRTIInfo().MeanAlignDiffExtLayer[0][1] > 35
      || chain.GetEventRTIInfo().IsInSAA()) cutmask|=P1;
    }

    if(!den(event)) cutmask|=P2;
    if(!num(event)) cutmask|=P3;

//eff_track
/*
    if(!trig(event)) cutmask|=P2;
    if(!binr(event)) cutmask|=P3;
    if(!ch2t(event)) cutmask|=P4;
    if(!ch2c(event)) cutmask|=P5;
    if(!zToF(event)) cutmask|=P6;
    if(!zLoF(event)) cutmask|=P7;
    if(!tinf(event)) cutmask|=P8;
    if(!tl1f(event)) cutmask|=P9;
    if(!exch(event)) cutmask|=P10;

    if(!c2IN(event)) cutmask|=P11;
    if(!zINN(event)) cutmask|=P12;
    if(!zRMS(event)) cutmask|=P13;
    if(!hitp(event)) cutmask|=P14;
    if(!hitn(event)) cutmask|=P15;
    if(!Ifid(event)) cutmask|=P16;
*/
//select_events
/*
    if(!trig(event)) cutmask|=P2;
    if(!binr(event)) cutmask|=P3;
    if(!zToF(event)) cutmask|=P4;
    if(!c2IN(event)) cutmask|=P5;
    if(!hitp(event)) cutmask|=P6;
    if(!hitn(event)) cutmask|=P7;
    if(!c2L1(event)) cutmask|=P8;
    if(!zINN(event)) cutmask|=P9;
    if(!zRMS(event)) cutmask|=P10;
    if(!zIL1(event)) cutmask|=P11;
    if(!Ifid(event)) cutmask|=P12;
    if(!Lfid(event)) cutmask|=P13;
    if(!Lres(event)) cutmask|=P14;
    if(!Mass(event)) cutmask|=P15;

    if(!isMC){
      if(hist_rig->GetBinLowEdge(hist_rig->FindBin(reco)) < cutoff)
        cutmask|=P16;
    }
    else{
      if(reco<=0) cutmask|=P16;
    }
*/

    tree->Fill();

  }

  auto outfile = new TFile(outfilename, "recreate");
  outfile->WriteTObject(tree, "reco_gen");
  outfile->Close();
}
