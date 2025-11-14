#ifndef DEFINITION_H
#define DEFINITION_H
#include "includes.h"

using namespace std;
namespace ns = NSL::Selections;
int fov = 1; // Fields of view: 0 = 25 deg, 1 = 30 deg, 2 = 35 deg, 3 = 40 deg
float inf = 9999999;

struct Sample { //ORIGINALE
  Short_t passed;
  Float_t reco_inn;
  Float_t reco_il1;
  Float_t reco_beta;
  Int_t   species;
};

auto CRT = NAIA::TrTrack::ChargeRecoType::YJ;
auto STD = NAIA::TrTrack::ChargeRecoType::STD;
auto FIT = NAIA::TrTrack::Fit::GBL;
auto INN = NAIA::TrTrack::Span::InnerOnly;
auto IL1 = NAIA::TrTrack::Span::InnerL1;
auto IL9 = NAIA::TrTrack::Span::InnerL9;
auto YSD = NAIA::TrTrack::Side::Y;
auto XSD = NAIA::TrTrack::Side::X;
auto UTC = NAIA::Tof::ChargeType::Upper;
auto LTC = NAIA::Tof::ChargeType::Lower;
auto TOT = NAIA::Tof::ChargeType::Total;
auto BTH = NAIA::Tof::BetaType::BetaH;
auto ON  = NAIA::Tof::OnTime;
auto OFF = NAIA::Tof::OnTime;
auto L1  = NAIA::TrTrack::FitPositionHeight::Layer1;
auto L2  = NAIA::TrTrack::FitPositionHeight::Layer2;
auto L3  = NAIA::TrTrack::FitPositionHeight::Layer3;
auto L4  = NAIA::TrTrack::FitPositionHeight::Layer4;
auto L5  = NAIA::TrTrack::FitPositionHeight::Layer5;
auto L6  = NAIA::TrTrack::FitPositionHeight::Layer6;
auto L7  = NAIA::TrTrack::FitPositionHeight::Layer7;
auto L8  = NAIA::TrTrack::FitPositionHeight::Layer8;
auto EL1 = NAIA::UnbExtHitBaseData::ExtHit::L1;
auto EL9 = NAIA::UnbExtHitBaseData::ExtHit::L9;
auto TRD_UP = NAIA::TrdK::ChargeType::Upper;
auto TRD_LOW = NAIA::TrdK::ChargeType::Lower;
auto TRD_TOT = NAIA::TrdK::ChargeType::Total;
auto TRD_ONTRACK = NAIA::TrdK::QualType::OnTrack;
auto TRD_OFFTRACK = NAIA::TrdK::QualType::OffTrack;
auto ONE_MM = NAIA::TrTrack::DistanceFromTrack::Onemm;
auto ONE_CM = NAIA::TrTrack::DistanceFromTrack::Onecm;
auto TWO_CM = NAIA::TrTrack::DistanceFromTrack::Twocm;
auto FIVE_CM = NAIA::TrTrack::DistanceFromTrack::Fivecm;
auto TEN_CM = NAIA::TrTrack::DistanceFromTrack::Tencm;
auto ALL_LAYER = NAIA::TrTrack::DistanceFromTrack::AllLayer;



constexpr double fid_R[] = {62, 62, 46, 46, 46, 46, 46, 46};
constexpr double fid_Y[] = {47, 40, 44, 44, 36, 36, 44, 44};

short bitcounter(short n) {
  if(n==0) return 0;
  else return (n&1)+bitcounter(n>>1);
}

namespace trig 	  = Efficiency::TriggerEff;
namespace track        = Efficiency::TrTrackEffSel;
namespace InnerTracker = Efficiency::InnerTracker;
namespace TofSt        = Efficiency::TofStandalone;
namespace MySel 	  = Efficiency::MySel;
namespace trd          = Efficiency::trd;

auto GoodSecTrack() {
  return MySel::SecondTrackYHitsGreaterThan(FIT,IL1,2) &&
	  MySel::SecondTrkChiSquareLessThan(20.0, YSD, FIT, IL1);	
}
auto MyselStoP(unsigned int charge) {
  return ns::Trigger::HasPhysicsTrigger() && 
        ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge-1.),UTC) && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) && 
        ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge-1.),CRT) && 
        InnerTracker::L1NormResidualLessThan(10.0,FIT) &&
        ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT) && 
        ns::Common::TrackerLayer::ChargeInRange(1,static_cast<float>(charge),CRT) &&
	 ns::TrackerLayer::ChargeStatus(1) && 
        ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, IL1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) && 
        ns::Track::HitCut(1) && 
        ns::Track::L1FiducialVolume(FIT, IL1) && 
	 ns::Track::InnerFiducialVolume(FIT,IL1,false); 
}
auto Mysel(unsigned int charge) {
  return ns::Trigger::HasPhysicsTrigger() && 
        ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),UTC) && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) && 
        ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT) && 
        InnerTracker::L1NormResidualLessThan(10.0,FIT) &&
        ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT) && 
        ns::Common::TrackerLayer::ChargeInRange(1,static_cast<float>(charge),CRT) &&
	 ns::TrackerLayer::ChargeStatus(1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) && 
        ns::Track::HitCut(1) && 
        ns::Track::L1FiducialVolume(FIT, IL1) && 
	 ns::Track::InnerFiducialVolume(FIT,IL1,false); 
}
auto MyselNoL1(unsigned int charge) {
  return ns::Trigger::HasPhysicsTrigger() && 
        ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),UTC) && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) && 
        ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT) && 
        InnerTracker::L1NormResidualLessThan(10.0,FIT) &&
        ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT) && 
        ns::Common::TrackerLayer::ChargeInRange(1,static_cast<float>(charge),CRT,true) &&
	 ns::TrackerLayer::ChargeStatus(1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) && 
        ns::Track::HitCut(1) && 
        ns::Track::L1FiducialVolume(FIT, IL1) && 
	 ns::Track::InnerFiducialVolume(FIT,IL1,false); 
}
auto Num_tof(unsigned int charge) {
 return ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),UTC); 
}
auto Den_tof(unsigned int charge) {
 return ns::Trigger::HasPhysicsTrigger() && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) && 
        ns::InnerTracker::ChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+0.5,CRT) && 
        InnerTracker::L1NormResidualLessThan(10.0,FIT) &&
	 ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT) && 
        ns::TrackerLayer::ChargeInRange(1,static_cast<float>(charge)-0.5,static_cast<float>(charge)+0.5,CRT) &&
	 ns::TrackerLayer::ChargeStatus(1) && 
        ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, IL1) && 
        ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN) && 
        ns::Track::HitCut(1) && 
        ns::Track::L1FiducialVolume(FIT, IL1) && 
        ns::Track::InnerFiducialVolume(FIT,IL1,false);
}
auto Den_tof_ltof(unsigned int  charge) {
 return Den_tof(charge) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4f,static_cast<float>(charge)+0.4f,LTC);
}
auto Den_tof_ltof_l9_05(unsigned int  charge) {
 return Den_tof(charge) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),LTC) &&
        track::ExtChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+100.,EL9,CRT);
}
auto Den_tof_ltof_l9_1(unsigned int  charge) {
 return Den_tof(charge) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),LTC) &&
        track::ExtChargeInRange(static_cast<float>(charge)-0.,static_cast<float>(charge)+100.,EL9,CRT);
}
auto Num_l1(unsigned int charge) {
 return ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT) && 
        ns::TrackerLayer::ChargeStatus(1) &&
        ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, IL1) && 
        ns::Track::L1FiducialVolume(FIT, IL1) && 
        ns::Track::InnerFiducialVolume(FIT, IL1) &&
        ns::Track::HitCut(1) && 
        InnerTracker::L1NormResidualLessThan(10.0,FIT); 
}
auto Den_l1(unsigned int charge) {
 return ns::Trigger::HasPhysicsTrigger() && 
        ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),UTC) && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
        ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT) && 
        ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, IL1) &&
        ns::Track::L1FiducialVolume(FIT, IL1) && 
        ns::Track::InnerFiducialVolume(FIT, IL1) &&
        track::IsHitPresent(EL1) &&
        track::UnbExtHitChargeStatus(EL1) &&
        track::UnbExtHitChargeInRange(EL1,static_cast<float>(charge),CRT,"qi");
}
auto Den_l1_ltof(unsigned int charge) {
 return Den_l1(charge) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC);
}
auto Den_l1_ltof_l9_05(unsigned int charge) {
 return Den_l1(charge) &&
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC) &&
        track::ExtChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+100.,EL9,CRT);
}
auto Den_l1_ltof_l9_1(unsigned int charge) {
 return Den_l1(charge) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC) &&
        track::ExtChargeInRange(static_cast<float>(charge)-0.,static_cast<float>(charge)+100.,EL9,CRT);
}
auto Num_l1Unb(unsigned int charge) {
 return track::IsHitPresent(EL1) &&
        track::UnbExtHitChargeStatus(EL1) &&
        track::UnbExtHitChargeInRange(EL1,static_cast<float>(charge),CRT,"UNB");
}
auto Den_l1Unb(unsigned int charge) {
 return ns::Trigger::HasPhysicsTrigger() && 
        ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),UTC) && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
        ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT) && 
        ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN) && 
        ns::Track::L1FiducialVolume(FIT, INN) && 
        ns::Track::InnerFiducialVolume(FIT, INN);
}
auto Den_l1Unb_ltof(unsigned int charge) {
 return Den_l1Unb(charge) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC);
}
auto Den_l1Unb_ltof_l9_05(unsigned int charge) {
 return Den_l1Unb(charge) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC) &&
        track::ExtChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+100.,EL9,CRT);
}
auto Den_l1Unb_ltof_l9_1(unsigned int charge) {
 return Den_l1Unb(charge) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC) &&
        track::ExtChargeInRange(static_cast<float>(charge)-0.,static_cast<float>(charge)+100.,EL9,CRT);
}
auto Num_track(unsigned int charge) {
 return ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) && 
        ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN); 
}
auto Den_track(unsigned int charge) {
 return ns::Trigger::HasPhysicsTrigger() && 
        TofSt::tofBetaInRange(0.4,inf,BTH) && 
        track::IsHitPresent(EL1) && 
        track::UnbExtHitChargeStatus(EL1) &&
        track::UnbExtHitChargeInRange(EL1,static_cast<float>(charge),CRT,"qi") &&
        track::InnerFiducialVolume(false) &&
        track::L1FiducialVolume() &&
        TofSt::tofChi2TimeLessThan(20) && 
        TofSt::tofChi2CooLessThan(20) && 
        TofSt::tofChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+0.5, UTC) && 
        TofSt::tofChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+0.5, LTC) &&
        TofSt::tofChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+0.5, TOT);
}
auto Num_trackCh(unsigned int charge) {
 return ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT);
}
auto Den_trackCh(unsigned int charge) {
 return ns::Trigger::HasPhysicsTrigger() && 
        ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) && 
        InnerTracker::L1NormResidualLessThan(10.0,FIT) &&
        ns::Track::L1FiducialVolume(FIT, IL1) && 
	 ns::Track::InnerFiducialVolume(FIT,IL1,false) &&
        ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT) && 
        ns::Common::TrackerLayer::ChargeInRange(1,static_cast<float>(charge),CRT) &&
        ns::TrackerLayer::ChargeStatus(1) &&
        ns::Track::HitCut(1) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, UTC);
}
auto Den_trackCh_ltof(unsigned int charge) {
 return Den_trackCh(charge) &&
	 ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC);
}
auto Den_trackCh_ltof_l9_05(unsigned int charge) {
 return Den_trackCh(charge) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC) &&
	 track::ExtChargeInRange(static_cast<float>(charge)-0.5,static_cast<float>(charge)+100.,EL9,CRT);
}
auto Den_trackCh_ltof_l9_1(unsigned int charge) {
 return Den_trackCh(charge) && 
        ns::Tof::ChargeInRange(static_cast<float>(charge)-0.4,static_cast<float>(charge)+0.4, LTC) &&
	 track::ExtChargeInRange(static_cast<float>(charge)-0.,static_cast<float>(charge)+100.,EL9,CRT);
}
auto Den_trig(unsigned int charge) {
 return ns::Tof::BetaInRange(0.4,inf,BTH) && 
        ns::Common::Tof::ChargeInRange(static_cast<float>(charge),UTC) && 
        ns::InnerTracker::HitPattern() && 
        ns::InnerTracker::NHitsGreaterThan(4, YSD) && 
        ns::Common::InnerTracker::ChargeInRange(static_cast<float>(charge),CRT) && 
        InnerTracker::L1NormResidualLessThan(10.0,FIT) &&
	 ns::Common::TrackerLayer::ChargeAsymmetry(1,CRT) && 
        ns::Common::TrackerLayer::ChargeInRange(1,static_cast<float>(charge),CRT) &&
	 ns::TrackerLayer::ChargeStatus(1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) && 
        ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) && 
        ns::Track::HitCut(1) && 
        ns::Track::L1FiducialVolume(FIT, IL1) && 
        ns::Track::InnerFiducialVolume(FIT,IL1,false); 
}
auto L1_temp(unsigned int charge) {
 return BuildTemplatesSel::InnerTrackerChargeInRange(static_cast<float>(charge), CRT) &&
	 BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge), UTC) &&       
	 BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge), LTC) &&
	 ns::TrackerLayer::ChargeStatus(1);
}
auto L2_temp(unsigned int charge) {
 return BuildTemplatesSel::TrackerLayerChargeInRange(1, static_cast<float>(charge), CRT) &&
	 BuildTemplatesSel::L3toL8ChargeSelection(static_cast<float>(charge), CRT) &&
        BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge),UTC) &&
        BuildTemplatesSel::TofChargeInRange(static_cast<float>(charge),LTC) &&
	 ns::TrackerLayer::ChargeStatus(1) &&
	 ns::TrackerLayer::ChargeStatus(2);
}

#endif
