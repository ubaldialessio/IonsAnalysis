#include "binning.h"

namespace ef = Efficiency::TriggerEffSel;

//selection
auto mysel =
	ns::Trigger::HasPhysicsTrigger() &&
	ns::Tof::BetaInRange(0.0,inf,BTH) &&
	ns::Tof::ChargeInRange(6.0-0.75,6.0+0.75,UTC) &&
    ns::InnerTracker::HitPattern() &&
    ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
    ns::InnerTracker::ChargeInRange(6.0-0.3,6.0+0.7,CRT) &&
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) &&
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) &&
    ns::Track::HitCut(1) &&
    ns::TrackerLayer::ChargeStatus(1) &&
	ns::TrackerLayer::ChargeInRange(1,0,6.0+0.8,CRT) &&
	ns::Track::InnerFiducialVolume(FIT, IL1) &&
	ns::Track::L1FiducialVolume(FIT, IL1)&&
	innTr::ChResoInnLessThan(CRT,0.2);

//tof eff
auto num_tof =
	ns::Tof::BetaInRange(0.0,inf,BTH) &&
	ns::Tof::ChargeInRange(6.0-0.75,6.0+0.75,UTC);
auto den_tof = 
	ns::Trigger::HasPhysicsTrigger() &&
	ns::InnerTracker::HitPattern() &&
	ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
	ns::InnerTracker::ChargeInRange(6.0-0.25,6.0+0.7,CRT) &&
	ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) &&
	ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) &&
	ns::Track::HitCut(1) &&
	ns::TrackerLayer::ChargeStatus(1) &&
	ns::TrackerLayer::ChargeInRange(1,0,6.0+0.8,CRT) &&
	ns::Track::ChiSquareLessThan(10.0,YSD, FIT, IL1) &&
	ns::Track::InnerFiducialVolume(FIT, IL1) &&
	ns::Track::L1FiducialVolume(FIT, IL1) && 
	innTr::ChResoInnLessThan(CRT,0.2);

//L1 eff
auto num_l1 =
	ns::TrackerLayer::ChargeInRange(1,0,6.0+0.8,CRT) &&
	ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) &&
    ns::Track::InnerFiducialVolume(FIT, IL1) &&
    ns::Track::L1FiducialVolume(FIT, IL1);
auto den_l1 = 
	ns::Trigger::HasPhysicsTrigger() &&
	ns::Tof::BetaInRange(0.0,inf,BTH) &&
	ns::Tof::ChargeInRange(6.0-0.75,6.0+0.75,UTC) &&
	ns::InnerTracker::HitPattern() &&
	ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
	ns::InnerTracker::ChargeInRange(6.0-0.3,6.0+0.7,CRT) &&
	ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN) &&
	ns::Track::L1FiducialVolume(FIT, INN) &&
	trd::trdUcharge(6.0-0.5,6.0+0.3) && 
	trd::trdLcharge(6.0-0.5,6.0+0.3) && 
	innTr::ChResoInnLessThan(CRT,0.2) &&
	TofSt::tofL1FiducialVolume();


//Inn trcker eff
auto num_track =
	ns::InnerTracker::HitPattern() &&
    ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
    ns::InnerTracker::ChargeInRange(6.0f-0.3f,6.0f+0.7,CRT) &&
    ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN);

auto den_track =
	ns::Trigger::HasPhysicsTrigger() &&
	TofSt::oneHitOnEachLayer() &&
	TofSt::tofBetaInRange(0.4,inf,BTH) &&
	TofSt::tofL1FiducialVolume() &&
	TrdSt::trdL1FiducialVolume() &&
	ns::Track::HitCut(1) &&
	ns::Track::HitCut(9) &&
	ns::TrackerLayer::ChargeInRange(1,6.0-0.4,6.0+0.7,CRT) &&
	ns::TrackerLayer::ChargeInRange(9,6.0-0.45,6.0+0.8,CRT) &&
	TofSt::tofChargeInRange(6.0-0.45,6.0+0.45, UTC) &&
	TofSt::tofChargeInRange(6.0-0.45,6.0+0.45, LTC) &&
	TofSt::tofChargeInRange(6.0-0.45,6.0+0.45, TOT) &&
	TrdSt::trdNHitsAtleast(12) &&
	TrdSt::trdUcharge(6.0-0.5,6.0+0.5) && 
	TrdSt::trdLcharge(6.0-0.5,6.0+0.5); 


//trigger
auto den_trig =
	ns::Tof::BetaInRange(0.0,inf,BTH) &&
	ns::Tof::ChargeInRange(6.0-0.75,6.0+0.75,UTC) &&
    ns::InnerTracker::HitPattern() &&
	ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
    ns::InnerTracker::ChargeInRange(6.0-0.3,6.0+0.7,CRT) &&
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, IL1) &&
    ns::Track::ChiSquareLessThan(10.0, YSD, FIT, INN) &&
    ns::Track::HitCut(1) &&
    ns::TrackerLayer::ChargeStatus(1) &&
	ns::TrackerLayer::ChargeInRange(1,0,6.0+0.8,CRT) &&
	ns::Track::InnerFiducialVolume(FIT, IL1) &&
	ns::Track::L1FiducialVolume(FIT, IL1)&&
	innTr::ChResoInnLessThan(CRT,0.2);
	
