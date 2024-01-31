#include "binning.h"

//RTI Default selection
namespace NSL {
namespace Selections {

bool DefaultRTISelection(const NAIA::RTIInfo &rti_info) {
  if (rti_info.LivetimeFraction < 0.5 || rti_info.Zenith > 40 || rti_info.nTrigger / rti_info.nEvent < 0.98 ||
      rti_info.nError < 0 || rti_info.nError / rti_info.nEvent > 0.1 || rti_info.MeanAlignDiffExtLayer[0][1] > 35 ||
      rti_info.MeanAlignDiffExtLayer[1][1] > 45)
    return false;
  return true;
}
} // namespace Selections
} // namespace NSL


//trigger selection
auto mysel =
	ns::Trigger::HasPhysicsTrigger() &&
	ns::Tof::BetaInRange(0.0f,1.0f,BTH) &&
	ns::Tof::ChargeInRange(6.0f-0.75f,6.0f+0.75f,UTC) &&
    ns::InnerTracker::HitPattern() &&
    ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
    ns::InnerTracker::ChargeInRange(6.0f-0.3f,6.0f+0.7,CRT) &&
    ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN) &&
	ns::TrackerLayer::ChargeInRange(1,0,6.0f+0.8f,CRT) &&
	ns::Track::ChiSquareLessThan(10.0f,YSD, FIT, IL1) &&
	ns::Track::InnerFiducialVolume(FIT, IL1) &&
	ns::Track::L1FiducialVolume(FIT, IL1);

//efficiencies
auto num_tof =
	ns::Tof::BetaInRange(0.0f,1.0f,BTH) &&
	ns::Tof::ChargeInRange(6.0f-0.75f,6.0f+0.75f,UTC);
auto den_tof = 
	ns::Trigger::HasPhysicsTrigger() &&
	ns::InnerTracker::HitPattern() &&
	ns::InnerTracker::NHitsGreaterThan(4, YSD) &&
	ns::InnerTracker::ChargeInRange(6.0f-0.25f,6.0f+0.6,CRT) && //stricter
	ns::Track::ChiSquareLessThan(10.0f, YSD, FIT, INN) &&
	ns::TrackerLayer::ChargeInRange(1,6.0f-0.4f,6.0f+0.4f,CRT) && //stricter
	ns::Track::ChiSquareLessThan(10.0f,YSD, FIT, IL1) &&
	ns::Track::InnerFiducialVolume(FIT, IL1) &&
	ns::Track::L1FiducialVolume(FIT, IL1);
