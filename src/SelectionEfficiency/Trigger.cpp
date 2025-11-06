#include "./../../include/Efficiencies/Trigger.h"

namespace Efficiency {
namespace TriggerEff {

bool IsInsideLayer(unsigned int jlayer, NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin,
                   NAIA::Event &event) {
  bool isinlayer = false;

  assert(jlayer > 0 && jlayer < 10);
  auto layer = static_cast<NAIA::TrTrack::FitPositionHeight>(jlayer - 1);

  static float tracker_planes_edges[9][4] = {
      {-62.14, -47.40, 62.14, 47.40}, {-62.14, -40.10, 62.14, 40.10}, {-49.70, -43.75, 49.70, 43.75},
      {-49.72, -43.75, 49.72, 43.75}, {-49.71, -36.45, 49.70, 36.45}, {-49.72, -36.45, 49.72, 36.45},
      {-49.72, -43.75, 49.71, 43.75}, {-49.72, -43.75, 49.71, 43.75}, {-45.62, -29.48, 45.55, 29.53}};

  float x = event.trTrackBase->TrTrackFitPos[layer][fitType][spanType][NAIA::TrTrack::Side::X];
  float y = event.trTrackBase->TrTrackFitPos[layer][fitType][spanType][NAIA::TrTrack::Side::Y];

  if ((x > tracker_planes_edges[layer][0] + margin) && (x < tracker_planes_edges[layer][2] - margin) &&
      (y > tracker_planes_edges[layer][1] + margin) && (y < tracker_planes_edges[layer][3] - margin)) {

    if (layer == NAIA::TrTrack::FitPositionHeight::Layer9) {
      // remove dead area center L9
      if (x < -0.5 || y > 0.5) {
        isinlayer = true;
      }
    } else {
      // circle
      if ((sqrt(x * x + y * y) < tracker_planes_edges[layer][2] - margin))
        isinlayer = true;
    }
  }
  return isinlayer;
}

TriggerEff::HasTrigger::HasTrigger(unsigned int triggerPatt) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto hasTrigger = false;
    if (event.evSummary->PhysBPatt == triggerPatt)
      hasTrigger = true;
    return (hasTrigger);
   });
 }

TriggerEff::InnerTrackerChargeInRange::InnerTrackerChargeInRange(float min, float max,
                                                             NAIA::TrTrack::ChargeRecoType recoType) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto charge = event.trTrackBase->InnerCharge[recoType];
    return (charge > min && charge < max);
  });
}

TriggerEff::TofChargeInRange::TofChargeInRange(float min, float max, NAIA::Tof::ChargeType type) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto charge = event.tofBase->Charge[type];
    return (charge > min && charge < max);
  });
}


TriggerEff::IsInsideL1::IsInsideL1(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    if (!event.trTrackBase->FitIDExists(fitType, spanType))
      return false;

    return IsInsideLayer(1, fitType, spanType, margin, event);
  });
}

} // namespace TriggerEffSel
} // namespace Efficiency

