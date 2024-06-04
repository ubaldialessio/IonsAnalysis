#include "BuildTemplatesSel.h"

namespace BuildTemplatesSel {

BuildTemplatesSel::HitPattern::HitPattern() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto pattern = event.trTrackBase->GetTrackPattern(NAIA::TrTrack::Side::X);

    return (pattern[1] && (pattern[2] || pattern[3]) && (pattern[4] || pattern[5]) && (pattern[6] || pattern[7]));
  });
}

BuildTemplatesSel::InnerTrackerChargeInRange::InnerTrackerChargeInRange(unsigned int charge,
                                                                        NAIA::TrTrack::ChargeRecoType recoType) {

  float min, max;
  switch (charge) {
  default:
    min = static_cast<float>(charge) - 0.3f;
    max = static_cast<float>(charge) + 0.3f;
  }

  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto charge = event.trTrackBase->InnerCharge[recoType];
    return (charge > min && charge < max);
  });
}

BuildTemplatesSel::InnerTrackerChargeRMSLessThan::InnerTrackerChargeRMSLessThan(
    float max, NAIA::TrTrack::ChargeRecoType recoType) {

  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
  	if (event.trTrackBase->InnerCharge[recoType]!=0)
    	return (event.trTrackBase->InnerChargeRMS[recoType] / event.trTrackBase->InnerCharge[recoType] < max);
  	return false; 
  });
}

BuildTemplatesSel::InnerTrackerNtrackLessThan::InnerTrackerNtrackLessThan(float max) {

  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return (event.evSummary->NTrTrack < max); });
}

BuildTemplatesSel::TrackerLayerChargeInRange::TrackerLayerChargeInRange(unsigned int layer, unsigned int charge,
                                                                        NAIA::TrTrack::ChargeRecoType recoType) {

  float min, max;

  if (layer == 1) {
    switch (charge) {
    default:
      min = static_cast<float>(charge) - 0.3f;
      max = static_cast<float>(charge) + 0.5f;
    }
  }

  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto jlayer = layer - 1;
    auto charge = event.trTrackBase->LayerChargeXY[jlayer][recoType];
    return (charge > min && charge < max);
  });
}

BuildTemplatesSel::TofChargeInRange::TofChargeInRange(unsigned int charge, NAIA::Tof::ChargeType type) {

  float min, max;

  switch (charge) {
  default:
    min = static_cast<float>(charge) - 0.5f;
    max = static_cast<float>(charge) + 0.5f;
  }

  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto charge = event.tofBase->Charge[type];
    return (charge > min && charge < max);
  });
}

BuildTemplatesSel::HasGoodSecondTrTrack::HasGoodSecondTrTrack(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType,
                                                              float min) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    bool good_sec_track = true;
    bool isAboveThreshold = false;
    bool hasMinumumHitsY = false;

    if (event.secondTrTrackBase->FitIDExists(fitType, spanType)) {

      if (event.secondTrTrackBase->Rigidity[fitType][spanType] > min)
        isAboveThreshold = true;

      unsigned int nHits = 0;
      for (unsigned int i = 1; i < 8; ++i) {
        if (NAIA::ContainsKeys(event.secondTrTrackBase->TrTrackHitPos, i, NAIA::TrTrack::Side::Y))
          ++nHits;
      }
      if (nHits > 3)
        hasMinumumHitsY = true;
    }

    if (isAboveThreshold && hasMinumumHitsY)
      good_sec_track = false;

    return good_sec_track;
  });
}

BuildTemplatesSel::HasGoodTRDHits::HasGoodTRDHits(unsigned int charge) {

  float thr1 = 1e6;
  float thr2 = 1;

  switch (charge) {
  case 3:
    thr1 = 200;
    thr2 = 0.6;
  case 4:
    thr1 = 220;
    thr2 = 0.68;
    break;
  case 5:
    thr1 = 250;
    thr2 = 0.7;
    break;
  case 6:
    thr1 = 320;
    break;
  case 7:
    thr1 = 250;
    thr2 = 0.75;
    break;
  case 8:
    thr1 = 400;
    break;
  case 9:
    thr1 = 320;
    thr2 = 0.83;
    break;
  case 10:
    thr1 = 350;
    thr2 = 0.84;
    break;
  case 11:
    thr1 = 420;
    thr2 = 0.85;
    break;
  case 12:
    thr1 = 420;
    thr2 = 0.86;
    break;
  case 13:
    thr1 = 420;
    thr2 = 0.86;
    break;
  case 14:
    thr1 = 520;
    thr2 = 0.87;
    break;
  default:
    thr1 = 400;
  }

  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    bool goodTRDhits = false;
    auto trdOffHits = event.trdKBase->NHits[NAIA::TrdK::QualType::OffTrack];
    auto trdOnHits = event.trdKBase->NHits[NAIA::TrdK::QualType::OnTrack];
    auto trdOffAmps = event.trdKBase->Amps[NAIA::TrdK::QualType::OffTrack];
    auto trdOnAmps = event.trdKBase->Amps[NAIA::TrdK::QualType::OnTrack];
    float check1, check2;

	if (trdOffHits!=0)
    	check1 = (float)trdOffAmps / trdOffHits;
    if ( (trdOnHits+trdOffHits)!=0 )
    	check2 = (float)trdOffHits / (trdOnHits + trdOffHits);

    if ((trdOffHits && check1 < thr1) && ((trdOnHits + trdOffHits) > 0 && check2 < thr2))
      goodTRDhits = true;

    return goodTRDhits;
  });
}

} // namespace BuildTemplatesSel
