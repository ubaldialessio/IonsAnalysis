#include "./../../include/Templates/BuildTemplatesSel.h"

namespace BuildTemplatesSel {

BuildTemplatesSel::HitPattern::HitPattern() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto pattern = event.trTrackBase->GetTrackPattern();

    return (pattern[1] > NAIA::HitClusterAssociation::None &&
            (pattern[2] > NAIA::HitClusterAssociation::None || pattern[3] > NAIA::HitClusterAssociation::None) &&
            (pattern[4] > NAIA::HitClusterAssociation::None || pattern[5] > NAIA::HitClusterAssociation::None) &&
            (pattern[6] > NAIA::HitClusterAssociation::None || pattern[7] > NAIA::HitClusterAssociation::None));
  });
}
// calculate L3-L8 charge and RMS for L2 template selection
BuildTemplatesSel::L3toL8ChargeSelection::L3toL8ChargeSelection(unsigned int charge,
                                     NAIA::TrTrack::ChargeRecoType recoType) {
        float min, max;
        switch (charge) {
        default:
            min = static_cast<float>(charge) - 0.5f;
            max = static_cast<float>(charge) + 0.5f;
        }
        m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
            // Calculate the average charge (L3-L8), removing the maximum value
            float In_Charge = 0.0f;
            float In_RMS = 0.0f;
            int validLayerCount = 0;
            float maxCharge = 0.0f;

            for (unsigned int jlayer = 2; jlayer < 8; ++jlayer) { // L3 to L8 (index 2 to 7)
                float chargeVal = event.trTrackBase->LayerChargeXY[jlayer][recoType];
                if (chargeVal < 0) continue; 
                if (chargeVal > maxCharge) maxCharge = chargeVal;
                In_Charge += chargeVal;
                In_RMS += std::pow(chargeVal, 2);
                validLayerCount++;
            }
            if (validLayerCount > 1) {
                In_Charge = (In_Charge - maxCharge) / (validLayerCount - 1); // Remove the maximum charge
                In_RMS = (In_RMS - std::pow(maxCharge, 2)) / (validLayerCount - 1);
                In_RMS = std::sqrt(std::fabs(In_RMS - In_Charge * In_Charge));
            } else {
                In_Charge = 0.0f;
                In_RMS = 0.0f;
            }
            
            // Add selection: In_RMS < 1
            //return (In_Charge > min && In_Charge < max && In_RMS < 1.0f);
            return (In_Charge > min && In_Charge < max);

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
      min = static_cast<float>(charge) - 0.5f;
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
