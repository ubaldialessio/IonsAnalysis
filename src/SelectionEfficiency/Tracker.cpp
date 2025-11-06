#include "./../../include/Efficiencies/Tracker.h"

namespace Efficiency{
namespace TrTrackEffSel{
TrTrackEffSel::IsInsideEcal::IsInsideEcal() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    float Ecal_fXmin = -32.270;
    float Ecal_fXmax = 32.530;
    float Ecal_fYmin = -32.470;
    float Ecal_fYmax = 32.330;
    auto Ecal_fZtop = NAIA::TrTrack::FitPositionHeight::EcalTop;
    auto Ecal_fZbot = NAIA::TrTrack::FitPositionHeight::EcalBottom;
    float xT, yT, xB, yB;
    if (event.header->CheckMask(NAIA::Category::HasTofStandalone)) {
      xT = event.tofBaseSt->InterpolateAtZ(Ecal_fZtop)[0];
      yT = event.tofBaseSt->InterpolateAtZ(Ecal_fZtop)[1];
      xB = event.tofBaseSt->InterpolateAtZ(Ecal_fZbot)[0];
      yB = event.tofBaseSt->InterpolateAtZ(Ecal_fZbot)[1];
    } else
      return false;
    if (xT < Ecal_fXmin || yT < Ecal_fYmin || xT > Ecal_fXmax || yT > Ecal_fYmax || xB < Ecal_fXmin ||
        yB < Ecal_fYmin || xB > Ecal_fXmax || yB > Ecal_fYmax || false)
      return false;
    return true;
  });
}

TrTrackEffSel::UnbExtHitChargeStatus::UnbExtHitChargeStatus(NAIA::UnbExtHitBaseData::ExtHit external_layer) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto result = false;
    auto chargestatus = event.extHitBase->ChargeStatus(external_layer);
    if ((chargestatus & 0x10013D) == 0)
      result = true;
    else
      result = false;
    return result;
  });
}

TrTrackEffSel::IsHitPresent::IsHitPresent(NAIA::UnbExtHitBaseData::ExtHit external_layer) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    return event.extHitBase->IsHitPresent(external_layer);
  });
}

TrTrackEffSel::UnbExtHitChargeInRange::UnbExtHitChargeInRange(NAIA::UnbExtHitBaseData::ExtHit external_layer,
                                                              unsigned int charge,
                                                              NAIA::TrTrack::ChargeRecoType recoType,
                                                              std::string option) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto hitcharge = event.extHitBase->Charge(external_layer, recoType);
    float min, max = 0;
    if (option=="UNB") {
      if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L1) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        default:
          min = (charge)/2.0f;
          max = (charge)+999.f;
        }
      } else if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L9) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        default:
          min = (charge)-0.5f;
          max = (charge) + 1.0f;
        }
      }
    } else if (option.empty() || option.substr(0, 2) == "PG") {
      if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L1) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        default:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        }
      } else if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L9) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        default:
          min = (charge)-0.5f;
          max = (charge) + 1.0f;
        }
      }
    } else if (option == "qi") {
      if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L1) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        case 9:
        case 11:
        case 13:
        case 10:
        case 12:
        case 14:
        case 16:
          min = (charge)-0.0585f * pow(charge, 1.15f) - 0.35f;
          max = (charge) + 0.0334f * pow(charge, 1.15f) + 0.20f;
        default:
              min = charge - 0.0585f * pow(charge, 1.15f) - 0.35f;
              max = charge + 0.0334f * pow(charge, 1.15f) + 0.20f;
        }
      }
      if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L9) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        case 9:
        case 11:
        case 13:
        case 10:
        case 12:
        case 14:
        case 16:
          min = (charge)-0.0284f * pow(charge, 1.15f) - 0.17f;
          max = (charge) + 0.0585f * pow(charge, 1.15f) + 0.35f;
        default:
          min = (charge)-0.5f;
          max = (charge) + 1.5f + (charge - 3.0f) * 0.06f;
        }
      }
    }
    return (hitcharge > min && hitcharge < max);
  });
}



TrTrackEffSel::ExtChargeInRange::ExtChargeInRange(double min, double max, NAIA::UnbExtHitBaseData::ExtHit layer, NAIA::TrTrack::ChargeRecoType recoType){
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      auto charge = event.extHitBase->Charge(layer, recoType);
      return (charge > min && charge < max);
    });
 }


bool IsInsideFiducial(unsigned int jlayer, NAIA::Event &event) {
  // Beware, for L9 the check is X only, not R
  constexpr std::array<float, 9> fid_R = {62, 62, 46, 46, 46, 46, 46, 46, 43};
  constexpr std::array<float, 9> fid_Y = {47, 40, 44, 44, 36, 36, 44, 44, 29};

  assert(jlayer > 0 && jlayer < 10);
  auto layer = static_cast<NAIA::TrTrack::FitPositionHeight>(jlayer - 1);

  constexpr std::array<float, 9> LayerPositionHeightZ = {
      158.920, 53.060, 29.228, 25.212, 1.698, -2.318, -25.212, -29.228, -135.882,
  };

  float x, y;
  if (event.header->CheckMask(NAIA::Category::HasTofStandalone)) {
    x = event.tofBaseSt->InterpolateAtZ(LayerPositionHeightZ.at(layer))[0];
    y = event.tofBaseSt->InterpolateAtZ(LayerPositionHeightZ.at(layer))[1];
  } else
    return false;

  if (layer == NAIA::TrTrack::FitPositionHeight::Layer9) {
    if (std::fabs(x) > fid_R[layer] || std::fabs(y) > fid_Y[layer])
      return false;
  } else {
    if (std::sqrt(x * x + y * y) > fid_R[layer] || std::fabs(y) > fid_Y[layer])
      return false;
  }

  return true;
}


bool IsInsideLayer(unsigned int jlayer, float margin, NAIA::Event &event) {
  bool isinlayer = false;

  assert(jlayer > 0 && jlayer < 10);
  auto layer = static_cast<NAIA::TrTrack::FitPositionHeight>(jlayer - 1);

  static float tracker_planes_edges[9][4] = {
      {-62.14, -47.40, 62.14, 47.40}, {-62.14, -40.10, 62.14, 40.10}, {-49.70, -43.75, 49.70, 43.75},
      {-49.72, -43.75, 49.72, 43.75}, {-49.71, -36.45, 49.70, 36.45}, {-49.72, -36.45, 49.72, 36.45},
      {-49.72, -43.75, 49.71, 43.75}, {-49.72, -43.75, 49.71, 43.75}, {-45.62, -29.48, 45.55, 29.53}};

  constexpr std::array<float, 9> LayerPositionHeightZ = {
      158.920, 53.060, 29.228, 25.212, 1.698, -2.318, -25.212, -29.228, -135.882,
  };

  float x, y;
  if (event.header->CheckMask(NAIA::Category::HasTofStandalone)) {
    x = event.tofBaseSt->InterpolateAtZ(LayerPositionHeightZ.at(layer))[0];
    y = event.tofBaseSt->InterpolateAtZ(LayerPositionHeightZ.at(layer))[1];
  } else
    isinlayer = false;

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

TrTrackEffSel::IsInsideL1::IsInsideL1(float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return IsInsideLayer(1, margin, event); });
}

TrTrackEffSel::IsInsideL9::IsInsideL9(float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return IsInsideLayer(9, margin, event); });
}

TrTrackEffSel::IsInsideInner::IsInsideInner(float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    if (IsInsideLayer(2, margin, event) && (IsInsideLayer(3, margin, event) || IsInsideLayer(4, margin, event)) &&
        (IsInsideLayer(5, margin, event) || IsInsideLayer(6, margin, event)) &&
        (IsInsideLayer(7, margin, event) || IsInsideLayer(8, margin, event)) && true)
      return true;
    else
      return false;
  });
}

TrTrackEffSel::InnerFiducialVolume::InnerFiducialVolume() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    int nhits = 0;
    for (auto layer = 2; layer < 9; layer++)
      if (IsInsideFiducial(layer, event))
        nhits++;
    if (nhits < 5)
      return false;

    if (!IsInsideFiducial(2, event))
      return false;
    if (!(IsInsideFiducial(3, event) || IsInsideFiducial(4, event)))
      return false;
    if (!(IsInsideFiducial(5, event) || IsInsideFiducial(6, event)))
      return false;
    if (!(IsInsideFiducial(7, event) || IsInsideFiducial(8, event)))
      return false;

    return true;
  });
}

TrTrackEffSel::L1FiducialVolume::L1FiducialVolume() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return IsInsideFiducial(1, event); });
}

TrTrackEffSel::L9FiducialVolume::L9FiducialVolume() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return IsInsideFiducial(9, event); });
}

TrTrackEffSel::TofStCharge::TofStCharge(float min, float max, NAIA::Tof::ChargeType type) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto charge = event.tofBaseSt->ChargeNoPL[type];
    return (charge > min && charge < max);
  });
}
 

}//TrTrackEffSel
}//Efficiency
