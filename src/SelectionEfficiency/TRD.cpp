#include "./../../include/Efficiencies/TRD.h"

namespace Efficiency {
namespace trd { //TRD

trd::trdUcharge::trdUcharge(float min, float max) {
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trdKBase->Charge[NAIA::TrdK::ChargeType::Upper] > min && event.trdKBase->Charge[NAIA::TrdK::ChargeType::Upper] < max);
    });
}


trd::trdLcharge::trdLcharge(float min, float max) {
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trdKBase->Charge[NAIA::TrdK::ChargeType::Lower] > min && event.trdKBase->Charge[NAIA::TrdK::ChargeType::Lower] < max);
    });
}


trd::HasGoodTRDHits::HasGoodTRDHits(unsigned int charge) {
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
  case 15:
    thr1 = 520;
    thr2 = 0.87;
    break;
  case 16:
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
    auto check1 = (float)trdOffAmps / trdOffHits;
    auto check2 = (float)trdOffHits / (trdOnHits + trdOffHits);
    if ((trdOffHits && check1 < thr1) && ((trdOnHits + trdOffHits) > 0 && check2 < thr2))
      goodTRDhits = true;
    return goodTRDhits;
  });
};

} //TRD

namespace TrdSt{ // TRD Standalone

TrdSt::trdL1FiducialVolume::trdL1FiducialVolume(){
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      if(!event.CheckMask(NAIA::Category::HasTrdStandalone)) return false;
      constexpr double fid_R[] = {62, 62, 46, 46, 46, 46, 46, 46};
      constexpr double fid_Y[] = {47, 40, 44, 44, 36, 36, 44, 44};
      auto pos = event.trdKBaseSt->InterpolateAtZ(NAIA::TrTrack::fitPositionHeightZ[NAIA::TrTrack::FitPositionHeight::Layer1]); //x,y
      if(sqrt(pos[0]*pos[0]+pos[1]*pos[1]) > fid_R[NAIA::TrTrack::FitPositionHeight::Layer1] || fabs(pos[1]) > fid_Y[NAIA::TrTrack::FitPositionHeight::Layer1]) return false;
      return true;
    });
}

TrdSt::trdNHitsAtleast::trdNHitsAtleast(float min) {
		m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
			return (event.trdKBaseSt->NHits[NAIA::TrdK::QualType::OffTrack] + event.trdKBaseSt->NHits[NAIA::TrdK::QualType::OnTrack] >= min);
		});	
}

TrdSt::trdUcharge::trdUcharge(float min, float max) {
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trdKBaseSt->Charge[NAIA::TrdK::ChargeType::Upper] > min && event.trdKBaseSt->Charge[NAIA::TrdK::ChargeType::Upper] < max);
    });
}

TrdSt::trdLcharge::trdLcharge(float min, float max) {
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trdKBaseSt->Charge[NAIA::TrdK::ChargeType::Upper] > min && event.trdKBaseSt->Charge[NAIA::TrdK::ChargeType::Lower] < max);
    });
}


} // TRD Standalone
} //Efficiency
