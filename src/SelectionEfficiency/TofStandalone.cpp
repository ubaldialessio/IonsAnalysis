#include "./../../include/Efficiencies/TofStandalone.h"

namespace Efficiency {
namespace TofStandalone{ // TofSt


TofStandalone::oneHitOnEachLayer::oneHitOnEachLayer() {
		m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
		if (!event.header->CheckMask(NAIA::Category::HasTofStandalone)) return false;
		return (event.tofPlusSt->NClusters[0][NAIA::Tof::OnTime] > 0 && event.tofPlusSt->NClusters[1][NAIA::Tof::OnTime] > 0 &&
			 event.tofPlusSt->NClusters[2][NAIA::Tof::OnTime] > 0 && event.tofPlusSt->NClusters[3][NAIA::Tof::OnTime] > 0 );
		});
}	


TofStandalone::goodPathAll::goodPathAll() {
		m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
			return (event.tofPlusSt->LayerGoodPathl[0b1111] );
		});	
}

TofStandalone::tofBetaInRange::tofBetaInRange(double min, double max, NAIA::Tof::BetaType type){
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      auto beta = event.tofBaseSt->Beta[type];
      return (beta > min && beta < max);
    });
}


TofStandalone::tofChi2TimeLessThan::tofChi2TimeLessThan(double max){
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.tofPlusSt->Chi2Time < max);
    });
}

TofStandalone::tofChi2CooLessThan::tofChi2CooLessThan(double max){
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.tofPlusSt->Chi2Coo < max);
    });
}

TofStandalone::tofChargeInRange::tofChargeInRange(double min, double max, NAIA::Tof::ChargeType type){
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      auto charge = event.tofBaseSt->ChargeNoPL[type];
      return (charge > min && charge < max);
    });
}

bool tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight layer, NAIA::Event &event) {
  constexpr double fid_R[] = {62, 62, 46, 46, 46, 46, 46, 46};
  constexpr double fid_Y[] = {47, 40, 44, 44, 36, 36, 44, 44};
  if (!event.header->CheckMask(NAIA::Category::HasTofStandalone)) return false;
  auto pos = event.tofBaseSt->InterpolateAtZ(NAIA::TrTrack::fitPositionHeightZ[layer]); //x,y
  if(sqrt(pos[0]*pos[0]+pos[1]*pos[1]) > fid_R[layer] || fabs(pos[1]) > fid_Y[layer]) return false;
  return true;
}

TofStandalone::tofInnerFiducialVolume::tofInnerFiducialVolume() {
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      short hits = 0;
      if(tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight::Layer2, event)) hits |= 1<<0;
      if(tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight::Layer3, event)) hits |= 1<<1;
      if(tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight::Layer4, event)) hits |= 1<<2;
      if(tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight::Layer5, event)) hits |= 1<<3;
      if(tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight::Layer6, event)) hits |= 1<<4;
      if(tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight::Layer7, event)) hits |= 1<<5;
      if(tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight::Layer8, event)) hits |= 1<<6;
      if(bitcounter(hits)     <5) return false; //check number of hits
      if(bitcounter(hits&1)  ==0) return false; //check L2 hit
      if(bitcounter(hits&6)  ==0) return false; //check L3||L4 hit
      if(bitcounter(hits&24) ==0) return false; //check L5||L6 hit
      if(bitcounter(hits&96) ==0) return false; //check L7||L8 hit
      return true;
    });
}

TofStandalone::tofL1FiducialVolume::tofL1FiducialVolume(){
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      return tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight::Layer1, event);
    });
}

short bitcounter(short n) {
  if(n==0) return 0;
  else return (n&1)+bitcounter(n>>1);
}


} // TofSt
} // Efficiency
