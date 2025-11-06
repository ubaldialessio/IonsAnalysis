#include "./../../include/Efficiencies/InnerTracker.h"
namespace Efficiency {
namespace InnerTracker { //innTr

InnerTracker::ChResoInnLessThan::ChResoInnLessThan(NAIA::TrTrack::ChargeRecoType charge, float max) {
    m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trTrackBase->InnerChargeRMS[charge]/(event.trTrackBase->InnerCharge[charge]) < max );
    });
}
  
InnerTracker::L1NormResidualLessThan::L1NormResidualLessThan(float max, NAIA::TrTrack::Fit fit) {
  	m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
  	    const auto &testContainer = event.trTrackBase->TrTrackHitPos; 	
  	    if (!event.trTrackBase->FitIDExists(fit, NAIA::TrTrack::Span::InnerL1))
  	      return false;
  	    int ndof_inner = -3;
  	    for (int il = 1; il < 8; il++) {
  	      if (NAIA::ContainsKeys(testContainer, il, NAIA::TrTrack::Y)) {
  	        ndof_inner++;
  	      }
  	    }
  	    float residual =
  	        event.trTrackBase->TrChiSq[fit][NAIA::TrTrack::Span::InnerL1][NAIA::TrTrack::Side::Y] * (ndof_inner + 1) -
  	        event.trTrackBase->TrChiSq[fit][NAIA::TrTrack::Span::InnerOnly][NAIA::TrTrack::Side::Y] * ndof_inner;
  	    return residual < max;
  	  });
  }


InnerTracker::NGoodClustersGreaterThan::NGoodClustersGreaterThan(unsigned int min, int mask) {
  m_matcher = std::make_unique<NSL::boolMatcher>([=](Event &event) {
    unsigned int nCl = 0;
    for (unsigned int i = 1; i < 8; ++i) {
      if (NAIA::ContainsKeys(event.trTrackBase->LayerChargeStatus, i) &&
          (event.trTrackBase->LayerChargeStatus[i] & mask) == 0)
        ++nCl;
    }
    return nCl > min;
  });
}

InnerTracker::ChargeWithExternal::ChargeWithExternal(float min, float max,NAIA::TrTrack::ChargeRecoType recoType) {
  m_matcher = std::make_unique<NSL::boolMatcher>([=](Event &event) {
  	auto charge = event.trTrackBase->Charge[recoType];
	return (charge > min && charge < max);
  });
}

} //InnerTracker
} //Efficiency
