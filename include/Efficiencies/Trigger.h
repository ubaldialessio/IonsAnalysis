#ifndef TRIGGER_H
#define TRIGGER_H
// NAIA headers
#include <Chain/NAIAChain.h>

// NSL headers
#include <NSL/AllSelections.h>

// ROOT headers
#include <TGraph.h>
#include <TH2D.h>
#include <TTimeStamp.h>

// c++ headers
#include <array>
#include <memory>
#include <numeric>

namespace Efficiency {
namespace TriggerEff{

inline bool IsInsideLayer(unsigned int jlayer, NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin,
                   NAIA::Event &event);

class HasTrigger : public NSL::Selection {
public:
  HasTrigger(unsigned int triggerPatt);
};

class IsInsideL1 : public NSL::Selection {
public:
  IsInsideL1(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin);
};

class InnerTrackerChargeInRange : public NSL::Selection {
public:
  InnerTrackerChargeInRange(float min, float max,NAIA::TrTrack::ChargeRecoType recoType);
};

class TofChargeInRange : public NSL::Selection {
public:
  TofChargeInRange(float min, float max, NAIA::Tof::ChargeType type);
};


}
}
#endif
