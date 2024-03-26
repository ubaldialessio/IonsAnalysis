#ifndef MYSEL_H
#define MYSEL_H
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
namespace MySel { //ms

inline bool IsInsideLayer(unsigned int jlayer, NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin,
                   NAIA::Event &event);

class IsInsideL1 : public NSL::Selection {
public:
  IsInsideL1(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin);
};

class HasGoodSecondTrTrack : public NSL::Selection {
public:
  HasGoodSecondTrTrack(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType,float min);
};

}//ms
}//Efficiency
#endif
