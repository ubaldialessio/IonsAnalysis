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

class HasNoGoodSecondTrTrack : public NSL::Selection {
public:
  HasNoGoodSecondTrTrack(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType,float min);
};

class SecondTrkChiSquareLessThan : public NSL::Selection {
public:
  SecondTrkChiSquareLessThan(float max, NAIA::TrTrack::Side side,
                                  NAIA::TrTrack::Fit fit, NAIA::TrTrack::Span span);
};

class SecondTrackYHitsGreaterThan : public NSL::Selection {
public:
  SecondTrackYHitsGreaterThan(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType,float min);
};

class IsSecondTrackOnDiagonal : public NSL::Selection {
public:
  IsSecondTrackOnDiagonal(NAIA::TrTrack::Fit fit, NAIA::TrTrack::Span span, double delta);
};

class L1ClusterCutBetween : public NSL::Selection {
public:
  L1ClusterCutBetween(int min, int max, NAIA::TrTrack::Side side, NAIA::TrTrack::DistanceFromTrack d);
};

}//ms
}//Efficiency
#endif
