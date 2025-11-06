#ifndef INNERTRACKER_H
#define INNERTRACKER_H
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
namespace InnerTracker{

class ChResoInnLessThan : public NSL::Selection { 
public:
  ChResoInnLessThan(NAIA::TrTrack::ChargeRecoType charge, float max);
};

class L1NormResidualLessThan : public NSL::Selection {
public:
  L1NormResidualLessThan(float max, NAIA::TrTrack::Fit fit);
};

class NGoodClustersGreaterThan : public NSL::Selection {
public:
  NGoodClustersGreaterThan(unsigned int min, int mask);
};

class ChargeWithExternal : public NSL::Selection {
public:
  ChargeWithExternal(float min,float max,NAIA::TrTrack::ChargeRecoType recoType);
};

}
}
#endif
