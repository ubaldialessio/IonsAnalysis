#ifndef TOFSTANDALONE_H
#define TOFSTANDALONE_H
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
namespace TofStandalone{ 

class oneHitOnEachLayer : public NSL::Selection {
public:
	oneHitOnEachLayer();
};

class goodPathAll : public NSL::Selection{
public:
	goodPathAll();
};

class tofBetaInRange : public NSL::Selection{
public:
  tofBetaInRange(double min, double max, NAIA::Tof::BetaType type);
};

class tofChi2TimeLessThan : public NSL::Selection{
public:
  tofChi2TimeLessThan(double max);
}
;
class tofChi2CooLessThan : public NSL::Selection{
public:
  tofChi2CooLessThan(double max);
};

class tofChargeInRange : public NSL::Selection{
public:
  tofChargeInRange(double min, double max, NAIA::Tof::ChargeType type);
};

bool tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight layer, NAIA::Event &event);

class tofInnerFiducialVolume : public NSL::Selection{
public:
  tofInnerFiducialVolume();
};

class tofL1FiducialVolume : public NSL::Selection{
public:
  tofL1FiducialVolume();
};

short bitcounter (short n);

} // TofStandalone
} //Efficiency
#endif
