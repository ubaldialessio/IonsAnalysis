
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
namespace TrTrackEffSel {

class TofChargeInRange : public NSL::Selection {
public:
  TofChargeInRange(unsigned int charge, NAIA::Tof::ChargeType type, std::string option);
};

class TofChi2TimeLessThan : public NSL::Selection {
public:
  TofChi2TimeLessThan(float max);
};

class TofChi2CooLessThan : public NSL::Selection {
public:
  TofChi2CooLessThan(float max);
};

class UnbExtHitChargeStatus : public NSL::Selection {
public:
  UnbExtHitChargeStatus(NAIA::UnbExtHitBaseData::ExtHit external_layer);
};

class UnbExtHitChargeInRange : public NSL::Selection {
public:
  UnbExtHitChargeInRange(NAIA::UnbExtHitBaseData::ExtHit external_layer, unsigned int charge,
                         NAIA::TrTrack::ChargeRecoType recoType, std::string option);
};

inline bool IsInsideLayer(unsigned int jlayer, float margin, NAIA::Event &event);

class IsInsideL1 : public NSL::Selection {
public:
  IsInsideL1(float margin);
};

class IsInsideL9 : public NSL::Selection {
public:
  IsInsideL9(float margin);
};

class IsInsideInner : public NSL::Selection {
public:
  IsInsideInner(float margin);
};

inline bool IsInsideFiducial(unsigned int jlayer, NAIA::Event &event);

class InnerFiducialVolume : public NSL::Selection {
public:
  InnerFiducialVolume();
};

class L1FiducialVolume : public NSL::Selection {
public:
  L1FiducialVolume();
};

class L9FiducialVolume : public NSL::Selection {
public:
  L9FiducialVolume();
};

class IsInsideEcal : public NSL::Selection {
public:
  IsInsideEcal();
};

} // namespace TrTrackEffSel
} // namespace Efficiency
