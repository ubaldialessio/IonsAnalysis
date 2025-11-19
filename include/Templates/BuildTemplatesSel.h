#ifndef BUILDTEMPLATESSEL_H
#define BUILDTEMPLATESSEL_H
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

namespace BuildTemplatesSel {

class ChargeSigma {
public:
    // imposta una sigma per un dato Z
    static void SetSigma(unsigned int Z, float sigma);

    // recupera la sigma di un dato Z
    static float GetSigma(unsigned int Z);
};

class HitPattern : public NSL::Selection {
public:
  HitPattern();
};

class L3toL8ChargeSelection : public NSL::Selection {
public:
  L3toL8ChargeSelection(unsigned int charge, NAIA::TrTrack::ChargeRecoType recoType, float n_sigma = 1.2f);
};

class InnerTrackerChargeInRange : public NSL::Selection {
public:
  InnerTrackerChargeInRange(unsigned int charge, NAIA::TrTrack::ChargeRecoType recoType, float n_sigma = 1.2f);
};

class InnerTrackerChargeRMSLessThan : public NSL::Selection {
public:
  InnerTrackerChargeRMSLessThan(float max, NAIA::TrTrack::ChargeRecoType recoType);
};

class InnerTrackerNtrackLessThan : public NSL::Selection {
public:
  InnerTrackerNtrackLessThan(float max);
};

class TrackerLayerChargeInRange : public NSL::Selection {
public:
  TrackerLayerChargeInRange(unsigned int layer, unsigned int charge, NAIA::TrTrack::ChargeRecoType recoType, float n_sigma = 1.2f);
};

class TofChargeInRange : public NSL::Selection {
public:
  TofChargeInRange(unsigned int charge, NAIA::Tof::ChargeType type, float n_sigma = 1.2f);
};

class HasGoodSecondTrTrack : public NSL::Selection {
public:
  HasGoodSecondTrTrack(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float min);
};  

class HasGoodTRDHits : public NSL::Selection {
public:
  HasGoodTRDHits(unsigned int charge);
}; 

} // namespace BuildTemplatesSel
#endif
