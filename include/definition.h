#include <Chain/NAIAChain.h>
#include <NSL/AllSelections.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <Math/Minimizer.h>
#include <TVirtualFitter.h>
#include <TEfficiency.h>
#include "TString.h"
#include "TTimeStamp.h"
#include "TKey.h"
#include <TROOT.h>
#include <TSystem.h>
#include <string>
#include <vector>

#include "Efficiencies/TrTrackEffSel.h"

using namespace std;
namespace ns = NSL::Selections;
int fov = 1; // Fields of view: 0 = 25 deg, 1 = 30 deg, 2 = 35 deg, 3 = 40 deg
float inf = 9999999;


auto CRT = NAIA::TrTrack::ChargeRecoType::YJ;
auto FIT = NAIA::TrTrack::Fit::GBL;
auto INN = NAIA::TrTrack::Span::InnerOnly;
auto IL1 = NAIA::TrTrack::Span::InnerL1;
auto YSD = NAIA::TrTrack::Side::Y;
auto XSD = NAIA::TrTrack::Side::X;
auto UTC = NAIA::Tof::ChargeType::Upper;
auto LTC = NAIA::Tof::ChargeType::Lower;
auto TOT = NAIA::Tof::ChargeType::Total;
auto BTH = NAIA::Tof::BetaType::BetaH;
auto ON  = NAIA::Tof::OnTime;
auto OFF = NAIA::Tof::OnTime;
auto L1  = NAIA::TrTrack::FitPositionHeight::Layer1;
auto L2  = NAIA::TrTrack::FitPositionHeight::Layer2;
auto L3  = NAIA::TrTrack::FitPositionHeight::Layer3;
auto L4  = NAIA::TrTrack::FitPositionHeight::Layer4;
auto L5  = NAIA::TrTrack::FitPositionHeight::Layer5;
auto L6  = NAIA::TrTrack::FitPositionHeight::Layer6;
auto L7  = NAIA::TrTrack::FitPositionHeight::Layer7;
auto L8  = NAIA::TrTrack::FitPositionHeight::Layer8;
auto EL1 = NAIA::UnbExtHitBaseData::ExtHit::L1;
auto TRD_UP = NAIA::TrdK::ChargeType::Upper;
auto TRD_LOW = NAIA::TrdK::ChargeType::Lower;
auto TRD_ONTRACK = NAIA::TrdK::QualType::OnTrack;
auto TRD_OFFTRACK = NAIA::TrdK::QualType::OffTrack;


constexpr double fid_R[] = {62, 62, 46, 46, 46, 46, 46, 46};
constexpr double fid_Y[] = {47, 40, 44, 44, 36, 36, 44, 44};

short bitcounter(short n) {
  if(n==0) return 0;
  else return (n&1)+bitcounter(n>>1);
}

namespace Efficiency{
namespace TrTrackEffSel{
TrTrackEffSel::IsInsideEcal::IsInsideEcal() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    float Ecal_fXmin = -32.270;
    float Ecal_fXmax = 32.530;
    float Ecal_fYmin = -32.470;
    float Ecal_fYmax = 32.330;
    auto Ecal_fZtop = NAIA::TrTrack::FitPositionHeight::EcalTop;
    auto Ecal_fZbot = NAIA::TrTrack::FitPositionHeight::EcalBottom;
    float xT, yT, xB, yB;
    if (event.header->CheckMask(NAIA::Category::HasTofStandalone)) {
      xT = event.tofBaseSt->InterpolateAtZ(Ecal_fZtop)[0];
      yT = event.tofBaseSt->InterpolateAtZ(Ecal_fZtop)[1];
      xB = event.tofBaseSt->InterpolateAtZ(Ecal_fZbot)[0];
      yB = event.tofBaseSt->InterpolateAtZ(Ecal_fZbot)[1];
    } else
      return false;
    if (xT < Ecal_fXmin || yT < Ecal_fYmin || xT > Ecal_fXmax || yT > Ecal_fYmax || xB < Ecal_fXmin ||
        yB < Ecal_fYmin || xB > Ecal_fXmax || yB > Ecal_fYmax || false)
      return false;
    return true;
  });
}

TrTrackEffSel::UnbExtHitChargeStatus::UnbExtHitChargeStatus(NAIA::UnbExtHitBaseData::ExtHit external_layer) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto result = false;
    auto chargestatus = event.extHitBase->ChargeStatus(external_layer);
    if ((chargestatus & 0x10013D) == 0)
      result = true;
    else
      result = false;
    return result;
  });
}

TrTrackEffSel::UnbExtHitChargeInRange::UnbExtHitChargeInRange(NAIA::UnbExtHitBaseData::ExtHit external_layer,
                                                              unsigned int charge,
                                                              NAIA::TrTrack::ChargeRecoType recoType,
                                                              std::string option) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto hitcharge = event.extHitBase->Charge(external_layer, recoType);
    float min, max = 0;
    if (option.empty() || option.substr(0, 2) == "PG") {
      if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L1) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        default:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        }
      } else if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L9) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        default:
          min = (charge)-0.5f;
          max = (charge) + 1.0f;
        }
      }
    } else if (option == "qi") {
      if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L1) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        case 9:
        case 11:
        case 13:
        case 10:
        case 12:
        case 14:
        case 16:
          min = (charge)-0.0585f * pow(charge, 1.15f) - 0.35f;
          max = (charge) + 0.0334f * pow(charge, 1.15f) + 0.20f;
        default:
          min = (charge)-0.46f - (charge - 3.0f) * 0.1f;
          max = (charge) + ((charge < 5.0f) ? 0.65f : 0.65f + (charge - 5.0f) * 0.03f);
        }
      }
      if (external_layer == NAIA::UnbExtHitBaseData::ExtHit::L9) {
        switch (charge) {
        case 1:
        case 2:
          min = (charge)-0.4f;
          max = (charge) + 0.9f;
        case 9:
        case 11:
        case 13:
        case 10:
        case 12:
        case 14:
        case 16:
          min = (charge)-0.0284f * pow(charge, 1.15f) - 0.17f;
          max = (charge) + 0.0585f * pow(charge, 1.15f) + 0.35f;
        default:
          min = (charge)-0.5f;
          max = (charge) + 1.5f + (charge - 3.0f) * 0.06f;
        }
      }
    }
    return (hitcharge > min && hitcharge < max);
  });
}



TrTrackEffSel::ExtChargeInRange::ExtChargeInRange(double min, double max, NAIA::UnbExtHitBaseData::ExtHit layer, NAIA::TrTrack::ChargeRecoType recoType){
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      auto charge = event.extHitBase->Charge(layer, recoType);
      return (charge > min && charge < max);
    });
 }


bool IsInsideFiducial(unsigned int jlayer, NAIA::Event &event) {
  // Beware, for L9 the check is X only, not R
  constexpr std::array<float, 9> fid_R = {62, 62, 46, 46, 46, 46, 46, 46, 43};
  constexpr std::array<float, 9> fid_Y = {47, 40, 44, 44, 36, 36, 44, 44, 29};

  assert(jlayer > 0 && jlayer < 10);
  auto layer = static_cast<NAIA::TrTrack::FitPositionHeight>(jlayer - 1);

  constexpr std::array<float, 9> LayerPositionHeightZ = {
      158.920, 53.060, 29.228, 25.212, 1.698, -2.318, -25.212, -29.228, -135.882,
  };

  float x, y;
  if (event.header->CheckMask(NAIA::Category::HasTofStandalone)) {
    x = event.tofBaseSt->InterpolateAtZ(LayerPositionHeightZ.at(layer))[0];
    y = event.tofBaseSt->InterpolateAtZ(LayerPositionHeightZ.at(layer))[1];
  } else
    return false;

  if (layer == NAIA::TrTrack::FitPositionHeight::Layer9) {
    if (std::fabs(x) > fid_R[layer] || std::fabs(y) > fid_Y[layer])
      return false;
  } else {
    if (std::sqrt(x * x + y * y) > fid_R[layer] || std::fabs(y) > fid_Y[layer])
      return false;
  }

  return true;
}


bool IsInsideLayer(unsigned int jlayer, float margin, NAIA::Event &event) {
  bool isinlayer = false;

  assert(jlayer > 0 && jlayer < 10);
  auto layer = static_cast<NAIA::TrTrack::FitPositionHeight>(jlayer - 1);

  static float tracker_planes_edges[9][4] = {
      {-62.14, -47.40, 62.14, 47.40}, {-62.14, -40.10, 62.14, 40.10}, {-49.70, -43.75, 49.70, 43.75},
      {-49.72, -43.75, 49.72, 43.75}, {-49.71, -36.45, 49.70, 36.45}, {-49.72, -36.45, 49.72, 36.45},
      {-49.72, -43.75, 49.71, 43.75}, {-49.72, -43.75, 49.71, 43.75}, {-45.62, -29.48, 45.55, 29.53}};

  constexpr std::array<float, 9> LayerPositionHeightZ = {
      158.920, 53.060, 29.228, 25.212, 1.698, -2.318, -25.212, -29.228, -135.882,
  };

  float x, y;
  if (event.header->CheckMask(NAIA::Category::HasTofStandalone)) {
    x = event.tofBaseSt->InterpolateAtZ(LayerPositionHeightZ.at(layer))[0];
    y = event.tofBaseSt->InterpolateAtZ(LayerPositionHeightZ.at(layer))[1];
  } else
    isinlayer = false;

  if ((x > tracker_planes_edges[layer][0] + margin) && (x < tracker_planes_edges[layer][2] - margin) &&
      (y > tracker_planes_edges[layer][1] + margin) && (y < tracker_planes_edges[layer][3] - margin)) {

    if (layer == NAIA::TrTrack::FitPositionHeight::Layer9) {
      // remove dead area center L9
      if (x < -0.5 || y > 0.5) {
        isinlayer = true;
      }
    } else {
      // circle
      if ((sqrt(x * x + y * y) < tracker_planes_edges[layer][2] - margin))
        isinlayer = true;
    }
  }
  return isinlayer;
}

TrTrackEffSel::IsInsideL1::IsInsideL1(float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return IsInsideLayer(1, margin, event); });
}

TrTrackEffSel::IsInsideL9::IsInsideL9(float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return IsInsideLayer(9, margin, event); });
}

TrTrackEffSel::IsInsideInner::IsInsideInner(float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    if (IsInsideLayer(2, margin, event) && (IsInsideLayer(3, margin, event) || IsInsideLayer(4, margin, event)) &&
        (IsInsideLayer(5, margin, event) || IsInsideLayer(6, margin, event)) &&
        (IsInsideLayer(7, margin, event) || IsInsideLayer(8, margin, event)) && true)
      return true;
    else
      return false;
  });
}

TrTrackEffSel::InnerFiducialVolume::InnerFiducialVolume() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    int nhits = 0;
    for (auto layer = 2; layer < 9; layer++)
      if (IsInsideFiducial(layer, event))
        nhits++;
    if (nhits < 5)
      return false;

    if (!IsInsideFiducial(2, event))
      return false;
    if (!(IsInsideFiducial(3, event) || IsInsideFiducial(4, event)))
      return false;
    if (!(IsInsideFiducial(5, event) || IsInsideFiducial(6, event)))
      return false;
    if (!(IsInsideFiducial(7, event) || IsInsideFiducial(8, event)))
      return false;

    return true;
  });
}

TrTrackEffSel::L1FiducialVolume::L1FiducialVolume() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return IsInsideFiducial(1, event); });
}

TrTrackEffSel::L9FiducialVolume::L9FiducialVolume() {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) { return IsInsideFiducial(9, event); });
}
 

}//TrTrackEffSel
}//Efficinecy


namespace trg { //trg
class HasTrigger : public NSL::Selection {
public:
 HasTrigger(unsigned int triggerPatt) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    auto hasTrigger = false;
    if (event.evSummary->PhysBPatt == triggerPatt)
      hasTrigger = true;
    return (hasTrigger);
   });
  }
};
} //trg



namespace innTr { //innTr
class ChResoInnLessThan : public NSL::Selection{ 
public:
  ChResoInnLessThan(NAIA::TrTrack::ChargeRecoType charge, float max) {
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trTrackBase->InnerChargeRMS[charge]/(event.trTrackBase->InnerCharge[charge]) < max );
    });
  }
};

class L1NormResidualLessThan : public NSL::Selection {
public:
  L1NormResidualLessThan(float max, NAIA::TrTrack::Fit fit) {
  	m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
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
};

class NGoodClustersGreaterThan : public NSL::Selection {
public:
  NGoodClustersGreaterThan(unsigned int min, int mask) {
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
};

} //innTr


namespace TofSt{ // TofSt

class oneHitOnEachLayer : public NSL::Selection {
public:
	oneHitOnEachLayer() {
		m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
		if (!event.header->CheckMask(NAIA::Category::HasTofStandalone)) return false;
		return (event.tofPlusSt->NClusters[0][ON] > 0 && event.tofPlusSt->NClusters[1][ON] > 0 &&
			 event.tofPlusSt->NClusters[2][ON] > 0 && event.tofPlusSt->NClusters[3][ON] > 0 );
		});
	}	
};


class goodPathAll : public NSL::Selection{
public:
	goodPathAll() {
		m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
			return (event.tofPlusSt->LayerGoodPathl[0b1111] );
		});	
	}
};
class tofBetaInRange : public NSL::Selection{
public:
  tofBetaInRange(double min, double max, NAIA::Tof::BetaType type){
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      auto beta = event.tofBaseSt->Beta[type];
      return (beta > min && beta < max);
    });
  }
};
class tofChi2TimeLessThan : public NSL::Selection{
public:
  tofChi2TimeLessThan(double max){
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.tofPlusSt->Chi2Time < max);
    });
  }
};
class tofChi2CooLessThan : public NSL::Selection{
public:
  tofChi2CooLessThan(double max){
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.tofPlusSt->Chi2Coo < max);
    });
  }
};
class tofChargeInRange : public NSL::Selection{
public:
  tofChargeInRange(double min, double max, NAIA::Tof::ChargeType type){
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      auto charge = event.tofBaseSt->Charge[type];
      return (charge > min && charge < max);
    });
  }
};
bool tofIsInsideFiducial(NAIA::TrTrack::FitPositionHeight layer, NAIA::Event &event) {
  if (!event.header->CheckMask(NAIA::Category::HasTofStandalone)) return false;
  auto pos = event.tofBaseSt->InterpolateAtZ(NAIA::TrTrack::fitPositionHeightZ[layer]); //x,y
  if(sqrt(pos[0]*pos[0]+pos[1]*pos[1]) > fid_R[layer] || fabs(pos[1]) > fid_Y[layer]) return false;
  return true;
}
class tofInnerFiducialVolume : public NSL::Selection{
public:
  tofInnerFiducialVolume() {
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      short hits = 0;
      if(tofIsInsideFiducial(L2, event)) hits |= 1<<0;
      if(tofIsInsideFiducial(L3, event)) hits |= 1<<1;
      if(tofIsInsideFiducial(L4, event)) hits |= 1<<2;
      if(tofIsInsideFiducial(L5, event)) hits |= 1<<3;
      if(tofIsInsideFiducial(L6, event)) hits |= 1<<4;
      if(tofIsInsideFiducial(L7, event)) hits |= 1<<5;
      if(tofIsInsideFiducial(L8, event)) hits |= 1<<6;
      if(bitcounter(hits)     <5) return false; //check number of hits
      if(bitcounter(hits&1)  ==0) return false; //check L2 hit
      if(bitcounter(hits&6)  ==0) return false; //check L3||L4 hit
      if(bitcounter(hits&24) ==0) return false; //check L5||L6 hit
      if(bitcounter(hits&96) ==0) return false; //check L7||L8 hit
      return true;
    });
  }
};
class tofL1FiducialVolume : public NSL::Selection{
public:
  tofL1FiducialVolume(){
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      return tofIsInsideFiducial(L1, event);
    });
  }
};
} // TofSt



namespace trd { //TRD
class trdUcharge : public NSL::Selection{ 
public:
  trdUcharge(float min, float max) {
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trdKBase->Charge[TRD_UP] > min && event.trdKBase->Charge[TRD_UP] < max);
    });
  }
};
class trdLcharge : public NSL::Selection{
public:
  trdLcharge(float min, float max) {
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trdKBase->Charge[TRD_LOW] > min && event.trdKBase->Charge[TRD_LOW] < max);
    });
  }
};
class HasGoodTRDHits : public NSL::Selection{
public:
  HasGoodTRDHits(unsigned int charge) {
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
}
};

} //TRD


namespace TrdSt{ // TRD Standalone
class trdL1FiducialVolume : public NSL::Selection{
public:
  trdL1FiducialVolume(){
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      if(!event.CheckMask(NAIA::Category::HasTrdStandalone)) return false;
      auto pos = event.trdKBaseSt->InterpolateAtZ(NAIA::TrTrack::fitPositionHeightZ[L1]); //x,y
      if(sqrt(pos[0]*pos[0]+pos[1]*pos[1]) > fid_R[L1] || fabs(pos[1]) > fid_Y[L1]) return false;
      return true;
    });
  }
};

class trdNHitsAtleast : public NSL::Selection{
	public:
	trdNHitsAtleast(float min) {
		m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
			return (event.trdKBaseSt->NHits[TRD_OFFTRACK] + event.trdKBaseSt->NHits[TRD_ONTRACK] >= min);
		});	
	}
};

class trdUcharge : public NSL::Selection{ 
public:
  trdUcharge(float min, float max) {
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trdKBaseSt->Charge[TRD_UP] > min && event.trdKBaseSt->Charge[TRD_UP] < max);
    });
  }
};
class trdLcharge : public NSL::Selection{
public:
  trdLcharge(float min, float max) {
    m_matcher = make_shared<NSL::boolMatcher>([=](Event &event) {
      return (event.trdKBaseSt->Charge[TRD_LOW] > min && event.trdKBaseSt->Charge[TRD_LOW] < max);
    });
  }
};

} // TRD Standalone


//RTI Default selection
namespace NSL {
namespace Selections {

bool DefaultRTISelection(const NAIA::RTIInfo &rti_info) {
  if (rti_info.LivetimeFraction < 0.5 || rti_info.Zenith > 40 || rti_info.nTrigger / rti_info.nEvent < 0.98 ||
      rti_info.nError < 0 || rti_info.nError / rti_info.nEvent > 0.1 || rti_info.MeanAlignDiffExtLayer[0][1] > 35 ||
      rti_info.MeanAlignDiffExtLayer[1][1] > 45)
    return false;
  return true;
}
} // namespace Selections
} // namespace NSL



namespace Efficiency {
namespace TriggerEffSel {

bool IsInsideLayer(unsigned int jlayer, NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin,
                   NAIA::Event &event) {
  bool isinlayer = false;

  assert(jlayer > 0 && jlayer < 10);
  auto layer = static_cast<NAIA::TrTrack::FitPositionHeight>(jlayer - 1);

  static float tracker_planes_edges[9][4] = {
      {-62.14, -47.40, 62.14, 47.40}, {-62.14, -40.10, 62.14, 40.10}, {-49.70, -43.75, 49.70, 43.75},
      {-49.72, -43.75, 49.72, 43.75}, {-49.71, -36.45, 49.70, 36.45}, {-49.72, -36.45, 49.72, 36.45},
      {-49.72, -43.75, 49.71, 43.75}, {-49.72, -43.75, 49.71, 43.75}, {-45.62, -29.48, 45.55, 29.53}};

  float x = event.trTrackBase->TrTrackFitPos[layer][fitType][spanType][NAIA::TrTrack::Side::X];
  float y = event.trTrackBase->TrTrackFitPos[layer][fitType][spanType][NAIA::TrTrack::Side::Y];

  if ((x > tracker_planes_edges[layer][0] + margin) && (x < tracker_planes_edges[layer][2] - margin) &&
      (y > tracker_planes_edges[layer][1] + margin) && (y < tracker_planes_edges[layer][3] - margin)) {

    if (layer == NAIA::TrTrack::FitPositionHeight::Layer9) {
      // remove dead area center L9
      if (x < -0.5 || y > 0.5) {
        isinlayer = true;
      }
    } else {
      // circle
      if ((sqrt(x * x + y * y) < tracker_planes_edges[layer][2] - margin))
        isinlayer = true;
    }
  }
  return isinlayer;
}

class HasTrigger : public NSL::Selection {
public:
  HasTrigger(unsigned int triggerPatt);
};

class IsInsideL1 : public NSL::Selection {
public:
  IsInsideL1(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin);
};

class IsInsideL9 : public NSL::Selection {
public:
  IsInsideL9(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin);
};

TriggerEffSel::IsInsideL1::IsInsideL1(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    if (!event.trTrackBase->FitIDExists(fitType, spanType))
      return false;

    return IsInsideLayer(1, fitType, spanType, margin, event);
  });
}

} // namespace TriggerEffSel
} // namespace Efficiency

namespace ms { //ms

bool IsInsideLayer(unsigned int jlayer, NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin,
                   NAIA::Event &event) {
  bool isinlayer = false;

  assert(jlayer > 0 && jlayer < 10);
  auto layer = static_cast<NAIA::TrTrack::FitPositionHeight>(jlayer - 1);

  static float tracker_planes_edges[9][4] = {
      {-62.14, -47.40, 62.14, 47.40}, {-62.14, -40.10, 62.14, 40.10}, {-49.70, -43.75, 49.70, 43.75},
      {-49.72, -43.75, 49.72, 43.75}, {-49.71, -36.45, 49.70, 36.45}, {-49.72, -36.45, 49.72, 36.45},
      {-49.72, -43.75, 49.71, 43.75}, {-49.72, -43.75, 49.71, 43.75}, {-45.62, -29.48, 45.55, 29.53}};

  float x = event.trTrackBase->TrTrackFitPos[layer][fitType][spanType][NAIA::TrTrack::Side::X];
  float y = event.trTrackBase->TrTrackFitPos[layer][fitType][spanType][NAIA::TrTrack::Side::Y];

  if ((x > tracker_planes_edges[layer][0] + margin) && (x < tracker_planes_edges[layer][2] - margin) &&
      (y > tracker_planes_edges[layer][1] + margin) && (y < tracker_planes_edges[layer][3] - margin)) {

    if (layer == NAIA::TrTrack::FitPositionHeight::Layer9) {
      // remove dead area center L9
      if (x < -0.5 || y > 0.5) {
        isinlayer = true;
      }
    } else {
      // circle
      if ((sqrt(x * x + y * y) < tracker_planes_edges[layer][2] - margin))
        isinlayer = true;
    }
  }
  return isinlayer;
}

class IsInsideL1 : public NSL::Selection {
public:
  IsInsideL1(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType, float margin) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    if (!event.trTrackBase->FitIDExists(fitType, spanType))
      return false;
    return IsInsideLayer(1, fitType, spanType, margin, event);
  });
}
};

class HasGoodSecondTrTrack : public NSL::Selection {
public:
  HasGoodSecondTrTrack(NAIA::TrTrack::Fit fitType, NAIA::TrTrack::Span spanType,
                                                      float min) {
  m_matcher = std::make_shared<NSL::boolMatcher>([=](Event &event) {
    bool good_sec_track = true;
    bool isAboveThreshold = false;
    bool hasMinumumHitsY = false;
    if (event.secondTrTrackBase->FitIDExists(fitType, spanType)) {
      if (event.secondTrTrackBase->Rigidity[fitType][spanType] > min)
        isAboveThreshold = true;
      unsigned int nHits = 0;
      for (unsigned int i = 1; i < 8; ++i) {
        if (NAIA::ContainsKeys(event.secondTrTrackBase->TrTrackHitPos, i, NAIA::TrTrack::Side::Y))
          ++nHits;
      }
      if (nHits > 3)
        hasMinumumHitsY = true;
    }
    if (isAboveThreshold && hasMinumumHitsY)
      good_sec_track = false;
    return good_sec_track;
  });
}
};

}//ms
