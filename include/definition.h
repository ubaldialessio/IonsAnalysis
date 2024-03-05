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

} // namespace TriggerEffSel
} // namespace Efficiency
