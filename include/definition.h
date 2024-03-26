#ifndef DEFINITION_H
#define DEFINITION_H
#include "includes.h"

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

namespace trig 			= Efficiency::TriggerEff;
namespace track 		= Efficiency::TrTrackEffSel;
namespace InnerTracker  = Efficiency::InnerTracker;
namespace TofSt 		= Efficiency::TofStandalone;
namespace MySel 	    = Efficiency::MySel;
namespace trd 		    = Efficiency::trd;
#endif
