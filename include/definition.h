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

using namespace std;
namespace ns = NSL::Selections;
int fov = 1; // Fields of view: 0 = 25 deg, 1 = 30 deg, 2 = 35 deg, 3 = 40 deg

auto CRT = NAIA::TrTrack::ChargeRecoType::YJ;
auto FIT = NAIA::TrTrack::Fit::GBL;
auto INN = NAIA::TrTrack::Span::InnerOnly;
auto IL1 = NAIA::TrTrack::Span::InnerL1;
auto YSD = NAIA::TrTrack::Side::Y;
auto XSD = NAIA::TrTrack::Side::X;
auto UTC = NAIA::Tof::ChargeType::Upper;
auto LTC = NAIA::Tof::ChargeType::Lower;
auto BTH = NAIA::Tof::BetaType::BetaH;
auto OT  = NAIA::Tof::OnTime;
auto L1  = NAIA::TrTrack::FitPositionHeight::Layer1;
auto L2  = NAIA::TrTrack::FitPositionHeight::Layer2;
auto L3  = NAIA::TrTrack::FitPositionHeight::Layer3;
auto L4  = NAIA::TrTrack::FitPositionHeight::Layer4;
auto L5  = NAIA::TrTrack::FitPositionHeight::Layer5;
auto L6  = NAIA::TrTrack::FitPositionHeight::Layer6;
auto L7  = NAIA::TrTrack::FitPositionHeight::Layer7;
auto L8  = NAIA::TrTrack::FitPositionHeight::Layer8;
auto EL1 = NAIA::UnbExtHitBaseData::ExtHit::L1;
