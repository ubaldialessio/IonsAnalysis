#include "RTI.h"
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
