#ifndef TRD_H
#define TRD_H
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
namespace trd { //TRD

class trdUcharge : public NSL::Selection{ 
public:
  trdUcharge(float min, float max);
};
  
class trdLcharge : public NSL::Selection{
public:
  trdLcharge(float min, float max);
};
  
class HasGoodTRDHits : public NSL::Selection{
public:
  HasGoodTRDHits(unsigned int charge);
};
  
} //TRD
}

namespace Efficiency {
namespace TrdSt{ // TRD Standalone

class trdL1FiducialVolume : public NSL::Selection{
public:
  trdL1FiducialVolume();
};
  
class trdNHitsAtleast : public NSL::Selection{
	public:
	trdNHitsAtleast(float min);
};
	
class trdUcharge : public NSL::Selection{ 
public:
  trdUcharge(float min, float max);
};
  
class trdLcharge : public NSL::Selection{
public:
  trdLcharge(float min, float max);
};

} // TRD Standalone
} //Efficiency
#endif
