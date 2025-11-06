#ifndef CHARGE_ANALYSIS_MANAGER_H
#define CHARGE_ANALYSIS_MANAGER_H

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
/*
#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"*/

double FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain,  NAIA::Event &event); 

/// Contesto per ogni carica (Amass, array DAQ, selezioni)
struct ChargeContext {
    double Amass;
    std::array<unsigned int, 24> daq;
    std::array<double, 8> daq_event_size;

    std::function<bool( NAIA::Event&)> mysel;
    std::function<bool( NAIA::Event&)> myselNoL1;

    std::function<bool( NAIA::Event&)> num_tof;
    std::function<bool( NAIA::Event&)> den_tof;
    std::function<bool( NAIA::Event&)> den_tof_ltof;
    std::function<bool( NAIA::Event&)> den_tof_ltof_l9_05;
    std::function<bool( NAIA::Event&)> den_tof_ltof_l9_1;

    std::function<bool( NAIA::Event&)> num_l1;
    std::function<bool( NAIA::Event&)> den_l1;
    std::function<bool( NAIA::Event&)> den_l1_ltof;
    std::function<bool( NAIA::Event&)> den_l1_ltof_l9_05;
    std::function<bool( NAIA::Event&)> den_l1_ltof_l9_1;

    std::function<bool( NAIA::Event&)> num_l1Unb;
    std::function<bool( NAIA::Event&)> den_l1Unb;
    std::function<bool( NAIA::Event&)> den_l1Unb_ltof;
    std::function<bool( NAIA::Event&)> den_l1Unb_ltof_l9_05;
    std::function<bool( NAIA::Event&)> den_l1Unb_ltof_l9_1;

    std::function<bool( NAIA::Event&)> num_track;
    std::function<bool( NAIA::Event&)> den_track;

    std::function<bool( NAIA::Event&)> num_trackCh;
    std::function<bool( NAIA::Event&)> den_trackCh;
    std::function<bool( NAIA::Event&)> den_trackCh_ltof;
    std::function<bool( NAIA::Event&)> den_trackCh_ltof_l9_05;
    std::function<bool( NAIA::Event&)> den_trackCh_ltof_l9_1;

    std::function<bool( NAIA::Event&)> den_trig;

    std::function<bool( NAIA::Event&)> l1_temp;
    std::function<bool( NAIA::Event&)> l2_temp;
};
// -------- TemplateSet helper --------
struct TemplateSet {
    // L1 templates
    std::vector<TH1D*> L1Templates;
    std::vector<TTree*> L1KernelTrees;
    std::vector<float>  L1charge; // must live for branches
    std::vector<double> L1weight;

    // L2 templates
    std::vector<TH1D*> L2Templates;
    std::vector<TTree*> L2KernelTrees;
    std::vector<float>  L2charge;
    std::vector<double> L2weight;

    // L1 charge distributions
    std::vector<TH1D*> L1Distributions;
    std::vector<TTree*> L1DistKernelTrees;
    std::vector<float>  L1Dcharge;
    std::vector<double> L1Dweight;

    void Clear() {
        for (auto p : L1Templates)    { if (p) delete p; }
        for (auto p : L2Templates)    { if (p) delete p; }
        for (auto p : L1Distributions){ if (p) delete p; }
        for (auto t : L1KernelTrees)  { if (t) delete t; }
        for (auto t : L2KernelTrees)  { if (t) delete t; }
        for (auto t : L1DistKernelTrees) { if (t) delete t; }

        L1Templates.clear(); L1KernelTrees.clear(); L1charge.clear(); L1weight.clear();
        L2Templates.clear(); L2KernelTrees.clear(); L2charge.clear(); L2weight.clear();
        L1Distributions.clear(); L1DistKernelTrees.clear(); L1Dcharge.clear(); L1Dweight.clear();
    }
};
struct McSample { //ORIGINALE
  double gen;
  unsigned int species;
  short passed;
  float reco_inn;
  float reco_il1;
  float reco_beta;
  void Reset() { gen = std::numeric_limits<double>::quiet_NaN(); species = 0; passed = 0; reco_inn = -999.0f; reco_il1 =-999.0f;  reco_beta = -999.0f; }
};

class ChargeAnalysisManager {
public:
    ChargeAnalysisManager(const std::vector<unsigned int>& charges,
                          int nTbins, const double* Tbins,
                          int nRbins, const double* Rbins,
                          bool isMC, TString exename);

    // Template filling helpers
    void FillTemplateL1(unsigned int charge, const std::string &geom, int ibin, float q, double w);
    void FillTemplateL2(unsigned int charge, const std::string &geom, int ibin, float q, double w);
    void FillTemplateL1Distribution(unsigned int charge, const std::string &geom, int ibin, float q, double w);

    // Write template distributions/trees to file
    void WriteTemplateDistributions(const std::map<unsigned int, std::string>& ionPaths,
                                    const std::string& outname,
                                    bool isMC);


    void Fill(unsigned int charge,
              const std::string& geom,
              const std::string& name,
              double x, double y, double weight);

    void Fill(unsigned int charge,
                                 const std::string& geom,
                                 const std::string& name,
                                 double x, double y);

    void FillDistributions(unsigned int charge,
                                 const std::string& geom,
                                 const std::string& name,
                                 double x);
    
    void FillMcTrees(NAIA::NAIAChain &chain, NAIA::Event &event, const std::string &geom);

    void FillMcSampleFromFileInfo(NAIA::NAIAChain& chain, const std::string &geom);

    void BookHistograms(unsigned int charge);

    void BookDistributionHistograms(unsigned int charge);

    void BookTemplateHistograms(unsigned int charge);

    void BookMcTrees(unsigned int charge);

    void BookMcHist(unsigned int charge);

    void Write(const std::map<unsigned int, std::string>& ionPaths,
               const std::string& outname,
               bool isMC,
               unsigned int charge);

    void WriteDistributions(const std::map<unsigned int, std::string>& ionPaths,
                                  const std::string& outname,
                                  bool isMC);
    
    void WriteMcTrees(const std::map<unsigned int, std::string>& ionPaths,
                                               const std::string& outname,
                                               bool isMC);

    TH2D* GetHist(unsigned int charge,
                  const std::string& geom,
                  const std::string& name);

    /// Accesso al contesto per una certa carica
    ChargeContext& GetContext(unsigned int charge) { return contexts[charge]; }

private:
    std::vector<unsigned int> fCharges;
    int fNTbins;
    int fNRbins;
    std::vector<double> fTbins;
    std::vector<double> fRbins;
    bool fIsMC = false;   
    TString fExename;

    // Map (charge → geom → name → histo)
    std::map<unsigned int, 
             std::map<std::string, 
                      std::map<std::string, TH2D*>>> histos;

    // Map of distribution histograms
    std::map<unsigned int, 
             std::map<std::string, 
                      std::map<std::string, TH1D*>>> distributionHistos;   

    // Map of template histograms (charge->geom->TemplateSet)
    std::map<unsigned int, std::map<std::string, TemplateSet>> templateSets;

    // Map of Mc Trees            
    std::map<unsigned int, 
             std::map<std::string, 
                      std::map<std::string, TTree*>>> mcTrees;   

    // Map of Mc Samples           
    std::map<unsigned int, 
             std::map<std::string, 
                      std::map<std::string, McSample>>> mcSamples;           
    // Map of Mc Histograms
    std::map<unsigned int, 
             std::map<std::string, 
                      std::map<std::string, TH1D*>>> mcHist; 

    std::map<unsigned int, ChargeContext> contexts;
};

#endif
