#include "selection.h"

TH1D* GetEff(TH1* hpass, TH1* htotal) {
  hpass->Sumw2();
  htotal->Sumw2();
  auto ratio = new TH1D("", "", nRbins-1, Rbins);
  auto eff = new TEfficiency(*hpass, *htotal);
  eff->SetStatisticOption(TEfficiency::kBJeffrey);
  
  for(int i=1; i<ratio->GetNbinsX(); ++i) {
	  ratio->SetBinContent(i, eff->GetEfficiency(i));
      double hierr = eff->GetEfficiencyErrorUp(i);
      double loerr = eff->GetEfficiencyErrorLow(i);
      double err = hierr>loerr ? hierr : loerr;
      ratio->SetBinError(i, err);
    }
  return ratio;
}
