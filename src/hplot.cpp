#include "binning.h"
#include "TLegend.h"
#include "utils.h"

void niceth2(TH2 *h);
void niceth1(TH1 *h, TString name, TString opt, unsigned int charge);
void print(TKey *key, TIter keyList, TString a, TString b, TString output, TString opt, unsigned int charge);


int main(int argc, char **argv) {
	if (argc < 2) {
		printf("Usage: \n");
		printf("%s <charge> <single/same> \n", argv[0]);
		return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString opt = argv[2];
	TString ionPath=getIonPath(charge);
	TString input  = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/"+ionPath+"/"+ionPath+"_output.root";
	TString output = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/output/"+ionPath+".pdf";
	TString a = output + "[", b = output + "]";
	TFile *f = new TFile(input.Data(),"read");
	TIter keyList(f->GetListOfKeys());
	TKey *key;
	//bool isMC = name_out.Contains("mc");
	auto style = effHistStyle(charge);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();
	print(key, keyList, a, b, output, opt,charge);
}

void niceth2(TH2 *h) {
	
}

void niceth1(TH1 *h, TString name, TString opt, unsigned int charge) {
	TString ion = getIonName(charge);
	TString TriggerEff = ion+" Trigger";
	TString L1Eff =ion+" L1";
	TString TofEff =ion+" Tof";
	TString TrackEff =ion+" Track";
	h->GetYaxis()->SetTitle("");
	if (opt.Contains("same"))
		h->Draw("SAME");
	if (opt.Contains("single"))
		h->Draw();
	h->GetXaxis()->SetTitle("R (GV)");
	if (name.Contains("d_")) {
		h->SetMarkerStyle(8);
		if (name.Contains("trig")) {
			h->GetYaxis()->SetRangeUser(0.5,1);
			h->SetTitle(TriggerEff.Data() );
			h->GetXaxis()->SetRangeUser(0,1000);
		} else if (name.Contains("l1")) {
			h->GetYaxis()->SetRangeUser(0.64,0.99);
			h->SetTitle(L1Eff.Data() );
			h->GetXaxis()->SetRangeUser(0,100);
		} else if (name.Contains("tof")) {
			h->GetYaxis()->SetRangeUser(0.80,1);
			h->SetTitle(TofEff.Data() );
			h->GetXaxis()->SetRangeUser(0,1000);
		} else if (name.Contains("tr") && !name.Contains("ig") ) {
			h->GetYaxis()->SetRangeUser(0.3,0.99);
			h->SetTitle(TrackEff.Data() );
			h->GetXaxis()->SetRangeUser(0,30);
		}
	}
	if (name.Contains("m_")) {
		h->SetMarkerStyle(32);
		if (name.Contains("trig")) {
			h->GetYaxis()->SetRangeUser(0.5,1);
			h->SetTitle(TriggerEff.Data() );
			h->GetXaxis()->SetRangeUser(0,1000);
		} else if (name.Contains("l1")) {
			h->GetYaxis()->SetRangeUser(0.6,0.9);
			h->SetTitle(L1Eff.Data() );
			h->GetXaxis()->SetRangeUser(0,100);
		} else if (name.Contains("tof")) {
			h->GetYaxis()->SetRangeUser(0.80,1);
			h->SetTitle(TofEff.Data() );
			h->GetXaxis()->SetRangeUser(0,1000);
		} else if (name.Contains("tr") && !name.Contains("ig") ) {
			h->GetYaxis()->SetRangeUser(0.1,0.99);
			h->SetTitle(TrackEff.Data() );
			h->GetXaxis()->SetRangeUser(0,30);
		}
	}
	
		if (name.Contains("lt")) {
			h->SetTitle("Livetime");
		} else if (name.Contains("counts")) {
			h->SetTitle("Counts");
		} else if (name.Contains("flux")) {
			h->SetTitle("Flux");
		} else if (name.Contains("old")) {
			h->SetTitle("OldFlux");
		}
		if (name.Contains("true_acc")) {
			h->SetTitle("True_Acc");
			h->GetXaxis()->SetRangeUser(0,3000);
		  	h->GetYaxis()->SetRangeUser(0,0.1);
		}
}

void print(TKey *key, TIter keyList, TString a, TString b, TString output, TString opt, unsigned int charge) {
	if (opt.Contains("single")) {
		TCanvas c1;
		c1.Print(a.Data() );
		while ((key = (TKey*)keyList())) {
			TClass *cl = gROOT->GetClass(key->GetClassName());
			if (cl->InheritsFrom("TH2")) {
				TH2 *h = (TH2*)key->ReadObj();
			}
			if (cl->InheritsFrom("TH1")) {
				gStyle->SetOptStat(0);
			   	TH1 *h = (TH1*)key->ReadObj();
			   	TString name = key->GetName();
			   	niceth1(h, name, opt, charge);
			   	if (!name.Contains("eff"))
			   		c1.SetLogy();
			}
			c1.SetLogx();
			c1.SetGridx();
			c1.SetGridy();
			c1.Print(output.Data() );
		}
		c1.Print(b.Data() );
	}
	
	if (opt.Contains("same")) {
		std::vector<TH1*> trig, tof, l1, tr;
		TCanvas c1,c2,c3,c4;
		TLegend *a1 = new TLegend(0.55,0.975,0.65,0.925);
		TLegend *a2 = new TLegend( 0.55,0.975,0.65,0.925);
		TLegend *a3 = new TLegend( 0.55,0.975,0.65,0.925);
		TLegend *a4 = new TLegend( 0.55,0.975,0.65,0.925);
		gROOT->ForceStyle();
		c1.SetLogx();
		c1.SetGridx();
		c1.SetGridy();
		c2.SetLogx();
		c2.SetGridx();
		c2.SetGridy();
		c3.SetLogx();
		c3.SetGridx();
		c3.SetGridy();
		c4.SetLogx();
		c4.SetGridx();
		c4.SetGridy();
		
		c1.Print(a.Data() );
		while ((key = (TKey*)keyList())) {
			TClass *cl = gROOT->GetClass(key->GetClassName());
			//if (cl->InheritsFrom("TH2"))  TH2 *h = (TH2*)key->ReadObj();
			if (cl->InheritsFrom("TH1")) {
				TH1 *h = (TH1*)key->ReadObj();
				h->GetYaxis()->SetTitle("");
				TString name = key->GetName();
				if (name.Contains("dataEff_tr") || name.Contains("mcEff_tr") )
					trig.push_back(h);
				if (name.Contains("dataEff_tf") || name.Contains("mcEff_tf") )
					tof.push_back(h);
				if (name.Contains("dataEff_l1") || name.Contains("mcEff_l1") )
					l1.push_back(h);
				if (name.Contains("dataEff_tk") || name.Contains("mcEff_tk") )
					tr.push_back(h);
			}
		}
		c1.cd();
		gROOT->ForceStyle();
		niceth1(trig[0],"d_trig",opt,charge);
		niceth1(trig[1],"m_trig",opt,charge);
		a1->AddEntry(trig[0],"Data");
		a1->AddEntry(trig[1],"MC");
		a1->Draw();
		c1.Print(output.Data() );
		c2.cd();
		niceth1(tof[0],"d_tof",opt,charge);
		niceth1(tof[1],"m_tof",opt,charge);
		a2->AddEntry(tof[0],"Data");
		a2->AddEntry(tof[1],"MC");
		a2->Draw();
		c2.Print(output.Data() );
		c3.cd();
		niceth1(l1[0],"d_l1",opt,charge);
		niceth1(l1[1],"m_l1",opt,charge);
		a3->AddEntry(l1[0],"Data");
		a3->AddEntry(l1[1],"MC");
		a3->Draw();
		c3.Print(output.Data() );
		c4.cd();
		niceth1(tr[0],"d_tr",opt,charge);
		niceth1(tr[1],"m_tr",opt,charge);
		a4->AddEntry(tr[0],"Data");
		a4->AddEntry(tr[1],"MC");
		a4->Draw();
		c4.Print(output.Data() );
		c4.Print(b.Data() );
	}
}
