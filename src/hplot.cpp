#include "binning.h"

void niceth2(TH2 *h) {
	
}

void niceth1(TH1 *h, TString name) {
	h->Draw();
	h->GetXaxis()->SetTitle("R (GV)");
	h->GetYaxis()->SetTitle("Efficiency");
	h->SetMarkerStyle(20);
	h->SetTitle(name.Data() );
	if (name.Contains("trig")) {
		h->GetYaxis()->SetRangeUser(0.99,1);
		h->GetXaxis()->SetRangeUser(0,1000);
	} else if (name.Contains("l1")) {
		h->GetYaxis()->SetRangeUser(0.5,1);
		h->GetXaxis()->SetRangeUser(0,100);
	} else if (name.Contains("tof")) {
		h->GetYaxis()->SetRangeUser(0.8,1);
		h->GetXaxis()->SetRangeUser(0,1000);
	} else if (name.Contains("tr") && !name.Contains("ig")) {
		h->GetYaxis()->SetRangeUser(0.7,1);
		h->GetXaxis()->SetRangeUser(0,50);
	} else if (name.Contains("lt")) {
		h->GetYaxis()->SetTitle("Livetime");
	} else if (name.Contains("counts")) {
		h->GetYaxis()->SetTitle("Counts");
	} else if (name.Contains("flux")) {
		h->GetYaxis()->SetTitle("Flux");
	}
	if (name.Contains("true_acc")) {
		h->GetYaxis()->SetTitle("True_Acc");
		h->GetXaxis()->SetRangeUser(0,3000);
	  	h->GetYaxis()->SetRangeUser(0,0.1);
	}
}

void print(TKey *key, TIter keyList, TString a, TString b, TString output, bool isMC) {
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
		   	if (isMC) {	name+="_mc";
		   	} else if (!isMC && name.Contains("eff")) { name +="_data"; }
		   	niceth1(h, name);
		   	if (!name.Contains("eff") && !isMC)
		   		c1.SetLogy();
		}
		c1.SetLogx();
		c1.SetGridx();
		c1.SetGridy();
		c1.Print(output.Data() );
	}
	c1.Print(b.Data() );
}

int main(int argc, char **argv) {
	if (argc < 2) {
		printf("Usage: \n");
		printf("%s <input.root> <output_name.pdf> \n", argv[0]);
		return 1;
	}
	TString input = argv[1], name_out = argv[2], path_out = "/storage/gpfs_ams/ams/users/aubaldi/Data/eff/", output = path_out+name_out;
	TString a = output + "[", b = output + "]";
	TFile *f = new TFile(input.Data(),"read");
	TIter keyList(f->GetListOfKeys());
	TKey *key;
	bool isMC = name_out.Contains("mc");
	gStyle->SetErrorX(0);

	print(key, keyList, a, b, output, isMC);
}
