#include "selection.h"

void niceth2(TH2 *h) {
	
}

void niceth1(TH1 *h, TString name) {
	h->Draw();
	h->GetXaxis()->SetTitle("R (GV)");
	h->GetYaxis()->SetTitle("Efficiency");
	h->SetMarkerStyle(4);
	h->SetTitle(name.Data() );
	if (name == "trig_eff") {
		h->GetYaxis()->SetRangeUser(0.99,1);
		h->GetXaxis()->SetRangeUser(0,1000);
	} else if (name == "l1_eff") {
		h->GetYaxis()->SetRangeUser(0.5,1);
		h->GetXaxis()->SetRangeUser(0,100);
	} else if (name == "tof_eff") {
		h->GetYaxis()->SetRangeUser(0.8,1);
		h->GetXaxis()->SetRangeUser(0,1000);
	} else if (name == "tr_eff") {
		h->GetYaxis()->SetRangeUser(0.7,1);
		h->GetXaxis()->SetRangeUser(0,50);
	}
}

void print(TKey *key, TIter keyList, TString a, TString b, TString output) {
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
		   	string name = key->GetName();
		   	niceth1(h, name);
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
		printf("%s <input.root> <output.pdf> \n", argv[0]);
		return 1;
	}
	TString input = argv[1], output = argv[2];
	TString a = output + "[", b = output + "]";
	TFile *f = new TFile(input.Data(),"read");
	TIter keyList(f->GetListOfKeys());
	TKey *key;
	gStyle->SetErrorX(0);

	print(key, keyList, a, b, output);
}
