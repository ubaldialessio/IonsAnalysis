#include "binning.h"
#include "utils.h"

template <typename T>
std::vector<T> getHistogram(std::vector<TFile*> v,std::vector<TString>t,bool normalize = false);

std::vector<TString> histogram = {"qInn_R_den","qInn_R_den_ltof","qInn_R_den_ltof_l9_05","qInn_R_den_ltof_l9_1"};
std::map<TString,TString> name = { {"qInn_R_den","Den"},
                                   {"qInn_R_den_ltof","Den+LTof"},
                                   {"qInn_R_den_ltof_l9_05","Den+LTof+UnbL9"},
                                   {"qInn_R_den_ltof_l9_1","Den+LTof+UnbL9>"} };
std::vector<TString> effHist      = {"sample_trch_den","pass_trch_den",
                                     "sample_trch_den_ltof","pass_trch_den_ltof",
                                     "sample_trch_den_ltof_l9_05","pass_trch_den_ltof_l9_05",
                                     "sample_trch_den_ltof_l9_1","pass_trch_den_ltof_l9_1"};
std::map<TString,TString> effName= { {"qInn_R_den","Den"},
                                   {"qInn_R_den_ltof","Den+LTof"},
                                   {"qInn_R_den_ltof_l9_05","Den+LTof+UnbL9"},
                                   {"qInn_R_den_ltof_l9_1","Den+LTof+UnbL9>"} };                                     

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("Usage: %s <mode: dat|mc|both> <charge1> <charge2> ... <charge n> [--norm]\n", argv[0]);
        return 1;
    }

    TString mode = argv[1];
    bool normalize = false;
    for (int i = 1; i < argc; ++i) {
        if (TString(argv[i]) == "--norm") normalize = true;
    }
    std::vector<unsigned int> charge_number_dat, charge_number_mc;
    std::vector<TFile*> dataFiles, mcFiles;
    std::vector<TH1D*> dataHist, mcHist;
    std::vector<TH2D*> dataEff, mcEff;

    if (mode=="dat" || mode=="both") {
        dataFiles = getDistributionFiles(argc,argv,charge_number_dat,"data","TrackCharge");
        dataHist = getHistogram<TH1D*>(dataFiles,histogram,normalize);
    }
    if (mode=="mc" || mode=="both") {
        mcFiles   = getDistributionFiles(argc,argv,charge_number_mc,"mc","TrackCharge");
        mcHist   = getHistogram<TH1D*>(mcFiles,histogram,normalize);
    }

    //auto dataEff  = getHistogram<TH2D*>(dataFiles,effHist);
    //auto mcEff    = getHistogram<TH2D*>(mcFiles,effHist);

    if (mode == "dat") {
        printDistributions(dataHist, {}, charge_number_dat, "innTrk", histogram, name);
    } else if (mode == "mc") {
        printDistributions({}, mcHist, charge_number_mc, "innTrk", histogram, name);
    } else if (mode == "both") {
        printDistributions(dataHist, mcHist, charge_number_dat, "innTrk", histogram, name);
    } else {
        printf("Unknown mode '%s'. Use dat, mc, or both.\n", argv[1]);
        return 1;
    }

    return 0;
}

template <typename T>
std::vector<T> getHistogram(std::vector<TFile*> v, std::vector<TString> t, bool normalize) {
    std::vector<T> vec;
    for (size_t i = 0; i < v.size(); i++) {
        for (const auto &hist : t) {
            if (!v[i]) {
                printf("Warning: null file pointer for file index %zu\n", i);
                vec.push_back((T)nullptr);
                continue;
            }
            auto h = (T)v[i]->Get(Form("IL1/%s", hist.Data()));
            if (h) {
                if (normalize && h->Integral() > 0) {
                    h->Scale(1.0 / h->Integral());
                }
                vec.push_back(h);
            } else {
                printf("Warning: histogram '%s' not found in file %s (file idx %zu). Pushing nullptr.\n",
                       hist.Data(), v[i]->GetName(), i);
                vec.push_back((T)nullptr);
            }
        }
    }
    return vec;
}