#include "binning.h"
#include "utils.h"


template <typename T>
std::vector<T> getHistogram(std::vector<TFile*> v,std::vector<TString> t);

std::vector<TString> histogram = {"qL1Unb_den","qL1Unb_den_ltof","qL1Unb_den_ltof_l9_05","qL1Unb_den_ltof_l9_1"};
std::map<TString,TString> name = { {"qL1Unb_den","Den"},
                                   {"qL1Unb_den_ltof","Den+LTof"},
                                   {"qL1Unb_den_ltof_l9_05","Den+LTof+UnbL9"},
                                   {"qL1Unb_den_ltof_l9_1","Den+LTof+UnbL9>"} };
std::vector<TString> effHist      = {"sample_l1Unb_den","pass_l1Unb_den",
                                     "sample_l1Unb_den_ltof","pass_l1Unb_den_ltof",
                                     "sample_l1Unb_den_ltof_l9_05","pass_l1Unb_den_ltof_l9_05",
                                     "sample_l1Unb_den_ltof_l9_1","pass_l1Unb_den_ltof_l9_1"};
std::map<TString,TString> effName= { {"qL1Unb_den","Den"},
                                   {"qL1Unb_den_ltof","Den+LTof"},
                                   {"qL1Unb_den_ltof_l9_05","Den+LTof+UnbL9"},
                                   {"qL1Unb_den_ltof_l9_1","Den+LTof+UnbL9>"} };                                     

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("Usage: %s <mode: dat|mc|both> <charge1> <charge2> ... <charge n>\n", argv[0]);
        return 1;
    }

    TString mode = argv[1];
    std::vector<unsigned int> charge_number_dat, charge_number_mc;
    std::vector<TFile*> dataFiles, mcFiles;
    std::vector<TH1D*> dataHist, mcHist;

    if (mode=="dat" || mode=="both") {
        dataFiles = getDistributionFiles(argc,argv,charge_number_dat,"data","L1Unb");
        dataHist = getHistogram<TH1D*>(dataFiles,histogram);
    }
    if (mode=="mc" || mode=="both") {
        mcFiles = getDistributionFiles(argc,argv,charge_number_mc,"mc","L1Unb");
        mcHist  = getHistogram<TH1D*>(mcFiles,histogram);
    }

    if (mode == "dat") {
        printDistributions(dataHist, {}, charge_number_dat, "l1unb", histogram, name);
    } else if (mode == "mc") {
        printDistributions({}, mcHist, charge_number_mc, "l1unb", histogram, name);
    } else if (mode == "both") {
        printDistributions(dataHist, mcHist, charge_number_dat, "l1unb", histogram, name);
    } else {
        printf("Unknown mode '%s'. Use dat, mc, or both.\n", argv[1]);
        return 1;
    }
    return 0;
}

template <typename T>
std::vector<T> getHistogram(std::vector<TFile*> v,std::vector<TString>t) {
    std::vector<T> vec;
    for (size_t i=0; i<v.size(); i++) {
        for (const auto &hist : t) {
            if (!v[i]) { vec.push_back((T)nullptr); continue; }
            auto h = (T)v[i]->Get(("IL1/"+hist).Data());
            if (h) vec.push_back(h);
            else {
                printf("Warning: histogram '%s' not found in file %s\n",
                       hist.Data(), v[i]->GetName());
                vec.push_back((T)nullptr);
            }
        }
    }
    return vec;
}