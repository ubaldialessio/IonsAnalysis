#include "binning.h"
#include "utils.h"

template <typename T>
std::vector<T> getHistogram(std::vector<TFile*> v,std::vector<TString> t);

std::vector<TString> histogram = {"qtof_den","qtof_den_ltof","qtof_den_ltof_l9_05","qtof_den_ltof_l9_1"};
std::map<TString,TString> name = { {"qtof_den","Den"},
                                   {"qtof_den_ltof","Den+LTof"},
                                   {"qtof_den_ltof_l9_05","Den+LTof+UnbL9"},
                                   {"qtof_den_ltof_l9_1","Den+LTof+UnbL9>"} };
std::vector<TString> effHist      = {"sample_tof_den","pass_tof_den",
                                     "sample_tof_den_ltof","pass_tof_den_ltof",
                                     "sample_tof_den_ltof_l9_05","pass_tof_den_ltof_l9_05",
                                     "sample_tof_den_ltof_l9_1","pass_tof_den_ltof_l9_1"};
std::map<TString,TString> effName= { {"qtof_den","Den"},
                                   {"qtof_den_ltof","Den+LTof"},
                                   {"qtof_den_ltof_l9_05","Den+LTof+UnbL9"},
                                   {"qtof_den_ltof_l9_1","Den+LTof+UnbL9>"} };                                     

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
        dataFiles = getDistributionFiles(argc,argv,charge_number_dat,"data","Tof");
        dataHist = getHistogram<TH1D*>(dataFiles,histogram);
    }
    if (mode=="mc" || mode=="both") {
        mcFiles = getDistributionFiles(argc,argv,charge_number_mc,"mc","Tof");
        mcHist  = getHistogram<TH1D*>(mcFiles,histogram);
    }

    if (mode == "dat") {
        printDistributions(dataHist, {}, charge_number_dat, "tof", histogram, name);
    } else if (mode == "mc") {
        printDistributions({}, mcHist, charge_number_mc, "tof", histogram, name);
    } else if (mode == "both") {
        printDistributions(dataHist, mcHist, charge_number_dat, "tof", histogram, name);
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