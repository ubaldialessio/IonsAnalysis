# IonsAnalysis
This is the repository taking care of collecting all the progress and all the stuff I make analyzing AMS-02 ion fluxes. This work is based on the GitLab PG-RM2-IonsAnalysis repo (https://gitlab.cern.ch/ams-italy/analysis/pg-rm2-ionsanalysis)

To install NSL:
git clone --branch v1.0.1 ssh://git@gitlab.cern.ch:7999/ams-italy/nsl.git and then follow the guide on https://nsl-readthedocs.readthedocs.io/en/latest/build-install.html. The -DCMAKE_INSTALL_PREFIX=${your-install-path-here} need to be -DCMAKE_INSTALL_PREFIX=/storage/gpfs_ams/ams/users/aubaldi/nsl.install.

To set the configuration for cmake, inside build do:
cmake .. -DNAIA_DIR=/cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/naia/v1.0.2/cmake -DNSL_DIR=/storage/gpfs_ams/ams/users/aubaldi/nsl.install/cmake

## selectEv
Is the executable involved in the selection of the given nuclei (charge). It computes the selection, the efficiencies, the mc acceptance, the livetime and the counts. Also produces the hooks histograms with the given selection.

## eff
Evaluate the total acceptance (data/mc corrections) and compute the actual flux.

## hplot
Usage: 
./hplot <charge> <single/same>
Produce a pdf file with the four efficiencies for a given nuclei (charge). With the option "single" it produces 1 efficiency per pages, with "both" overlaps (SAME) data and mc efficiencies.

## fragmentation
Usage: 
./fragmentation <charge> <path/to/input.root> <output.root> <AboveL1/BelowL1> 
Is involved in the production of the purity for a given nuclei. Actually only BelowL1 is implemented.

## CreateTemplatesBelowL1
Usage: ./CreateTemplatesBelowL1 <charge>
Once fragmentation has run over all the ntuples, this executable produce 3 files: one .root and two .pdf. The .root file contains two histogram that are the L1 charge distribution (to be fitted) and the L2 templates. The two pdf files are:
a) a pdf file containg the distribution of L1 and L2 for every rigity bin (22 in total), one bin per page);
b) a pdf file containg the summed of the distribution in every rigidity bin for both L1 and L2 drawn overlapped.

## CompareEff
Usage: ./CompareEff <charge1> <charge2> ... <charge n>
Create a pdf file comparing the efficiencies for the given nuclei (charges). For example by submitting ./CompareEff 5 6 8 it will produce some nice plots with the efficiencies for both data and mc

## UnfoldAcceptance
Usage: 
./UnfoldAcceptance <charge> <all/single> <reweigh/no> <rebin/no> <rigmap/no> <mc single/global> <normal/emulate> <period> -sec_track <y/n>
The logic for unfolding P acceptance is:
    1) Run UnfoldAcceptance on Si and P
    2) Run InterpolateCorrections
    3) Run UnfoldAcceptance on P
