# IonsAnalysis
This is the repository taking care of collecting all the progress and all the stuff I make analyzing AMS-02 ion fluxes. This work is based on the GitLab PG-RM2-IonsAnalysis repo (https://gitlab.cern.ch/ams-italy/analysis/pg-rm2-ionsanalysis)

To install NSL:
git clone --branch v1.0.1 ssh://git@gitlab.cern.ch:7999/ams-italy/nsl.git and then follow the guide on https://nsl-readthedocs.readthedocs.io/en/latest/build-install.html. The -DCMAKE_INSTALL_PREFIX=${your-install-path-here} need to be -DCMAKE_INSTALL_PREFIX=/storage/gpfs_ams/ams/users/aubaldi/nsl.install.

To set the configuration for cmake, inside build do:
cmake .. -DNAIA_DIR=/cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/naia/v1.0.2/cmake -DNSL_DIR=/storage/gpfs_ams/ams/users/aubaldi/nsl.install/cmake


- - - In the getData directory there is all the stuff involved on taking data from pass8. The logical order is the following:
	  - the script find_and_sort.sh take all the n-tuples (in this case from v.1.0.0) and split them in 100 .txt files. The script saves them on the 'ntuples' dir along with the file 'data_filelist.txt' ready to be submitted to condor;
- - - The Data directory contains the outputs from the executables: the logic is Data/<name_of_executable>/ and here there is:
	  - condor (with all the stuff and the submission code);
	  - actually ouptut;
- - - checkVar is involved in plotting all the variables (charge, beta , rigidity etc..) to justify the paramters used in the selection (cuts);
- - - In the data_test/ dir there are some test.
