# IonsAnalysis
This is the repository taking care of collecting all the progress and all the stuff I make analyzing AMS-02 ion fluxes. This work is based on the GitLab PG-RM2-IonsAnalysis repo (https://gitlab.cern.ch/ams-italy/analysis/pg-rm2-ionsanalysis)

- - - In the getData directory there is all the stuff involved on taking data from pass8. The logical order is the following:
	  - the script find_and_sort.sh take all the n-tuples (in this case from v.1.0.0) and split them in 100 .txt files. The script saves them on the 'ntuples' dir along with the file 'data_filelist.txt' ready to be submitted to condor;
- - - The checkVar dir is involved in plotting all the variables (charge, beta , rigidity etc..) to justify the paramters used in the selection (cuts);
- - - In the data_test/ dir there are some test.
