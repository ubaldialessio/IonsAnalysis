#!/bin/bash

find /storage/gpfs_ams/ams/groups/AMS-Italy/ntuples/v1.0.0/ISS.B1236/pass8/ -name *.root | sort > ./ntuples/data_list.txt
cd ntuples/
split -l 100 -d -a 5 --additional-suffix .txt data_list.txt file_
ls -d $PWD/* > data_filelist.txt
