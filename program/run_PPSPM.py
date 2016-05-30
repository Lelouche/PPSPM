# Python program script for privacy-presereving similar patient matching  
# (PP-SPM) using Bloom filters
#
# This is the program controls the overall experiments and starts the programs
#
# DV and PC, Jun 2015
# ----------------------------------------------------------------------------

# Standard Python modules
import os
import subprocess  
import sys
import time

# Generate a time string to be used for the log file
#
curr_time_tuple = time.localtime()
curr_time_str = str(curr_time_tuple[0]) + str(curr_time_tuple[1]).zfill(2) + \
                str(curr_time_tuple[2]).zfill(2) + '-' + \
                str(curr_time_tuple[3]).zfill(2) + \
                str(curr_time_tuple[4]).zfill(2)

# ----------------------------------------------------------------------------
# Start the program 
#

datasets = ['Cleveland','Breast cancer diagnostic data','Pima Indians Diabetes data','MQIC']

block_attrs = {'Cleveland':'236', 'Breast cancer diagnostic data':'2', 'Pima Indians Diabetes data':'9', 'MQIC':'134567'}
match_attrs = {'Cleveland':'145', 'Breast cancer diagnostic data':'456', 'Pima Indians Diabetes data':'268', 'MQIC':'28'}

m_list = ['1','5','10'] 

encode = ['abf','clk','rbf']

for dataset in datasets:
  for m in m_list:
    for enc in encode:
        
        print
        print '----------------------------------------------------------------'
        print '  Starting experiment on:', time.asctime(time.localtime())
        print '  Data sets used:', dataset
        print '  Blocking attributes used: %s' % block_attrs[dataset]
        print '  Matching attributes used: %s' % match_attrs[dataset]
        print '  m: %s' % m
        print '  Encode: %s' % enc
        print

        print 'Starting processes.'
        
        #database = '../datasets/database_recs_%s.csv' % dataset
        database = '..'+os.sep+'datasets'+os.sep+'database_recs_%s.csv' % dataset
        #query = '../datasets/query_recs_%s.csv' % dataset
        query = '..'+os.sep+'datasets'+os.sep+'query_recs_%s.csv' % dataset

        cmd_line_list = ['python', 'PPSPM.py', database, query, block_attrs[dataset], match_attrs[dataset], m, enc]

        party_proc = subprocess.Popen(cmd_line_list) 


        print '  Waiting for processes to complete...'
        print

        ret_code = party_proc.wait()
        print 'finished.'

        if (ret_code != 0):
          print 'returned code:', ret_code
          sys.exit()  # Stop experiment
