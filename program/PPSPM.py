# Python program for similar patient matching (SPM) 
# (Privacy-preserving and non privacy-preserving).
#
# DV and PC, Jun 2015
# ----------------------------------------------------------------------------

# imports
#

# Febrl modules
import auxiliary
import stringcmp

# Numerical Bloom filter module
from NBF import NBF

# Standard Python modules
import random
import os
import math
import sys
import time
import gzip
import string
import gc
import shutil
import decimal
import operator
import hashlib

import numpy
from bitarray import bitarray

#--------------------------------------------------------------------------------

class SPM:

  def __init__(self, m, length, use_attr_index, blk_attr_index, encode):  
    """Constructor. Initialise an index, set input parameters.

       Arguments:
       - m 		       Number of similar records	
       - length		       Number of bits used in the Bloom filter
       - use_attr_index        A list of the index numbers of the attributes
                               that are used in the linking.
       - blk_attr_index        A list of the index numbers of the attributes
                               that are used in the blocking.
       - encode 	       Bloom filter encoding method (field, clk, rbf, clkrbf)
   """

    self.m	 	        =     	 m
    self.length 	        =     	 length
    self.use_attr_index         =   	 use_attr_index
    self.blk_attr_index         =  	 blk_attr_index
    self.encode 	        = 	 encode

    self.h1                     = hashlib.sha1	# hash functions used for Bloom filter masking
    self.h2                     = hashlib.md5	# hash functions used for Bloom filter masking

    # Define the main data structures to be used

    self.rec_dict 		= {}  		# To hold the data set loaded from a file. Keys in
                        			# this dictionary are the record identifiers and
                        			# values are the cleaned records (without the
                        			# identifiers)

    self.query_dict		= {}		# To hold the query records loaded from a file. Keys in
                        			# this dictionary are the record identifiers and
                        			# values are the cleaned records (without the
                        			# identifiers)

    self.attr_info              = {}            # For storing attribute parameter values

    self.block_index 		= {}		# The block index data structure. Keys are the blocking key 
                                                # values and values are list of attribute values

    self.candidate_dict		= {}		# To store the true candidate matches for each query record along
						# with their similarities (SPM)

    self.results_dict		= {}		# To store the true top m matches for each query record (SPM)
  
    self.mcandidate_dict	= {}    	# masked candidate dict - for storing candidate matches
 						# calculated in a privacy-preserving context for each query record (PP-SPM)

    self.mresults_dict		= {}		# masked results dict - for storing the top m matches
 						# calculated in a privacy-preserving context for each query record (PP-SPM)

    self.sim_dict 		= {}		# For storing the similarity values of matches required for rank calculation
    self.rank_dict 		= {}		# For ranking results
    self.mrank_dict 		= {}		# For ranking masked results

  # ---------------------------------------------------------------------------

  def read_database(self, file_name):
    """Load the database and store in memory as a dictionary with record
       identifiers as keys.
    """

    rec_dict = self.rec_dict  

    print 'Load data file: %s' % (file_name)

    if (file_name.lower().endswith('.gz')):
      in_file =  gzip.open(file_name)  # Open gzipped file
    else:
      in_file =  open(file_name)       # Open normal file

    for rec in in_file:
      rec = rec.lower().strip()
      rec = rec.split(',')
      clean_rec = map(string.strip, rec) # Remove all surrounding whitespaces        
        
      rec_id = clean_rec[0]
      assert rec_id not in rec_dict, ('Record ID not unique:', rec_id)

      # Save the original record without the record identifier
      #
      rec_dict[rec_id] = clean_rec[1:]

    print '    Read %d records' % (len(rec_dict))
    print 

  # ---------------------------------------------------------------------------

  def analyse_database(self,b):
    """Analyse matching attribute values in the database to set the parameters.
    """
    
    rec_dict        = self.rec_dict  
    use_attr_index  = self.use_attr_index
    attr_info       = self.attr_info
    encode          = self.encode
    length 	    = self.length

    if encode == 'abf':		  # Encoding is ABF
      for attr in use_attr_index: # loop over matching attributes
        this_attr_vals = []
        attr_type = ''
        for rec in rec_dict.values():
          this_attr_val = rec[attr]
          if self.isfloat(this_attr_val) or self.isint(this_attr_val):
            attr_type = 'num'
            this_attr_val = float(this_attr_val)
            this_attr_vals.append(this_attr_val)

        # Numerical attribute
        #
        if attr_type == 'num':
          nei_list = []
          d = 0.05*(max(this_attr_vals) - min(this_attr_vals))
          # calculate n - len(set(nei_list))
          uniq_elem_list = set(this_attr_vals)
          for uniq_elem in uniq_elem_list:
            nei_list.append(uniq_elem)
            rem_uniq_elem_val = uniq_elem % (d/b)  # Convert into values within same interval
            if rem_uniq_elem_val >= 0.5*(d/b):
              nei_list.append(uniq_elem + ((d/b) - rem_uniq_elem_val))
            else:
              nei_list.append(uniq_elem - rem_uniq_elem_val)
          # calculate k - num_hash
          num_hash = math.ceil(float(length/len(set(nei_list))) * numpy.log(2))
          #num_hash = 10
          attr_info[attr] = ['num',min(this_attr_vals), max(this_attr_vals), d, int(num_hash), length] 

        # Categorical attribute
        #
        elif len(set(this_attr_vals)) <= 0.05 * len(this_attr_vals):  # assumed to have at most 0.05% values of all values
          # calculate k - num_hash based on n = len(set(this_attr_vals))
          num_hash = math.ceil(float(length/len(set(this_attr_vals))) * numpy.log(2))
          #num_hash = 10
          attr_info[attr] = ['cate', int(num_hash), length]

        # String attribute
        #  
        else:
          # calculate k - num_hash based on n = len(set(this_attr_vals))
          num_hash = math.ceil(float(length/len(set(this_attr_vals))) * numpy.log(2))
          #num_hash = 10
          attr_info[attr] = ['string', int(num_hash), length]

    elif encode == 'clk':	  # Encoding is CLK
      nei_list = []
      len_vals = 0
      for attr in use_attr_index: # loop over matching attributes
        this_attr_vals = []
        attr_type = ''
        for rec in rec_dict.values():
          this_attr_val = rec[attr]
          if self.isfloat(this_attr_val) or self.isint(this_attr_val):
            attr_type = 'num'
            this_attr_val = float(this_attr_val)
            this_attr_vals.append(this_attr_val)

        if attr_type == 'num':
          d = 0.05*(max(this_attr_vals) - min(this_attr_vals))

          uniq_elem_list = set(this_attr_vals)
          for uniq_elem in uniq_elem_list:
            nei_list.append(uniq_elem)
            rem_uniq_elem_val = uniq_elem % (d/b)  
            if rem_uniq_elem_val >= 0.5*(d/b):
              nei_list.append(uniq_elem + ((d/b) - rem_uniq_elem_val))
            else:
              nei_list.append(uniq_elem - rem_uniq_elem_val)
          attr_info[attr] = ['num',min(this_attr_vals), max(this_attr_vals), d] 

        elif len(set(this_attr_vals)) <= 0.05 * len(this_attr_vals):  
          nei_list += this_attr_vals
          len_vals += len(this_attr_vals)
          attr_info[attr] = ['cate']
        else:
          nei_list += this_attr_vals
          len_vals += len(this_attr_vals)
          attr_info[attr] = ['string']

      # calculate k - num_hash based on record-level n = len(set(nei_list))
      num_hash = math.ceil(float(length/len(set(nei_list))) * numpy.log(2))
      #num_hash = 10
      for atr in attr_info:
        attr_info[atr] += [int(num_hash), length]

    elif encode == 'rbf':	  # Encoding is RBF
      lengths = []
      nei_list = []
      len_vals = 0
      for attr in use_attr_index: # loop over matching attributes
        this_attr_vals = []
        attr_type = ''
        for rec in rec_dict.values():
          this_attr_val = rec[attr]
          if self.isfloat(this_attr_val) or self.isint(this_attr_val):
            attr_type = 'num'
            this_attr_val = float(this_attr_val)
            this_attr_vals.append(this_attr_val)

        if attr_type == 'num':
          d = 0.05*(max(this_attr_vals) - min(this_attr_vals))

          uniq_elem_list = set(this_attr_vals)
          for uniq_elem in uniq_elem_list:
            nei_list.append(uniq_elem)
            rem_uniq_elem_val = uniq_elem % (d/b)  
            if rem_uniq_elem_val >= 0.5*(d/b):
              nei_list.append(uniq_elem + ((d/b) - rem_uniq_elem_val))
            else:
              nei_list.append(uniq_elem - rem_uniq_elem_val)
          attr_info[attr] = ['num',min(this_attr_vals), max(this_attr_vals), d] 

        elif len(set(this_attr_vals)) <= 0.05 * len(this_attr_vals):  
          nei_list += this_attr_vals
          len_vals += len(this_attr_vals)
          attr_info[attr] = ['cate']
        else:
          nei_list += this_attr_vals
          len_vals += len(this_attr_vals)
          attr_info[attr] = ['string']

        lengths.append(len(set(this_attr_vals)))

      # calculate k - num_hash based on record-level n = len(set(nei_list))
      num_hash = math.ceil(float(length/len(set(nei_list))) * numpy.log(2))
      #num_hash = 10
      for atr in attr_info:
        attr_info[atr] += [int(num_hash)]

      # calculate length based on weight of attributes
      for qid in range(len(lengths)):
        #ratio = (float(lengths[qid])/sum(lengths)) * length
        ratio = (float(length)/min(lengths)) * lengths[qid]
        attr_info[use_attr_index[qid]] += [int(math.ceil(ratio))]
        
    #print attr_info

  # ---------------------------------------------------------------------------

  def isfloat(self,x):
    """A method that returns true if an input value is a floating point data
    """

    try:
        a = float(x)
    except ValueError:
        return False
    else:
        return True

  # ---------------------------------------------------------------------------

  def isint(self,x):
    """A method that returns true if an input value is a integer data
    """

    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b

  # ---------------------------------------------------------------------------

  def stddev(self,x):
    """Compute the standard deviation of a list of numbers.
    """

    assert len(x) > 0

    if (len(x) == 1):  # Only one element in list
      return 0.0
    else:
      sum = 0.0
      for i in x:
        sum += i
      avrg = sum / float(len(x))
      sum = 0.0
      for i in x:
        sum = sum + (i - avrg) * (i - avrg)

      return (math.sqrt(sum / float(len(x))))

  # ---------------------------------------------------------------------------

  def build_BI(self):
    """Build block_index data structure to store the BKVs and
       the corresponding list of record identifiers in the database.
    """

    rec_dict        = self.rec_dict  
    block_index     = self.block_index
    blk_attr_index  = self.blk_attr_index

    print 'Build Block Index for attributes:', blk_attr_index

    for (rec_id, clean_rec) in rec_dict.iteritems():
      compound_bkv = ""

      if blk_attr_index == []:        # No blocking keys
        compound_bkv = 'No_blk'
      else:
        for attr in blk_attr_index:   # Process selected blocking attributes
          attr_val = clean_rec[attr]
          attr_encode = attr_val      # Actual categorical value as BKV
          compound_bkv += attr_encode
      
      if (compound_bkv in block_index): 	 # Block value in index, only add attribute value
        rec_id_list = block_index[compound_bkv]
        rec_id_list.append(rec_id)
      else:					 # A new block, add block value and attribute value
        rec_id_list = [rec_id]
        block_index[compound_bkv] = rec_id_list

    print '    Generated %d blocks' % (len(block_index))
    print

  # ---------------------------------------------------------------------------

  def read_query(self, file_name):
    """Load the query set and store in memory as a dictionary with record
       identifiers as keys.
    """

    query_dict = self.query_dict  

    print 'Load query file: %s' % (file_name)

    if (file_name.lower().endswith('.gz')):
      in_file =  gzip.open(file_name)  # Open gzipped file
    else:
      in_file =  open(file_name)  # Open normal file

    for rec in in_file:
      rec = rec.lower().strip()
      rec = rec.split(',')
      clean_rec = map(string.strip, rec) # Remove all surrounding whitespaces        
        
      rec_id = clean_rec[0]
      assert rec_id not in query_dict, ('Query ID not unique:', rec_id)

      # Save the original record without the record identifier
      #
      query_dict[rec_id] = clean_rec[1:]

    print '    Read %d queries' % (len(query_dict))
    print 

  # ---------------------------------------------------------------------------

  def compare(self,b):
    """Compare query rec with database records for SPM.
    """
    
    rec_dict        = self.rec_dict
    query_dict      = self.query_dict  
    block_index     = self.block_index
    candidate_dict  = self.candidate_dict
    blk_attr_index  = self.blk_attr_index
    use_attr_index  = self.use_attr_index
    attr_info       = self.attr_info
    m	 	    = self.m
    length 	    = self.length


    for (rec_id, clean_rec) in query_dict.iteritems(): # Iterate over query records
      bkv = ''
      if blk_attr_index == []:
        bkv = 'No_blk'
      else:
        for attr in blk_attr_index:
          attr_val = clean_rec[attr]
          bkv += attr_val

      if bkv in block_index:
        cand_recs = block_index[bkv]  # candidate records for the query
        for cand_rec in cand_recs:    # Iterate over candidate records
          tot_sim_val = 0.0
          for qid in use_attr_index:  # Loop over attributes
            query_val = clean_rec[qid]
            database_val = rec_dict[cand_rec][qid]
            #qid_type = attr_type[qid]
            qid_type = attr_info[qid][0]

            if qid_type == 'num':     # Numerical attribute
              num_attr_info = attr_info[qid]
	      nbf = NBF(length, num_attr_info[-1], b, num_attr_info[3]/(2*b), num_attr_info[3], num_attr_info[1], num_attr_info[2]) 
              abs_sim = nbf.calc_abs_diff(float(database_val), float(query_val))

            elif qid_type == 'cate':  # Categorical attribute
              if query_val == database_val:
                abs_sim = 1.0
              else:
                abs_sim = 0.0
           
            else:  		      # String attribute
              abs_sim = stringcmp.editdist(query_val,database_val)
 
            tot_sim_val += abs_sim	     # Sum similarity over attributes

          tot_sim_val /= len(use_attr_index) # Calculate overall similarity
           
          # Store similarity results in candidate_dict
          if rec_id in candidate_dict:
            this_rec_dict = candidate_dict[rec_id] 
            this_rec_dict[cand_rec] = tot_sim_val
          else:
            this_rec_dict = {}
            this_rec_dict[cand_rec] = tot_sim_val
            candidate_dict[rec_id] = this_rec_dict
             
      else:
	candidate_dict[rec_id] = {} # no similar records for the query
	
      if rec_id not in candidate_dict:
        candidate_dict[rec_id] = {}
 
    #print 'actual candidates:', candidate_dict

  # ---------------------------------------------------------------------------

  def compare_masked(self,b):
    """Compare query rec with database records in a privacy-preserving manner (PP-SPM).
    """
    
    rec_dict        = self.rec_dict
    query_dict      = self.query_dict  
    block_index     = self.block_index
    mcandidate_dict = self.mcandidate_dict
    blk_attr_index  = self.blk_attr_index
    use_attr_index  = self.use_attr_index
    encode          = self.encode
    attr_info       = self.attr_info
    m	 	    = self.m
    length 	    = self.length
    
    for (rec_id, clean_rec) in query_dict.iteritems(): # Iterate over query records
      bkv = ''
      if blk_attr_index == []:
        bkv = 'No_blk'
      else:
        for attr in blk_attr_index:
          attr_val = clean_rec[attr]
          bkv += attr_val

      if bkv in block_index:
        cand_recs = block_index[bkv]  # candidate records for the query
        for cand_rec in cand_recs:    # Iterate over candidate records
          tot_sim_val = 0.0
          rbf1 = bitarray(length)     # record-level Bloom filters 		
          rbf1.setall(False)
          rbf2 = bitarray(length)  		
          rbf2.setall(False)

          for qid in use_attr_index:  # Loop over attributes
            query_val = clean_rec[qid]
            database_val = rec_dict[cand_rec][qid]
            qid_type = attr_info[qid][0]

            if qid_type == 'num':    # Numerical attribute
              num_attr_info = attr_info[qid]
	      nbf = NBF(num_attr_info[-1], num_attr_info[-2], b, num_attr_info[3]/(2 *b), num_attr_info[3], num_attr_info[1], num_attr_info[2]) 
              val1_set, val2_set  = nbf.nbf_calc_abs_diff(float(database_val), float(query_val))
            elif qid_type == 'cate': # Categorical attribute
              if query_val == database_val:
	        nbf = NBF(num_attr_info[-1], num_attr_info[-2], 0, 0, 0, 0, 0) 
                val1_set = [database_val]
                val2_set = [query_val]
            else:  		     # String attribute
	      nbf = NBF(num_attr_info[-1], num_attr_info[-2], 0, 0, 0, 0, 0) 
              val1_set = [query_val[i:i+q] for i in range(len(query_val) - (q-1))]
              val2_set = [database_val[i:i+q] for i in range(len(database_val) - (q-1))]

            bf1 = nbf.__set_to_bloom_filter__(val1_set)  	# Convert lists to Bloom filters
            bf2 = nbf.__set_to_bloom_filter__(val2_set)
            comm_bits = (bf1 & bf2).count()		 	# Calculate common 1-bits
            nbf_abs_sim = 2 * float(comm_bits) / (int(bf1.count())+int(bf2.count()))  # Calculate ABF similarity
            assert nbf_abs_sim >= 0.0 and nbf_abs_sim <= 1.0, (val1, val2, nbf_abs_sim)

            if encode == 'abf':
              tot_sim_val += nbf_abs_sim	# Sum ABF similarities
            elif encode == 'clk':
              rbf1 |= bf1			# Union over ABF to generate CLK
              rbf2 |= bf2
            elif encode == 'rbf':		# Concatenate over ABF to generate RBF
              rbf1 += bf1
              rbf2 += bf2

          if encode == 'abf':				# Calculate overall ABF similarity
            tot_sim_val /= len(use_attr_index)

          elif encode == 'clk' or encode == 'rbf': 	# Calculate overall record-level Bloom filter similarity
            comm_bits = (rbf1 & rbf2).count()
            tot_sim_val = 2 * float(comm_bits) / (int(rbf1.count())+int(rbf2.count()))

	  # Store similarity results in mcandidate_dict
          if rec_id in mcandidate_dict:
            this_rec_dict = mcandidate_dict[rec_id] 
            this_rec_dict[cand_rec] = tot_sim_val
          else:
            this_rec_dict = {}
            this_rec_dict[cand_rec] = tot_sim_val
            mcandidate_dict[rec_id] = this_rec_dict
               
      else:
	mcandidate_dict[rec_id] = {} # no similar records for the query

      if rec_id not in mcandidate_dict:
        mcandidate_dict[rec_id] = {}
 
    #print 'masked candidates:', mcandidate_dict

  # ---------------------------------------------------------------------------

  def find_m_similar(self):
    """Find m most similar records for each query based on the overall similarity.
    """

    candidate_dict  = self.candidate_dict
    m	 	    = self.m
    results_dict    = self.results_dict
    sim_dict        = self.sim_dict

    for query in candidate_dict:	# loop over query records
      this_query_res = []
      this_sim_res = []
      this_query_dict = candidate_dict[query]
      if this_query_dict != {}:
        # Sort by similarity and retrieve m records with highest similarity
        #
        sorted_this_query_dict = sorted(this_query_dict.items(), key=operator.itemgetter(1))
        if len(sorted_this_query_dict) >= m:
          for x in range(m):
            this_query_res.append(sorted_this_query_dict[-(x+1)][0])
            this_sim_res.append(sorted_this_query_dict[-(x+1)][1])
        else:
          for x in range(len(sorted_this_query_dict)):
            this_query_res.append(sorted_this_query_dict[-(x+1)][0])
            this_sim_res.append(sorted_this_query_dict[-(x+1)][1])

        results_dict[query] = this_query_res
        sim_dict[query] = this_sim_res 

      else:			        # No similar records
        results_dict[query] = this_query_res
        sim_dict[query] = this_sim_res

    #print 'matches:', results_dict
    #print 'matches similarities:', sim_dict

  # ---------------------------------------------------------------------------

  def find_m_similar_masked(self):
    """Find m most similar records for each query based on the overall similarity
       in a privacy-preserving setting.
    """

    mcandidate_dict = self.mcandidate_dict
    m	 	    = self.m
    mresults_dict   = self.mresults_dict
    mrank_dict      = self.mrank_dict

    results_dict    = self.results_dict
    sim_dict        = self.sim_dict
    rank_dict       = self.rank_dict

    for query in mcandidate_dict:	# loop over query records
      this_query_res = []
      this_sim_res = []
      this_query_dict = mcandidate_dict[query]
      if this_query_dict != {}:
        # Sort by similarity and retrieve m records with highest similarity
        #
        sorted_this_query_dict = sorted(this_query_dict.items(), key=operator.itemgetter(1))
        if len(sorted_this_query_dict) >= m:
          for x in range(m):
            this_query_res.append(sorted_this_query_dict[-(x+1)][0])
            this_sim_res.append(sorted_this_query_dict[-(x+1)][1])
        else:
          for x in range(len(sorted_this_query_dict)):
            this_query_res.append(sorted_this_query_dict[-(x+1)][0])
            this_sim_res.append(sorted_this_query_dict[-(x+1)][1])

        # Rank similarities
        actual_rank_res = []
        masked_rank_res = []
        intersect_items = list(set(results_dict[query]).intersection(this_query_res))
        for inter in intersect_items:
          pos = results_dict[query].index(inter)
          actual_rank_res.append(sim_dict[query][pos])
          pos = this_query_res.index(inter)
          masked_rank_res.append(this_sim_res[pos])

        mresults_dict[query] = this_query_res
        rank_dict[query] = list(self.argsortdup(actual_rank_res))
        mrank_dict[query] = list(self.argsortdup(masked_rank_res))

      else:			     # No similar records
        mresults_dict[query] = this_query_res
        mrank_dict[query] = []
        rank_dict[query] = []

    #print 'masked matches:', mresults_dict
    #print 'actual matches ranking:', rank_dict
    #print 'masked matches ranking:', mrank_dict

  # ---------------------------------------------------------------------------

  def argsortdup(self,a):
    """A function to rank a list of values.
    """
    
    sorted = numpy.sort(a)
    ranked = []
    for item in a:
      ranked.append(sorted.searchsorted(item)+1)
    return numpy.array(ranked)

  # ---------------------------------------------------------------------------

  def calculate_accuracy(self):
    """Calculate accuracy of privacy-preserving comparison
       using actual values comparison results as the truth data.
    """

    results_dict  = self.results_dict
    mresults_dict = self.mresults_dict
    rank_dict     = self.rank_dict
    mrank_dict    = self.mrank_dict
    accuracy_dict = {}

    for query in mresults_dict:
      query_res = mresults_dict[query]
      actual_res = results_dict[query]

      tot_m = len(actual_res)
      m = len(query_res)
      tm = 0
      for res in query_res:
        if res in actual_res:
          tm += 1

      # Calculate precision, recall, and F1 measures
      if m >0:
        prec = tm / float(m)
      else:
        prec = 0
      if tot_m > 0.0:
        rec = tm / float(tot_m)
      else:
        rec = 0.0
      if (prec+rec) > 0.0:
        fsco = (2 * prec * rec) / float(prec + rec)
      else:
        fsco = 0.0

      # Calculate Spearman's rank correlation
      query_rank = mrank_dict[query]
      actual_rank = rank_dict[query]
      assert len(actual_rank) == len(query_rank)
      n = len(actual_rank)

      if n==1:
        if actual_rank[0] == query_rank[0]:
          rank_cor = 1.0
        else:
          rank_cor = 0.0
      elif n==0:
        rank_cor = 0.0
      else:
        dist_sqr = 0
        for x in range(n):
          dist = actual_rank[x] - query_rank[x]
          dist_sqr += dist**2
        dist_sqr *= 6
        tot_elem = float(n*(n**2 -1))
        rank_cor = 1.0 - (dist_sqr/tot_elem)

      accuracy_dict[query] = [prec,rec,fsco, rank_cor]

    print 'accuracy_dict:', accuracy_dict
    return accuracy_dict    

  # ---------------------------------------------------------------------------

# Input parameter setting for the framework
#
database_file =     sys.argv[1]
query_file =     sys.argv[2]
block_atr_str = sys.argv[3]
attr_str =       sys.argv[4]
m = int(sys.argv[5]) 
encode = sys.argv[6]

length = 2000
b = 5 #10
#num_hash_func = 10
#d_max - 5% of (max-min)

# Check and convert the linkage attribute list
#
attr_list_str = list(attr_str)
attr_list = []
for attr in attr_list_str:
  attr_list.append(int(attr)-1)

# Check and convert the block attribute list
#
blk_attr_list_str = list(block_atr_str)
blk_attr_list = []
for attr in blk_attr_list_str:
  blk_attr_list.append(int(attr)-1)

# Initialize the framework
#
spm = SPM(m, length, attr_list, blk_attr_list, encode) 

# Log file to write the results
#
log_file_name = '../results/'+database_file.split('_')[2]+'_'+encode+'_'+str(m)
log_file = open(log_file_name, 'w')
log_file.write(block_atr_str+','+attr_str+','+str(m)+ ','+encode + os.linesep)

# Blocking
#
start_time = time.time()
spm.read_database(database_file)
spm.analyse_database(b)
spm.build_BI()
spm.read_query(query_file)
blocking_phase_time = time.time() - start_time 

# Matching
#
start_time = time.time()
spm.compare(b)
spm.find_m_similar()
matching_phase_time = time.time() - start_time

# Privacy-preserving matching
#
start_time = time.time()
spm.compare_masked(b)
spm.find_m_similar_masked()
masked_matching_phase_time = time.time() -start_time

# Calculate total runtime for PP-SPM
tot_time = blocking_phase_time + masked_matching_phase_time
str_tot_time = '%.4f' %(tot_time)

# Calculate total memory usage for PP-SPM
memo_usage =     auxiliary.get_memory_usage() 
memo_usage_val = auxiliary.get_memory_usage_val()
memo_usage_val = memo_usage_val if memo_usage_val else 0.0
str_mem_usage = '%.4f' %(memo_usage_val)

# write efficiency results into the log file
log_file.write(str_tot_time +','+ str_mem_usage + os.linesep)

# write effectiveness results (precision, recall, F1, and rank) 
# into the log file (one line per query record)
accuracy_dict = spm.calculate_accuracy()
for query in accuracy_dict:
  res_list = accuracy_dict[query]
  log_file.write(query)
  for res in res_list:
    log_file.write(','+str(res))
  log_file.write(os.linesep)

log_file.close()
