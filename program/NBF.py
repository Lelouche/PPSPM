# NBF - Numerical Bloom filter program to convert numerical 
# (integer, floating point, and modulus) values into Bloom filters to 
# allow privacy-preserving similarity calculations.
#
# DV and PC, Mar 2015
# -----------------------------------------------------------------------------

# imports
#

# Standard Python modules
import math
import random
import hashlib

from bitarray import bitarray

# -----------------------------------------------------------------------------

class NBF():

  def __init__(self, bf_len, bf_num_hash_func, bf_num_inter, bf_step, 
               max_abs_diff, min_val, max_val):
    """Initialisation, set class parameters:
       - bf_len            Length of Bloom filters
       - bf_num_hash_func  Number of hash functions
       - bf_num_interval   Number of intervals to use for BF based similarities

       - max_abs_diff      Maximum absolute difference allowed
       - min_val           Minimum value 
       - max_val           Maximum value
    """

    self.bf_len =           bf_len
    self.bf_num_hash_func = bf_num_hash_func
    self.bf_num_inter =     bf_num_inter
    self.bf_step = 	    bf_step

    self.max_abs_diff =  max_abs_diff
    self.min_val =       min_val
    self.max_val =       max_val

    assert max_val > min_val

    # Bloom filter shortcuts
    #
    self.h1 = hashlib.sha1
    self.h2 = hashlib.md5

  # ---------------------------------------------------------------------------

  def __set_to_bloom_filter__(self, val_set):
    """Convert an input set of values into a Bloom filter.
    """

    k = self.bf_num_hash_func
    l = self.bf_len

    bloom_set = bitarray(l)  # Convert set into a bit array		
    bloom_set.setall(False)
         
    for val in val_set:
      hex_str1 = self.h1(val).hexdigest()
      int1 =     int(hex_str1, 16)
      hex_str2 = self.h2(val).hexdigest()
      int2 =     int(hex_str2, 16)

      for i in range(k):
        gi = int1 + i*int2
        gi = int(gi % l)
        bloom_set[gi] = True
      
    return bloom_set

  # ---------------------------------------------------------------------------

  def calc_abs_diff(self, val1, val2):
    """Calculate absolute difference similarity between two values based on the
       approach described in:

       Data Matching, P Christen, Springer 2012, page 121, equations (5.28).
    """

    max_abs_diff = self.max_abs_diff

    if (val1 == val2):
      return 1.0

    abs_val_diff = abs(float(val1) - float(val2))

    if (abs_val_diff >= max_abs_diff):
      return 0.0  # Outside allowed maximum difference

    abs_sim = 1.0 - abs_val_diff / max_abs_diff

    assert abs_sim > 0.0 and abs_sim < 1.0, (val1, val2, abs_sim)

    return abs_sim

  # ---------------------------------------------------------------------------

  def nbf_calc_abs_diff(self, val1, val2):
    """Calculate absolute difference similarity between two 
       numerical values encoded into Bloom filters.
    """

    # Number of intervals and their sizes (step) to consider
    #
    bf_num_inter = self.bf_num_inter
    bf_step =           self.bf_step

    val1_set = set()
    val2_set = set()

    rem_val1 = val1 % bf_step  # Convert into values within same interval
    if rem_val1 >= bf_step/2:
      use_val1 = val1 + (bf_step - rem_val1)
    else:
      use_val1 = val1 - rem_val1

    rem_val2 = val2 % bf_step
    if rem_val2 >= bf_step/2:
      use_val2 = val2 + (bf_step - rem_val2)
    else:
      use_val2 = val2 - rem_val2

    val1_set.add(str(float(use_val1)))  # Add the actual value
    val2_set.add(str(float(use_val2)))  # Add the actual value

    # Add variations larger and smaller than the actual value
    #
    for i in range(bf_num_inter+1):
      diff_val = (i+1)*bf_step
      val1_set.add(str(use_val1 - diff_val))
      val2_set.add(str(use_val2 - diff_val))

      diff_val = (i)*bf_step
      val1_set.add(str(use_val1 + diff_val))
      val2_set.add(str(use_val2 + diff_val))

    return val1_set, val2_set

  # ---------------------------------------------------------------------------

  def calc_mod_diff(self, val1, val2):
    """Calculate difference similarity between two modulus values.
    """

    max_abs_diff = self.max_abs_diff
    min_val = self.min_val
    max_val = self.max_val

    if (val1 == val2):
      return 1.0

    mod_val_diff = float((max_val - max(val1,val2)) + (min(val1,val2)-min_val)+1)
    print mod_val_diff
    if (mod_val_diff >= max_abs_diff):
      return 0.0  # Outside allowed maximum difference

    mod_sim = 1.0 - mod_val_diff / max_abs_diff

    assert mod_sim > 0.0 and mod_sim < 1.0, (val1, val2, mod_sim)

    return mod_sim

  # ---------------------------------------------------------------------------

  def nbf_calc_mod_diff(self, val1, val2):
    """Calculate absolute modulus difference similarity between two 
       numerical modulus values encoded into Bloom filters.
    """

    # Number of intervals and their sizes (step) to consider
    #
    bf_num_inter = self.bf_num_inter
    bf_step =           self.bf_step

    min_val = self.min_val
    max_val = self.max_val

    val1_set = set()
    val2_set = set()

    rem_val1 = val1 % bf_step  # Convert into values within same interval
    if rem_val1 >= bf_step/2:
      use_val1 = val1 + (bf_step - rem_val1)
    else:
      use_val1 = val1 - rem_val1

    rem_val2 = val2 % bf_step
    if rem_val2 >= bf_step/2:
      use_val2 = val2 + (bf_step - rem_val2)
    else:
      use_val2 = val2 - rem_val2

    val1_set.add(str(float(use_val1)))  # Add the actual values
    val2_set.add(str(float(use_val2)))  # Add the actual values

    # Add variations larger and smaller than the actual value
    #
    for i in range(bf_num_inter+1):

      diff_val = (i+1)*bf_step
      prev_val1 = use_val1 - diff_val
      if prev_val1 < min_val:
        val1_set.add(str(prev_val1 + (max_val-min_val+1)))
      else:
        val1_set.add(str(prev_val1))
      prev_val2 = use_val2 - diff_val
      if prev_val2 < min_val:
        val2_set.add(str(prev_val2 + (max_val-min_val+1)))
      else:
        val2_set.add(str(prev_val2))

      diff_val = (i)*bf_step
      next_val1 = use_val1 + diff_val
      if next_val1 > max_val:
        val1_set.add(str(next_val1%(max_val-min_val)))
      else:
        val1_set.add(str(next_val1))
      next_val2 = use_val2 + diff_val
      if next_val2 > max_val:
        val2_set.add(str(next_val2%(max_val-min_val)))
      else:
        val2_set.add(str(next_val2))

    return val1_set, val2_set

  # ---------------------------------------------------------------------------

  def gen_rand_val_pairs(self):
    """Generating random pairs of values for testing.
    """

    max_abs_diff = self.max_abs_diff
    min_val =      self.min_val
    max_val =      self.max_val

    rand_val_pairs = []

    for i in range(10000):
      v1 = random.randint(min_val, max_val)

      # Make second value to be close (within 2 * maximum difference)
      #
      v2_min = max(min_val, v1-int(2*max_abs_diff*(max_val-min_val)))
      v2_max = min(max_val, v1+int(2*max_abs_diff*(max_val-min_val)))
      v2 = random.randint(v2_min, v2_max)

      assert v1 >= 1
      assert v2 >= 1

      rand_val_pairs.append([v1,v2])

    #print rand_val_pairs
    return rand_val_pairs

