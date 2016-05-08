Proof-of-concept software programs written in the 
Python programming language (version 2.7.3) for the
privacy-preserving similar patient matching (PP-SPM) framework.

Authors - Dinusha Vatsalan and Peter Christen,
The Australian National University,
23/06/2015

Python command to run experiments:
==================================================================
python run_PPSPM.py


The folder structure:
==================================================================
The PPSPM software contains three folders (directories), which are
1. datasets
2. program
3. results

1. datasets:
Four publicly available real-world medical datasets used by the authors 
for the experimental study of the PP-SPM framework are included - 
three from UCI repository 
http://archive.ics.uci.edu/ml/datasets.html: Heart Disease (Cleveland),
Breast Cancer, and Diabetes,
and one from Medical Quality Improvement Consortium 
http://visualizing.org/

These four datasets are split into two each, containing database and 
query records (randomly sampled 20%), respectively.

2. programs:
The script run_PPSPM.py includes experiment runs of PP-SPM framework
on the four datasets with different Bloom filter encoding
methods of ABF, CLK, and RBF and different values of the number
of similar records m [1,5,10].

This script calls the main program PPSPM.py which consists of several
functions required for SPM and PP-SPM. PP-SPM for numerical data is 
conducted using our numerical Bloom filter module NBF.pyc (the source
code for this module is also included - NBF.py).

The software uses some standard Python modules as well as two
non-standard Python modules that need to be installed. They are:
a. numpy
b. bitarray
In addition, it uses some Febrl modules which are auxiliary.pyc, 
encode.pyc, mymath.pyc, and stringcmp.pyc. These modules are 
included in the software.

3. results:
Finally, the efficiency (runtime and memory size) and effectiveness
(precision, recall, F-measure, and Spearman correlation coefficient)
results are written into a log file under the results folder - named as
'%dataset_%encode_%m' %(dataset name, encoding method, m).
The first line writes the parameter values used for the experiments
and the efficiency results are written in the second line of the file.
Effectiveness results are written in the third line onwards -
one line of results per query record.

License:
==================================================================
The PPSPM software is licensed under the GNU General Public License,
version 3 (GPL-3.0).
