## split filename, sort by trial, and make dictionary with ultrasound ##

import os
import re


# make list of tuples

tuplist = []
for dirs, subdir, files in os.walk(acoustics_dir):
   for wav in subdir:
      tup = re.split('_|.wav',wav)
      tuplist.append(tup)

# sort tuples by trial number, index 3

ordtrials = []
ordtrials = sorted(tuplist, key = lambda x: x[3])

# make list of ultrasounds

ultras = []
for dirs, subdir, files in os.walk(ultrasound_dir):
   for tp in subdir:
# split up timestamp
      tptup = re.split('-',tp)
      ultras.append(tp)

# split order by 3rd element, time

ordultras = []
ordultras = sorted(ultras, key = lambda x: x[2])

# make dictionary

tamildict = dict(zip(ordtrials,ordultras))


