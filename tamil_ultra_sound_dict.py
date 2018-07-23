## split filename, sort by trial, and make dictionary with ultrasound ##

import os
import sys
import re
import operator

acoustics_dir = sys.argv[1]
ultrasound_dir = sys.argv[2]

# make list of tuples

tuplist = []
for dirs, subdir, files in os.walk(acoustics_dir):
   for wav in files:
      if os.path.splitext(wav)[1]=='.wav':
         tup = re.split('_|.wav',wav)
         tuplist.append(tuple(tup))
#print(tuplist)
# sort tuples by trial number, index 3

ordtrials = []
ordtrials = sorted(tuplist, key = lambda x: int(x[3]))

# make list of ultrasounds

ultras = []
for dirs, subdir, files in os.walk(ultrasound_dir):
   for tp in subdir:
# split up timestamp
      tptup = re.split('-',tp)
      ultras.append(tuple(tptup))

# split order by 3rd element, time

ordultras = []
ordultras = sorted(ultras, key = lambda x: x[2])

# make dictionary

tamildict = dict(zip(ordtrials,ordultras))

dictitems = tamildict.items()
tamtam = (sorted(dictitems, key= lambda x:int(x[0][3])))

print(tamtam)

for entry in tamtam:
   soundname = ('_'.join(entry[0]))[:-1]
   wav = os.path.join(acoustics_dir, (soundname + '.wav'))
   tg = os.path.join(acoustics_dir, (soundname + '.TextGrid'))
   # now for the timestamp
   timedir = os.path.join(ultrasound_dir,('-'.join(entry[1])))



