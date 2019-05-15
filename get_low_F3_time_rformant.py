# PCA with bmps
# based off Susan's PCA (reading in bmps) and Matt's but also my frame-subtracting technique for finding R steady state
# assumes bmps have already been made

import sys
import os
import re
import argparse
import subprocess
from subprocess import call
import audiolabel
import numpy as np
import shutil
import datetime
from PIL import Image
from scipy import misc
from scipy import ndimage
from scipy import sparse
import glob
import matplotlib.pyplot as plt
from ultratils.exp import Exp
from itertools import *

from sklearn import decomposition
from sklearn.decomposition import PCA



parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Experiment directory containing all subjects' bmps and jpgs")
parser.add_argument("otherdir", help="Experiment directory containing all subjects' tgs and sync files and wavs etc")
parser.add_argument("wordfi", help = "Txt file containing all words to search for")
parser.add_argument("TARG", help = "target phone in arpabet")
parser.add_argument("subject", help="subjnumber")

# parser.add_argument("-n", "n_components", help="Number of principal components to output", action="store_true")
parser.add_argument("-v", "--visualize", help="Produce plots of PC loadings on fan",action="store_true")
parser.add_argument("-p", "--pca", help="run pca",action="store_true")

args = parser.parse_args()

wordfi = args.wordfi
TARG = args.TARG


w = open(wordfi)
WORDS = str.split(w.read(),"\n")
print(WORDS)


# PARENTDIR = sys.argv[1]
# DESTDIR = sys.argv[2]
# SUB = sys.argv[3]

subbmpdir = os.path.join(args.directory,args.subject)
utildir = os.path.join(args.otherdir,args.subject)

# required vectors:
# vowel?
# frame (array)
framenos = [] # frame number
frmtimes = [] # when the frame occurs from rformant
stimulus = [] # stimulus word
ts = []
frmf3 = []# value of F3
frmf2= []
frmf1 = []
# First, for each directory of bmps, read in the stim.txt.

rpca = None
mrpca = None


for dirs, times, files in os.walk(subbmpdir):
# When you see this a year from now, you can remember sitting in the Illini union. ;)
    for tindex, timestamp in enumerate(times):
# append some kind of placeholder for when you skip a file
        subjrframelist = []
        utt = os.path.join(utildir,str(timestamp))
        syncfile = os.path.join(utt,(timestamp+'.bpr.sync.txt'))
       # print(syncfile)
        if not os.path.isfile(syncfile):
            print("can't find syncfile")
#        print(syncfile)
#        if not syncfile:
            continue
        stimfile = os.path.join(utt, 'stim.txt')
        stim = open(stimfile)
        stimtext = stim.read()
        print (stimtext)
        stimmy = stimtext[6:]
        print(stimmy)
        if stimmy in WORDS:
#            print('yay!')
#        if any(substring in stimtext for substring in WORDS): 
            tg = os.path.join(utt,(timestamp+'.TextGrid')) # may need to change to just TextGrid depending on when data is from
            wav = os.path.join(utt,(timestamp+'.wav'))
            if not os.path.isfile(tg):
                tg = os.path.join(utt,(timestamp+'.bpr.ch1.TextGrid'))
            if not os.path.isfile(wav):
                wav = os.path.join(utt,(timestamp+'.bpr.ch1.wav')) # may need to change to bpr.ch1.wav depending on when data is from
            # syncfile = os.path.join(utt,(timestamp+'.bpr.sync.txt'))
       #     if not syncfile:
       #         continue

            # get midpoint time
            
            pm = audiolabel.LabelManager(from_file = tg, from_type = 'praat')
            for lab in pm.tier('phone') :
                if (re.match(TARG,lab.text)) :
                    print(lab.text)
                    label = lab.text
                    rt1 = lab.t1
                    print(rt1)
                    rt2 = lab.t2
                    print(rt2)
                    rt_frmtime = ((rt1 + rt2)/2)
            subjrframelist = re.compile('.*\.bmp')
            regex = re.compile('(pc[0-9])|(mean).jpg')
            stimex = re.compile(stimtext) 
            bmpdir = os.path.join(subbmpdir, timestamp)
            imlist = [i for i in os.listdir(bmpdir) if (subjrframelist.search(i) and not regex.search(i) and not stimex.search(i))]
            # print(imlist)
#            try:
#                im = np.array(Image.open(os.path.join(bmpdir,(imlist[0])))) #open one image to get the size
#            except IndexError as e:
#                myrframe = 'NA'
#                continue
            if len(imlist) > 3:
                # run rformant    
                ts.append(timestamp)
                stimulus.append(str(stimmy))
                
                    
                fbfilename = (os.path.splitext(wav)[0]+'.fb') # make fbfile name that rformant will write to
                fbfile = os.path.join(utt,fbfilename) 

                try:
                    formant_proc = subprocess.check_call(["rformant", wav])#, stdin=rdc_proc.stdout) # also remove 20
                except subprocess.CalledProcessError as e:
                    print(e)
                    print(e.stdout)
                ppl_proc = subprocess.Popen(
                    ["pplain", fbfile],
                    stdout=subprocess.PIPE)
                print(fbfile)   
 			
                lm = audiolabel.LabelManager(
                   from_file=ppl_proc.stdout,   # read directly from pplain output
                   from_type='table',
                   sep=" ",
                   fields_in_head=False,
                   fields="f1,f2,f3,f4",
                   t1_col=None,                 # esps output doesn't have a t1 column
                   t1_start=0.0,          # autocreate t1 starting with this value and
                   t1_step=0.01)             # increase by this step
                f3_list = []
                vt1_list = []
                for v in lm.tier('f3'):
                    if rt1 <= v.t1 <= rt2:
                        print(v.t1)
                        print(v.text)
                        f3_list.append(float(v.text))
                        vt1_list.append(v.t1)
                f3ind = np.argmin(f3_list)
                print(f3ind)
                print(vt1_list)
                f3time = vt1_list[f3ind]
                print(f3_list)
                f3min = np.min(f3_list)
                if (stimmy == 'beer' or stimmy == 'ream'):
                    if f3min > 2400 or f3min < 1500:
                        f3min = (float('nan'))
                elif f3min > 2300 or f3min < 1300:
                    f3min = (float('nan'))
                
                frmf3.append(f3min)
                frmtimes.append(f3time)

if np.isnan(frmf3).any():
    datafiname = TARG+'acoustic_data_requires_attention.txt'
else:                 
    datafiname = TARG+'acoustic_data.txt'
datafi = os.path.join(subbmpdir,datafiname)

#da = open(datafi, 'w')
da_headers = ["timestamp","frmf3","frmtime","stim"]
h = np.row_stack((da_headers,np.column_stack((ts,frmf3,frmtimes,stimulus))))
np.savetxt(datafi, h, fmt = "%s", delimiter = ",") 
                
