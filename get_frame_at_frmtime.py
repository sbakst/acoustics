'''
Use acoustics file to find ultrasound frame.
First run get_low_f3_



'''
import sys
import os
import re
import argparse
import subprocess
from subprocess import call
import audiolabel
import csv
import pandas
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


# from stackoverflow

# define masks

def create_circular_mask(h, w, center=None, radius=None, shape = None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])
    if shape is None:
        shape = 'circle'

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    if shape == 'circle':
        mask = dist_from_center <= radius
    elif shape == 'semi':
        mask = dist_from_center >=radius
    return mask




parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Experiment directory containing all subjects' bmps and jpgs")
parser.add_argument("acoustdif", help="Path to TARGacoustic_data.txt")
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


subbmpdir = os.path.join(args.directory,args.subject)
acoustfi = TARG+'acoustic_data.txt'
afi = os.path.join(args.acoustdir,acoustfi)

# open one image in subbmpdir to get size

q = 0
dirind = 0
while q == 0: # account for possibility that dirs might be empty won't be any good
    try:
        somebmptime = os.listdir(subbmpdir)[dirind]
        somebmpdir = os.path.join(subbmpdir,somebmptime)
        subjframelist = re.compile('.*\.jpg')
        dirjpgs = [i for i in os.listdir(somebmpdir) if (subjframelist.search(i)]       
        somebmp = dirjpgs[0]
        im = np.array(Image.open(somebmp))
        q,s = im.shape[0:2]
     except IndexError as e:
        dirind = dirind+1    

# with open(afi) as afi_csv:
#     afi_reader = csv.reader(afi_csv, delimiter=',')
#     d = list(afi_reader)
#     rows = sum(1 for row in d)
times = []
stims = []
framenum = []
frames = None
frames = np.empty([len(os.listdir(subbmpdir))]+list(im.shape[0:2])) * np.nan  

ac_df = pandas.read_csv(afi,
    index_col = 'timestamp'
    )    
# get length of csv
    
for index, row in ac_df.iterrows():
    print(index, row)
    ts = row['timestamp']
    frmtime = row['frmtime']
    stim = row['stim']
    # get sync file
    utt = os.path.join(args.acoustdir, ts)
    syncfile = os.path.join(utt,(ts+'.bpr.sync.txt'))
    if not os.path.isfile(syncfile):
        print("can't find syncfile")
        continue
    try:
        sm = audiolabel.LabelManager(from_file = syncfile, from_type = 'table', sep = '\t', fields_in_head = False, fields = 't1,frameidx')
        frmnum = int(sm.tier('frameidx').label_at(rt_frmtime).text)
    except ValueError:
        syncfile = os.path.join(utt, (timestamp + '.bpr.sync.TextGrid')) # god you're so smart
        # print(syncfile)
        sm = audiolabel.LabelManager(from_file = syncfile, from_type='praat')#,sep = '\t', fields_in_head = True, fields = 'seconds,pulse_idx,raw_data_idx')
        frmnum = int(sm.tier('pulse_idx').label_at(rt_frmtime).text)
    # add frame to array    
    bmpdir = os.path.join(subbmpdir,ts)
    subjframelist = re.compile('.*\.jpg')
    dirjpgs = [i for i in os.listdir(bmpdir) if (subjframelist.search(i)]
    frmind = frmnum - 1
    testinds = [frmind, frmind-1, frmind+1]# one-frame tolerance
    for t = 1:len(testinds):
        theframe = dirjpgs[]
    
    
    # check for a bad frame
        
        
        