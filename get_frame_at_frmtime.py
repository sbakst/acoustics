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
    for t in range (0,len(testinds)):
        thisind = testinds[t]
        theframe = dirjpgs[thisind]
        openframe = np.array(Image.open(os.path.join(bmpdir,theframe)))
        brightness = np.linalg.norm(openframe)
        if brightness > 25000: # bad frame
            theframe = float(nan)
            continue
        else:
            break
    frames[arrind,:,:] = theframe   
    times.append[ts]
    stims.append[stim]
    framenum.append[thisind+1]     
    arrind = arrind + 1

keep_indices = np.where(~np.isnan(frames).any(axis=(1,2)))[0]
kept_frames = frames[keep_indices]
kept_times = times[keep_indices]
kept_stims = stims[keep_indices]
kept_framenums = framenums[keep_indices]            
        

# test and add masks

testframe = kept_frames[5]
rds = 65 # test radius
cx = w/2 # start center for lower mask


# first find lower mask
success = 0
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [cx,h], radius = rds)
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
# ask opinion upon closing
    resp = input('Is the mask centered?')
    if resp == 'y' | if resp == 'gdx':
        success = 1
    elif resp 'tr':
        cx = cx -5
    elif resp == 'tl':
        cx = cx + 5
    elif resp == 'sug':
        cx = input('Enter an x-coordinate. Currently displayed is ' + str(cx))
        
success = 0
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [cx,h], radius = rds)
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    resp = input('Is the mask too big, too small, or juuuuust right? (tb, ts, gdx)')
    if resp == 'tb':
        rds = rds - 3
    elif resp == 'ts':
        rds = rds + 7
    elif resp == 'sug':
        rds = input('Enter a number for the radius. On the screen is ' + str(rds))
        rds = int(radius)        
    elif resp == 'gdx':
        success = 1
    print(rds)
    plt.close()


# add lower mask to all images

for i in range(0,(kept_frames.shape[0])):
    maskframe = kept_frames[i]
    mask = create_circular_mask(h, w, center = [cx,h], radius = rds)
    masked_img = maskframe.copy()
    masked_img[mask] = 0
    kept_frames[i] = masked_img

testframe = kept_frames[5]

# semi-circular upper mask
radius = 180
success = 0
while success == 0 :
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [s/2,q], radius = radius, shape = 'semi')
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    resp = input('Is the mask too wide, too narrow, or juuuust right?(tw, tn, gdx, sug)')
    if resp == 'tw':
        radius = radius - 10
    elif resp == 'tn':
        radius = radius + 10
    elif resp == 'sug':
        radius = input('Enter a number for the radius. On the screen is ' + str(radius))
        radius = int(radius)
    elif resp == 'gdx':
        success = 1
    plt.close()
plt.close()
print(radius)
testframe = masked_img

success = 0
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [s/2, q], radius = radius, shape = 'semi')
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    ht = input('Is the mask too high, too low, or juuust right?(th, tl, gdx, sug)')
    if ht == 'th':
       q = q+10
    elif ht == 'tl':
        q = q-10
    elif ht == 'sug':
        q = input('Enter a number to add to the height; larger = lower. On the screen is ' +str (q))
        q = int(q)
    elif ht == 'gdx':
        success = 1

    plt.close()

# add semi-circular mask
for i in range(0,(kept_frames.shape[0])):
    maskframe = kept_frames[i]
    mask = create_circular_mask(h, w, center = [s/2,q], radius = radius, shape = 'semi')
    masked_img = maskframe.copy()
    masked_img[mask] = 0
    if i == 1:
        plt.imshow(masked_img, cmap = "Greys_r")
        plt.show()
    kept_frames[i] = masked_img
plt.close()

# now do the subtraction

# average frame
avgfrm = np.mean(kept_frams, axis = 0)
# norm of the average frame
normavgfrm = np.linalg.norm(avgfrm) # slight digression from original

for i in np.arange(len(kept_frames)):
    rawdiff = kept_frames[i,:,:] - avgfrm # matrix of difference
    raw_list.append(rawdiff)
    normdiff = np.linalg.norm(rawdiff) # norm of matrix of difference
    norm_list.append(normdiff)
    normalizeddiff = normdiff/normavgfrm # norm of matrix of difference, normalized by the norm of the mean frame
    normalized_list.append(normalizeddiff)

# take averages of raw diff before norming
rawdiff_avg = np.mean(raw_list)
rawdiff_avg_norm = np.linalg.norm(rawdiff_avg)
rawdiff_avg_norm_normed = rawdiff_avg_norm/normavgfrm

# norm avg
norm_avg = np.mean(norm_list)
# -alized
normalized_avg = np.mean(normalized_list)


outfiname = TARG + '_mask_diffs_frmtime.txt'
outdiffs = os.path.join(subbmpdir, outfiname)
od = open(outdiffs, 'w')
od.write('\t'.join(['subject','norm_avg','normalized_avg','raw_normavg','raw_normavg_normed','avg_brightness'])+'\n')
od.write('\t'.join([args.subject,str(norm_avg),str(normalized_avg),str(rawdiff_avg_norm),str(rawdiff_avg_norm_normed),str(normavgfrm))+'\n')
od.close()

# separate file with by-timestamp info
byfiname = TARG + '_mask_diffs_byts_frmtime.txt'
byts = os.path.join(subbmpdir, byfiname)
data_headers = ["timestamp","framenum","normdiff","normeddiff"] 
b = np.row_stack((data_headers,np.column_stack((kept_ts,kept_framenums,norm_list,normalized_list))))
np.savetxt(byts,b,fmt="%s",delimiter = ",")


# mask params
maskfi = TARG + 'maskparams.txt'
maskparams = os.path.join(subbmpdir, maskfi)
mp = open(maskparams, 'w')
mp.write('\t'.join(['low_ctr_x','low_radius','high_ctr_y','high_radius'])+'\n')
mp.write('\t'.join([str(cx), str(rds), str(q),str(radius)]))
mp.close()
