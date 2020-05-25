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
parser.add_argument("acoustdir", help="Path to TARGacoustic_data.txt")
parser.add_argument("wordfi", help = "Txt file containing all words to search for")
parser.add_argument("TARG", help = "target phone in arpabet")
parser.add_argument("bigbmpdir", help = "where to find all the bmps for this subject, if any are missing")
parser.add_argument("subject", help="subjnumber")

# parser.add_argument("-n", "n_components", help="Number of principal components to output", action="store_true")
parser.add_argument("-v", "--visualize", help="Produce plots of PC loadings on fan",action="store_true")
parser.add_argument("-p", "--pca", help="run pca",action="store_true")

args = parser.parse_args()

wordfi = args.wordfi
TARG = args.TARG
bigbmpdir = args.bigbmpdir

w = open(wordfi)
WORDS = str.split(w.read(),"\n")
print(WORDS)


subbmpdir = os.path.join(args.directory,args.subject)
acoustfi = TARG+'acoustic_data.txt'
afi = os.path.join(os.path.join(args.acoustdir,args.subject),acoustfi)
omaskfi = os.path.join(subbmpdir, TARG + 'mbw_maskparams.txt')

cx = None
rds = None
q = None
hrad = None

if os.path.isfile(omaskfi):
    mf = open(omaskfi)
    mfr = csv.read(mf)
    rows = [r for r in mfr]
    mparams = rows[1]
    cx = mparams[0]
    rds = mparams[1]
    yht = mparams[2]
    hrad = mparams[3]
#mp.write('\t'.join(['low_ctr_x','low_radius','high_ctr_y','high_radius'])+'\n')
#mp.write('\t'.join([str(cx), str(rds), str(q),str(radius)]))



# open one image in subbmpdir to get size

q = 0
dirind = 0
while q == 0: # account for possibility that dirs might be empty won't be any good
    try:
        utterances = (glob.glob(os.path.join(subbmpdir,'201*')))
        somebmptime = utterances[dirind]
        # print(somebmptime)
        somebmpdir = os.path.join(subbmpdir,somebmptime)
        subjframelist = re.compile('.*\.jpg')
        # print(subjframelist)
        dirjpgs = [i for i in os.listdir(somebmpdir) if (subjframelist.search(i))]       
        # print(dirjpgs)
        somebmp = dirjpgs[0]
        im = np.array(Image.open(os.path.join(os.path.join(subbmpdir,somebmptime),somebmp)))
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
frames = np.empty([len(utterances)]+list(im.shape[0:2])) * np.nan  

ac_df = pandas.read_csv(afi)
# get length of csv
    
for index, row in ac_df.iterrows():
    # print(index, row)
    ts = row['timestamp']
    print(ts)
    frmtime = row['frmtime']
    stim = row['stim']
    # get sync file
    utt =os.path.join(os.path.join(args.acoustdir,args.subject), ts)
    syncfile = os.path.join(utt,(ts+'.bpr.sync.txt'))
    if not os.path.isfile(syncfile):
        print("can't find syncfile")
        times.append(ts)
        stims.append(stim)
        framenum.append(float('nan'))
        continue
    try:
        sm = audiolabel.LabelManager(from_file = syncfile, from_type = 'table', sep = '\t', fields_in_head = False, fields = 't1,frameidx')
        frmnum = int(sm.tier('frameidx').label_at(frmtime).text)
    except ValueError:
        print('try this')
        syncfile = os.path.join(utt, (ts + '.bpr.sync.TextGrid')) # god you're so smart
        print(syncfile)
        print(frmtime)
        sm = audiolabel.LabelManager(from_file = syncfile, from_type='praat')#,sep = '\t', fields_in_head = True, fields = 'seconds,pulse_idx,raw_data_idx')
        frmnum = int(sm.tier('pulse_idx').label_at(frmtime).text)
    # add frame to array    
    bmpdir = os.path.join(subbmpdir,ts)
#    subjframelist = re.compile('.*\.jpg')
#    dirjpgs = [i for i in os.listdir(bmpdir) if (subjframelist.search(i))]
    frmind = frmnum# - 1
     
    stims.append(stim)
    times.append(ts)
    # create name of file
    testinds = [frmind, frmind-1, frmind+1]# one-frame tolerance
    for t in range (0,len(testinds)):
        thisind = testinds[t]
        theframe = ts + '.' + str(thisind) + '.jpg'
        # theframe = dirjpgs[thisind]
        framepath = os.path.join(bmpdir,theframe)
        if not os.path.isfile(framepath):
            frompath = os.path.join(os.path.join(bigbmpdir,args.subject),ts)
            frame2copy = os.path.join(frompath,(ts + '.' + str(thisind) + '.bmp'))
            if not os.path.isfile(frame2copy):
                break
            else:
                shutil.copy(frame2copy, bmpdir)
                newbmp = os.path.join(bmpdir,(ts+'.'+str(thisind)+'.bmp'))
                proc = subprocess.check_output(['mogrify', '-format', 'jpg', newbmp]) #, shell=True])
        openframe = np.array(Image.open(framepath))
        brightness = np.linalg.norm(openframe)
        if brightness > 25000: # bad frame
            openframe = float('nan')
            continue
        else:
            break
    frames[index,:,:] = openframe   
    framenum.append(thisind+1)     
#    arrind = arrind + 1

keep_indices = np.where(~np.isnan(frames).any(axis=(1,2)))[0]
# print(keep_indices)
print(len(times))
frames = np.squeeze(frames)
times = np.squeeze(times)
# print(times)
stims = np.squeeze(stims)
framenum = np.squeeze(framenum)
kept_frames = frames[keep_indices]
kept_times = times[keep_indices]
kept_stims = stims[keep_indices]
kept_framenum = framenum[keep_indices]            
        

# test and add masks

testframe = kept_frames[5]

if rds is None:
    rds = 65 # test radius
if cx is None:
    cx = s/2 # start center for lower mask


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
    resp = input('Is the mask centered? (y, tr, tl, sug)')
    if (resp == 'y'):#  | resp == 'gdx'):
        success = 1
    elif resp == 'tr':
        cx = cx -5
    elif resp == 'tl':
        cx = cx + 5
    elif resp == 'sug':
        cx = int(input('Enter an x-coordinate. Currently displayed is ' + str(cx)))
        print(cx)
        
success = 0
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [cx,h], radius = rds)
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    resp = input('Is the mask too big, too small, or juuuuust right? (tb, ts, gdx, sug)')
    if resp == 'tb':
        rds = rds - 3
    elif resp == 'ts':
        rds = rds + 7
    elif resp == 'sug':
        rds = input('Enter a number for the radius. On the screen is ' + str(rds))
        rds = int(rds)        
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
success = 0
if hrad is None:
    hrad = 180
while success == 0 :
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [s/2,q], radius = hrad, shape = 'semi')
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    resp = input('Is the mask too wide, too narrow, or juuuust right?(tw, tn, gdx, sug)')
    if resp == 'tw':
        hrad = hrad - 10
    elif resp == 'tn':
        hrad =hrad + 10
    elif resp == 'sug':
        hrad = input('Enter a number for the radius. On the screen is ' + str(radius))
        hrad = int(radius)
        print(hrad)
    elif resp == 'gdx':
        success = 1
    plt.close()
plt.close()
print(hrad)
testframe = masked_img

success = 0
if yht is None:
    yht = 180
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [s/2, yht], radius = hrad, shape = 'semi')
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    ht = input('Is the mask too high, too low, or juuust right?(th, tl, gdx, sug)')
    if ht == 'th':
       yht = yht+10
    elif ht == 'tl':
        yht = yht-10
    elif ht == 'sug':
        yht = input('Enter a number to add to the height; larger = lower. On the screen is ' +str (yht))
        yht = int(yht)
    elif ht == 'gdx':
        success = 1

    plt.close()

# add semi-circular mask
for i in range(0,(kept_frames.shape[0])):
    maskframe = kept_frames[i]
    mask = create_circular_mask(h, w, center = [s/2,yht], radius = hrad, shape = 'semi')
    masked_img = maskframe.copy()
    masked_img[mask] = 0
    if i == 1:
        plt.imshow(masked_img, cmap = "Greys_r")
        plt.show()
    kept_frames[i] = masked_img
plt.close()

# make six means
kept_stims = [item.lower() for item in kept_stims]
# iteratively do means and subtractions
# bwlists= [rahlist, romelist, reamlist, barlist, borelist, beerlist]  # by word lists

words = [item.lower() for item in WORDS] 
words=list(set(words))
numwords = len(words)
print('here are the words')
print(words)
print(numwords)
for n in np.arange(0,numwords):
    bword = words[n]
    wordinds = ([i for i, x in enumerate(kept_stims) if x == bword])
    print(wordinds)
    b = kept_frames[wordinds]
#    print(b)
    avgfrm = np.mean(b, axis = 0)
    normavgfrm = np.linalg.norm(avgfrm)
    



# now do the subtraction

# average frame
#avgfrm = np.mean(kept_frames, axis = 0)
# norm of the average frame
#normavgfrm = np.linalg.norm(avgfrm) # slight digression from original

    raw_list = []
    norm_list = []
    normalized_list = []

    for i in np.arange(len(b)):
        rawdiff = b[i,:,:] - avgfrm # matrix of difference
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
# TODO: ADD TO LIST

    

    outfiname = bword + '_mbw_mask_diffs_frmtime.txt'
    outdiffs = os.path.join(subbmpdir, outfiname)
    od = open(outdiffs, 'w')
    od.write('\t'.join(['subject','norm_avg','normalized_avg','raw_normavg','raw_normavg_normed','avg_brightness'])+'\n')
    od.write('\t'.join([args.subject,str(norm_avg),str(normalized_avg),str(rawdiff_avg_norm),str(rawdiff_avg_norm_normed),str(normavgfrm)])+'\n')
    od.close()
    
    wordtimes = kept_times[wordinds]
    wordframenum = kept_framenum[wordinds]

    # separate file with by-timestamp info
    byfiname = bword + '_mbw_mask_diffs_byts_frmtime.txt'
    byts = os.path.join(subbmpdir, byfiname)
    data_headers = ["timestamp","framenum","normdiff","normeddiff"] 
    b = np.row_stack((data_headers,np.column_stack((wordtimes,wordframenum,norm_list,normalized_list))))
    np.savetxt(byts,b,fmt="%s",delimiter = ",")



# mask params
maskfi = TARG + 'mbw_maskparams.txt'
maskparams = os.path.join(subbmpdir, maskfi)
mp = open(maskparams, 'w')
mp.write('\t'.join(['low_ctr_x','low_radius','high_ctr_y','high_radius'])+'\n')
mp.write('\t'.join([str(cx), str(rds), str(yht),str(hrad)]))
mp.close()
