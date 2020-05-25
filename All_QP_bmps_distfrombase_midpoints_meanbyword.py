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
import matplotlib
import matplotlib.pyplot as plt
from ultratils.exp import Exp
import csv
from itertools import *
from numpy import nan

from sklearn import decomposition
from sklearn.decomposition import PCA


# from stackoverflow

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

#def semi_circ_mask(h, w, center = None, radius = None, side = None):
#    if center is None:
#        center = [int(w/2), int(h/2)]
#    if radius is None: # use the smallest distance between the center and image walls
#        radius = min(center[0], center[1], w-center[0], h-center[1])


    # for hollow mask
#     mask = dist_from_center >= radius
#     return mask    



# this could be made a separate file I read in



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

matplotlib.use('macosx',warn=False, force=True)


wordfi = args.wordfi
TARG = args.TARG


w = open(wordfi)
WORDS = str.split(w.read(),"\n")
print(WORDS)

subbmpdir = os.path.join(args.directory,args.subject)
utildir = os.path.join(args.otherdir,args.subject)

cx = None
rds = None
ux = None
hrad = None
hq = None


omaskfi = os.path.join(subbmpdir,'Smaskparams.txt')
if os.path.isfile(omaskfi):
    mf = open(omaskfi)
    mfr = csv.reader(mf,delimiter='\t')
    rows = [r for r in mfr]
    mparams = rows[1]
    print(mparams)
    cx = int(float(mparams[0]))
    rds= int((mparams[3]))
    ux = int(float(mparams[1]))
    hrad = int((mparams[4]))
    hq = int((mparams[2])) 
    
nmaskfi = os.path.join(subbmpdir,TARG + 'mbw_maskparams.txt')    
if os.path.isfile(nmaskfi):
    mf = open(nmaskfi)
    mfr = csv.reader(mf,delimiter='\t')
    rows = [r for r in mfr]
    mparams = rows[1]
    print(mparams)
    cx = int(float(mparams[0]))
    rds= int(float(mparams[1]))
    ux = int((mparams[2]))
    hrad = int((mparams[3]))
    hq = int((mparams[4])) 
    

framenos = None # frame number
frmtimes = [] # when the frame occurs from sync.txt
stimulus = [] # stimulus word
ts = []



# First, for each directory of bmps, read in the stim.txt.

# rpca = None
mpca = None


for dirs, times, files in os.walk(subbmpdir):
# When you see this a year from now, you can remember sitting in the Illini union. ;)
    for tindex, timestamp in enumerate(times):
# append some kind of placeholder for when you skip a file
        print(timestamp)
        ts.append(timestamp)
        subjrframelist = []
        utt = os.path.join(utildir,str(timestamp))
        syncfile = os.path.join(utt,(timestamp+'.bpr.sync.txt'))
        stimfile = os.path.join(utt, 'stim.txt')
        if os.path.isfile(stimfile):
            stim = open(stimfile)
            stimtext = stim.read()
            print (stimtext)
            stimmy = stimtext[6:]
            stimmy = re.sub("[^a-zA-Z]+", "", stimmy)        
            print(stimmy)
            print('no'+ stimmy +'yes')
            if int(args.subject) > 120:
                stimmy = stimtext
            stimulus.append(str(stimmy))
        else:
            stimulus.append('doh!')
            continue
       # print(syncfile)
        if not os.path.isfile(syncfile):
            print("can't find syncfile")
            if not os.path.isfile(os.path.join(utt, (timestamp + '.bpr.sync.TextGrid'))):
#        print(syncfile)
#        if not syncfile:
                continue
        if stimmy in WORDS:
            print('yay!')
#        if any(substring in stimtext for substring in WORDS): 
            tg = os.path.join(utt,(timestamp+'.TextGrid')) # may need to change to just TextGrid depending on when data is from
            wav = os.path.join(utt,(timestamp+'.wav'))
            if not os.path.isfile(tg):
                tg = os.path.join(utt,(timestamp+'.bpr.ch1.TextGrid'))
            if not os.path.isfile(wav):
                wav = os.path.join(utt,(timestamp+'.bpr.ch1.wav')) # may need to change to bpr.ch1.wav depending on when data is from
            # print(wav)
            # syncfile = os.path.join(utt,(timestamp+'.bpr.sync.txt'))
       #     if not syncfile:
       #         continue

            # get midpoint time
            pm = audiolabel.LabelManager(from_file = tg, from_type = 'praat')
            for lab in pm.tier('phone') :
                if (re.match(TARG,lab.text)) :
                    print(lab.text)
                    label = lab.text
                    t1 = lab.t1
                    t2 = lab.t2
                    mid_frmtime = ((t1 + t2)/2)
            subjrframelist = re.compile('.*\.jpg')
            regex = re.compile('(pc[0-9])|(mean).jpg')
            stimex = re.compile(stimtext) 
            bmpdir = os.path.join(subbmpdir, timestamp)
            imlist = [i for i in os.listdir(bmpdir) if (subjrframelist.search(i) and not regex.search(i) and not stimex.search(i))]
            # print(imlist)
            try:
                im = np.array(Image.open(os.path.join(bmpdir,(imlist[0])))) #open one image to get the size
            except IndexError as e:
                mymidframe = 'NA'
                continue
            if len(imlist) > 3:
                q,s = im.shape[0:2] #get the size of the images
                imnbr = len(imlist) #get the number of images
                if mpca is None:
                    mpca = np.empty([len(os.listdir(subbmpdir))]+list(im.shape[0:2])) * np.nan
                    framenos = np.empty([len(os.listdir(subbmpdir))])*np.nan
                try:
                    sm = audiolabel.LabelManager(from_file = syncfile, from_type = 'table', sep = '\t', fields_in_head = False, fields = 't1,frameidx')
                    midfrm = int(sm.tier('frameidx').label_at(mid_frmtime).text)
                    print('oonoo')
                except ValueError:
                    sm = audiolabel.LabelManager(from_file = syncfile, from_type='table',sep = '\t', fields_in_head = True, fields = 'seconds,pulse_idx,raw_data_idx')
                    midfrm=int(sm.tier('pulse_idx').label_at(mid_frmtime).text)
                    print('eenee')
                except ValueError:
                    syncfile = os.path.join(utt, (timestamp + '.bpr.sync.TextGrid'))
                    os.path.isfile(syncfile)
                    sm = audiolabel.LabelManager(from_file = syncfile, from_type='praat')#,sep = '\t', fields_in_head = True, fields = 'seconds,pulse_idx,raw_data_idx')
                    midfrm = int(sm.tier('pulse_idx').label_at(mid_frmtime).text)
                     
#                except IndexError:
#                        
#                    syncfile = os.path.join(utt, (timestamp + '.bpr.sync.TextGrid'))
#                    print(syncfile)
#
#                    sm = audiolabel.LabelManager(from_file = syncfile, from_type='praat')#,sep = '\t', fields_in_head = True, fields = 'seconds,pulse_idx,raw_data_idx')
#                    midfrm=int(sm.tier('pulse_idx').label_at(mid_frmtime).text)
                print(midfrm)
                frmind = midfrm
                
                testinds = [frmind, frmind-1, frmind+1]# one-frame tolerance
                for t in range (0,len(testinds)):
                    thisind = testinds[t]
                    theframe = timestamp + '.' + str(thisind) + '.jpg'
                    # theframe = dirjpgs[thisind]
                    try:
                        openframe = np.array(Image.open(os.path.join(bmpdir,theframe)))
                    except FileNotFoundError:
                        continue
                    brightness = np.linalg.norm(openframe)
                    if brightness > 25000: # bad frame
                        theframe = float(nan)
                        openframe = float(nan)
                        continue
                    else:
                        break
                mpca[tindex,:,:] = openframe   
                framenos[tindex]=(thisind+1)     
#                
#                print('that was the midframe')
#                framename = timestamp + '.' + str(midfrm) + '.jpg'
#                theframe = os.path.join(bmpdir,framename)
#                print(theframe)
#                if not (os.path.isfile(theframe)):
#                    midfrm = midfrm+1
#                    print(midfrm)
#                    framename = timestamp + '.' + str(midfrm) + '.jpg'
#                    theframe = os.path.join(bmpdir,framename)
#                if np.linalg.norm(np.array(Image.open(theframe)))> 25000:
#                    midfrm = midfrm + 1
#                    if np.linalg.norm(np.array(Image.open(midfrm)))> 25000:
#                        midfrm = midfrm - 2
#                        if np.linalg.norm(np.array(Image.open(midfrm)))> 25000:
#                            mpca[tindex,:,:] = nan
#                            continue
#                framename = timestamp + '.' + str(midfrm) + '.jpg'
#                framenos.append(midfrm)
#                openframe = np.array(Image.open(os.path.join(bmpdir,framename)))
#                # framenos.append(np.linalg.norm(openframe))
#                # mymidframe = ndimage.median_filter(openframe,5)
#                mpca[tindex,:,:] = openframe #mymidframe
#

print(len(mpca))
print(len(framenos))

# print(keep_indices)
print(len(ts))
frames = np.squeeze(mpca)
keep_indices = np.where(~np.isnan(frames).any(axis=(1,2)))
naninds = np.where(np.isnan(frames).any(axis=(1,2)))
print(naninds)
print(keep_indices)
times = np.squeeze(ts)
# print(times)
stims = np.squeeze(stimulus)
print(stims)
print(len(stims))
framenums = np.squeeze(framenos)
#print(framenums)
#print(framenos)

kept_frames = frames[keep_indices]
kept_times = times[keep_indices]
kept_stims = stims[keep_indices]
kept_framenum = framenums[keep_indices]            
        

# test and add masks

testframe = kept_frames[5]
print(kept_times[5])
if rds is None:
    rds = 65 # test radius
if cx is None:
    cx = s/2 # start center for lower mask
print(cx)
print(rds)
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

if ux is None:
    ux = s/2
if hrad is None:
    hrad = 180
success = 0
while success == 0 :
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h,w, center = [ux,q], radius = hrad, shape = 'semi')
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    resp = input('Is the mask centered? (y, tl, tr, sug)')

    if resp == 'y' : 
        success = 1
    elif resp == 'tr':
        ux = ux - 5
    elif resp == 'tl':
        ux = ux + 5
    elif resp == 'sug':
        ux = input('Enter an x-coordinate. Currently displayed is ' + str(ux))
        ux = int(ux)
plt.close()



success = 0
while success == 0 :
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [ux,q], radius = hrad, shape = 'semi')
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    resp = input('Is the mask too wide, too narrow, or juuuust right?(tw, tn, gdx, sug)')
    if resp == 'tw':
        hrad = hrad - 10
    elif resp == 'tn':
        hrad = hrad + 10
    elif resp == 'sug':
        hrad = input('Enter a number for the radius. On the screen is ' + str(hrad))
        hrad = int(hrad)
        print(hrad)
    elif resp == 'gdx':
        success = 1
    plt.close()
plt.close()
print(hrad)
testframe = masked_img

success = 0
if hq is None:
    hq = 180
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [ux, hq], radius = hrad, shape = 'semi')
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    ht = input('Is the mask too high, too low, or juuust right?(th, tl, gdx, sug)')
    if ht == 'th':
       hq = hq+10
    elif ht == 'tl':
        hq = hq-10
    elif ht == 'sug':
        hq = input('Enter a number to add to the height; larger = lower. On the screen is ' +str(hq))
        hq = int(hq)
    elif ht == 'gdx':
        success = 1

    plt.close()

# add semi-circular mask
for i in range(0,(kept_frames.shape[0])):
    maskframe = kept_frames[i]
    mask = create_circular_mask(h, w, center = [ux,hq], radius = hrad, shape = 'semi')
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

### save files


    outfiname = bword + '_mbw_mask_diffsfrombase_mid.txt'
    outdiffs = os.path.join(subbmpdir, outfiname)
    od = open(outdiffs, 'w')
    od.write('\t'.join(['subject','norm_avg','normalized_avg','raw_normavg','raw_normavg_normed','avg_brightness'])+'\n')
    od.write('\t'.join([args.subject,str(norm_avg),str(normalized_avg),str(rawdiff_avg_norm),str(rawdiff_avg_norm_normed),str(normavgfrm)])+'\n')
    od.close()
    
    wordtimes = kept_times[wordinds]
    wordframenum = kept_framenum[wordinds]

    # separate file with by-timestamp info
    byfiname = bword + '_midpt_data.txt'
    dfi = os.path.join(subbmpdir, byfiname)
    data_headers = ["timestamp","framenum","normdiff","normeddiff"] 
    h = np.row_stack((data_headers,np.column_stack((wordtimes,wordframenum,norm_list,normalized_list))))
    np.savetxt(dfi,h,fmt="%s",delimiter = ",")

# next time add word >.>
    
maskfi = TARG + 'mbw_maskparams.txt'
maskparams = os.path.join(subbmpdir, maskfi)
mp = open(maskparams, 'w')
mp.write('\t'.join(['posx','low_radius','highx','high_radius','height'])+'\n')
mp.write('\t'.join([str(cx), str(rds), str(ux),str(hrad),str(hq)]))
mp.close()


# do this to print midpoint frame and some diff frames

# pic=avgfrm
# mag = np.max(pic)-np.min(pic)
# pic = (pic-np.min(pic))/mag*255
# printy = misc.imsave(os.path.join(subbmpdir,'AverageS.png'),pic)
# 
# 
# # raw_list
# 
# # image_shape = (q,s)
# for r in range(0, len(raw_list):
#    pic = raw_list[r]
#    mag = np.max(pic) - np.min(pic)
#    pic = (pic-np.min(pic))/mag*255
#    pcn = misc.imsave(os.path.join(subbmpdir,'dfb'+str(r+1)+'.png'),pic)
# 




# print(framenos)

if args.pca:
####################################################################################################


# remove any indices for all objects generated above where frames have NaN values (due to skipping or otherwise)
# ooooo she's squeezin'!


    n_components = 5
    
    pca = PCA(n_components = n_components)
    
    print(kept_mpca.shape[0]) # 48ish?
    print(kept_mpca.shape[1]) # q?
    print(kept_mpca.shape[2]) # s?
    
    
    frames_reshaped = kept_mpca.reshape([kept_mpca.shape[0], kept_mpca.shape[1]*kept_mpca.shape[2]])
    
    sparse_frames = sparse.csr_matrix(frames_reshaped)
    
    # pca with sparse matrix
    pca.fit(frames_reshaped)
    analysis = pca.transform(frames_reshaped)
    #covariance = pca.get_covariance()
    
    
    meta_headers = ["timestamp","framenum","f1","f2","f3"]
    pc_headers = ["pc"+str(i+1) for i in range(0,n_components)] # determine number of PC columns; changes w.r.t. n_components
    headers = meta_headers + pc_headers
    print(analysis)
    print(headers)
    #print('covariance' + covariance)
    
    #print(kept_stimulus)
    print(kept_ts)
    print(kept_framenos)
    
    outfile = os.path.join(subbmpdir, 'pca.csv')
    d = np.row_stack((headers,np.column_stack((kept_ts,kept_framenos,kept_f1,kept_f2,kept_f3,analysis))))
    np.savetxt(outfile, d, fmt="%s", delimiter =',')
    
    expl = os.path.join(subbmpdir, 'explain.txt')
    f = open(expl, 'w')
    f.write('explained variance'+ str(pca.explained_variance_)+'\n')
    f.write('explained ratio: '+ str(pca.explained_variance_ratio_))
    f.write('\n')
    f.write('mean: ' + str(pca.mean_))
    f.write('\n')
    #f.write(str(covariance))
    f.close()
    
    # make pc pictures
    if args.visualize:
        image_shape = (q,s)
        for n in range(0, n_components):
            pic = pca.components_[n].reshape(image_shape)
            print ('look tis pic')
            print (pic)
            mag = np.max(pic) - np.min(pic)
            pic = (pic-np.min(pic))/mag*255
            pcn = misc.imsave(os.path.join(subbmpdir,'pc'+str(n+1)+'.png'),pic)
    


#r_mean = np.kept_mpca.mean(0)
#for i in np.arange(length(kept_mpca[0])):
#   std=(x-mu)


#standev = 
#f = open(





