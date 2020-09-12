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
from itertools import *

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


framenos = [] # frame number
frmtimes = [] # when the frame occurs from sync.txt
stimulus = [] # stimulus word
ts = []
frmf3 = []
frmf2= []
frmf1 = []
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
        if int(args.subject) > 120:
            stimmy = stimtext
        stimulus.append(str(stimtext))
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
                try:
                    sm = audiolabel.LabelManager(from_file = syncfile, from_type = 'table', sep = '\t', fields_in_head = False, fields = 't1,frameidx')
                    midfrm = int(sm.tier('frameidx').label_at(mid_frmtime).text)
                except ValueError:
                    syncfile = os.path.join(utt, (timestamp + '.bpr.sync.TextGrid'))
                    # print(syncfile)
                    sm = audiolabel.LabelManager(from_file = syncfile, from_type='praat')#,sep = '\t', fields_in_head = True, fields = 'seconds,pulse_idx,raw_data_idx')
#                    print(sm)
                    midfrm = int(sm.tier('pulse_idx').label_at(mid_frmtime).text)
                print(midfrm)
                print('that was the midframe')
                framename = timestamp + '.' + str(midfrm) + '.jpg'
                theframe = os.path.join(bmpdir,framename)
                if not (os.path.isfile(theframe)):
                    midfrm = midfrm+1
                    framename = timestamp + '.' + str(midfrm) + '.jpg'
                if np.linalg.norm(np.array(Image.open(theframe)))> 25000:
                    midfrm = midfrm + 1
                    framename = timestamp + '.' + str(midfrm) + '.jpg'
                framenos.append(midfrm)
                openframe = np.array(Image.open(os.path.join(bmpdir,framename)))
                print(np.linalg.norm(openframe))
                mymidframe = ndimage.median_filter(openframe,5)
                mpca[tindex,:,:] = mymidframe


print(len(mpca))




# now we are going to take all the differences.
raw_list=[]
norm_list = []
normalized_list = []


mpca = np.squeeze(mpca)
framenos = np.squeeze(np.array(framenos))
ts = np.squeeze(np.array(ts))

# frmtimes = np.squeeze(np.array(frmtimes))
# frmf3 = np.squeeze(np.array(frmf3))
# frmf2 = np.squeeze(np.array(frmf2))
# frmf1 = np.squeeze(np.array(frmf1))


print(mpca)


keep_indices = np.where(~np.isnan(mpca).any(axis=(1,2)))[0]
print(keep_indices)

kept_mpca = mpca[keep_indices]
kept_ts = ts[keep_indices]

# try mask
# masking is being done on mpca

testframe = kept_mpca[5] # 
rds = 65 # test radius
cx = s/2

success = 0
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [cx,h], radius = rds)
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()

    resp = input('Is the mask centered? (y, tr, tl, sug)')
    if resp == 'y' : 
        success = 1
    elif resp == 'tr':
        cx = cx - 5
    elif resp == 'tl':
        cx = cx + 5
    elif resp == 'sug':
        cx = input('Enter an x-coordinate. Currently displayed is ' + str(cx))
        cx = int(cx)
plt.close()

success = 0
while success == 0:
    h,w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [cx,h], radius = rds)
    masked_img = masked_img.copy() #testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()

    resp = input('Is the mask too big, too small, or juuuuust right? (tb, ts, gdx, sug)')
    if resp == 'tb':
        rds = rds - 3
    elif resp == 'ts':
        rds = rds + 7
    elif resp == 'gdx':
        success = 1
    elif resp == 'sug':
        rds = input('Enter a radius. Currently displayed is ' + str(rds))
        rds = int(rds)
    print(rds)
plt.close()



for i in range(0,(kept_mpca.shape[0])):
    maskframe = kept_mpca[i]
    mask = create_circular_mask(h, w, center = [cx,h], radius = rds)
    masked_img = maskframe.copy()
    masked_img[mask] = 0
    kept_mpca[i] = masked_img

testframe = kept_mpca[5]

# semi-circular mask

cx = s/2
radius = 180
success = 0
while success == 0 :
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h,w, center = [cx,q], radius = radius, shape = 'semi')
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()
    resp = input('Is the mask centered? (y, tl, tr)')

    if resp == 'y' : 
        success = 1
    elif resp == 'tr':
        cx = cx - 5
    elif resp == 'tl':
        cx = cx + 5
    elif resp == 'sug':
        cx = input('Enter an x-coordinate. Currently displayed is ' + str(cx))
        cx = int(cx)
plt.close()

success = 0
while success == 0:    
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [cx,q], radius = radius, shape = 'semi')
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
print(radius)
testframe = masked_img

success = 0
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [cx, q], radius = radius, shape = 'semi')
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


for i in range(0,(kept_mpca.shape[0])):
    maskframe = kept_mpca[i]
    mask = create_circular_mask(h, w, center = [cx,q], radius = radius, shape = 'semi')
    masked_img = maskframe.copy()
    masked_img[mask] = 0
    if i == 1:
        plt.imshow(masked_img, cmap = "Greys_r")
        plt.show()
    kept_mpca[i] = masked_img
plt.close()


avgfrm = np.mean(kept_mpca,axis=0) # mean frame
subavg = np.linalg.norm(np.mean(kept_mpca,axis=0)) # average brightness

# mpca_seq = list(range((len(kept_mpca)-1)))


for i in np.arange(len(kept_mpca)):
#    print(len(kept_mpca))
    rawdiff = kept_mpca[i,:,:] - avgfrm
    raw_list.append(rawdiff)
    normdiff = np.linalg.norm(rawdiff)
    norm_list.append(normdiff)
    normalizeddiff = normdiff/subavg
    normalized_list.append(normalizeddiff)

print(len(kept_mpca))
print(len(normalized_list))

norm_avg = np.mean(norm_list)
normalized_avg = np.mean(normalized_list)
sd_norm = np.std(norm_list)
sd_normalized = np.std(normalized_list)


#####       # now we are going to take all the differences using midpoint frames.
#####       mid_raw_list=[]
#####       mid_norm_list = []
#####       mid_normalized_list = []
#####       
#####       
#####       mpca = np.squeeze(mpca)
#####       #mframenos = np.squeeze(np.array(mframenos))
#####       #ts = np.squeeze(np.array(mts))
#####       #frmtimes = np.squeeze(np.array(frmtimes))
#####       #frmf3 = np.squeeze(np.array(frmf3))
#####       #frmf2 = np.squeeze(np.array(frmf2))
#####       #frmf1 = np.squeeze(np.array(frmf1))
#####       # deal with acoustics later
#####       
#####       
#####       print(mpca)
#####       
#####       mid_keep_indices = np.where(~np.isnan(mpca).any(axis=(1,2)))[0]
#####       print(mid_keep_indices)
#####       
#####       kept_mpca = mpca[mid_keep_indices]
#####       #kept_framenos = framenos[keep_indices]
#####       #kept_ts = ts[keep_indices]
#####       #kept_frmtimes = frmtimes[keep_indices] # when the frame occurs from sync.txt
#####       #kept_stimulus = np.array(stimulus,str)[keep_indices] # stimulus word
#####       #kept_f3 = frmf3[keep_indices]
#####       #kept_f2 = frmf2[keep_indices]
#####       #kept_f1 = frmf1[keep_indices]
#####       
#####       mid_subavg = np.linalg.norm(np.mean(kept_mpca,axis=0))# find average along axis 0
#####       midavgfrm = np.mean(kept_mpca,axis=0)
#####       
#####       mpca_seq = list(range((len(kept_mpca)-1)))
#####       
#####       
#####       
#####       for i in np.arange(len(kept_mpca)):
#####       #    print(len(kept_mpca))
#####           rawdiff = kept_mpca[i,:,:] - avgfrm
#####           mid_raw_list.append(rawdiff)
#####           normdiff = np.linalg.norm(rawdiff)
#####           mid_norm_list.append(normdiff)
#####           normalizeddiff = normdiff/subavg
#####           mid_normalized_list.append(normalizeddiff)
#####       
#####       
#####       # for double in combinations(mpca_seq,2):
#####       #    minuend = double[0]
#####       #    subtrahend = double[1]    
#####       #    rawdiff = kept_mpca[minuend,:,:]-kept_mpca[subtrahend,:,:]
#####       #    normdiff = np.linalg.norm(rawdiff)
#####       #    normalizeddiff = normdiff/subavg
#####       #    mid_raw_list.append(rawdiff)
#####       #    mid_norm_list.append(normdiff)
#####       #    mid_normalized_list.append(normalizeddiff)
#####       #    #print(double)
#####       mid_norm_avg = np.mean(mid_norm_list)
#####       mid_normalized_avg = np.mean(mid_normalized_list)
#####       mid_sd_norm = np.std(mid_norm_list)
#####       mid_sd_normalized = np.std(mid_normalized_list)



outfiname = TARG + '_mask_diffs_frombase_2020.txt'
outdiffs = os.path.join(subbmpdir, outfiname)
# middiffs = os.path.join(subbmpdir,'r_middiffs.txt')

od = open(outdiffs, 'w')
# od.write('\t'.join(['subject','norm_avg','normalized_avg','sd_norm','sd_normalized','mid_norm_avg','mid_normalized_avg','mid_sd_norm','mid_normalized'])+'\n')
od.write('\t'.join(['subject','norm_avg','normalized_avg','sd_norm','sd_normalized'])+'\n')
od.write('\t'.join([args.subject,str(norm_avg),str(normalized_avg),str(sd_norm),str(sd_normalized)])+'\n') 
# ,str(mid_norm_avg),str(mid_normalized_avg),str(mid_sd_norm),str(mid_sd_normalized)])+'\n')
#f.write(str(covariance))
od.close()

datafiname = TARG+'_data.txt'
datafi = os.path.join(subbmpdir, datafiname)
data_headers = ["timestamp","framenum","normdiff","normalizeddiff"]#"midnormdiff","midnormalizeddiff"] 


h = np.row_stack((data_headers,np.column_stack((kept_ts,framenos,norm_list,normalized_list))))#mid_norm_list,mid_normalized_list))))
np.savetxt(datafi, h, fmt = "%s", delimiter = ",") 

maskfi = TARG + 'maskparams.txt'
maskparams = os.path.join(subbmpdir, maskfi)
mp = open(maskparams, 'w')
mp.write('\t'.join(['posq','low_radius','high_radius'])+'\n')
mp.write('\t'.join([str(q), str(rds), str(radius)]))
mp.close()


# do this to print midpoint frame and some diff frames

pic=avgfrm
mag = np.max(pic)-np.min(pic)
pic = (pic-np.min(pic))/mag*255
printy = misc.imsave(os.path.join(subbmpdir,'AverageS.png'),pic)


# raw_list

# image_shape = (q,s)
for r in range(0, len(raw_list):
   pic = raw_list[r]
   mag = np.max(pic) - np.min(pic)
   pic = (pic-np.min(pic))/mag*255
   pcn = misc.imsave(os.path.join(subbmpdir,'dfb'+str(r+1)+'.png'),pic)





print(framenos)

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
    print(kept_f1)
    print(kept_f2)
    print(kept_f3)
    
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





