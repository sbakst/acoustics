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



# this could be made a separate file I read in

WORDS = ['rah', 'rome', 'Rome', 'ream', 'bar', 'bore', 'beer', 'RAH', 'ROME', 'REAM','BAR','BORE','BEER']


parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Experiment directory containing all subjects' bmps and jpgs")
parser.add_argument("otherdir", help="Experiment directory containing all subjects' tgs and sync files and wavs etc")
parser.add_argument("subject", help="subjnumber")

# parser.add_argument("-n", "n_components", help="Number of principal components to output", action="store_true")
parser.add_argument("-v", "--visualize", help="Produce plots of PC loadings on fan",action="store_true")
parser.add_argument("-p", "--pca", help="run pca",action="store_true")

args = parser.parse_args()


# PARENTDIR = sys.argv[1]
# DESTDIR = sys.argv[2]
# SUB = sys.argv[3]

subbmpdir = os.path.join(args.directory,args.subject)
utildir = os.path.join(args.otherdir,args.subject)

# required vectors:
# vowel?
# frame (array)
framenos = [] # frame number
frmtimes = [] # when the frame occurs from sync.txt
stimulus = [] # stimulus word
ts = []
frmf3 = []
frmf2= []
frmf1 = []
# First, for each directory of bmps, read in the stim.txt.

rpca = None



for dirs, times, files in os.walk(subbmpdir):
# When you see this a year from now, you can remember sitting in the Illini union. ;)
    for tindex, timestamp in enumerate(times):
# append some kind of placeholder for when you skip a file
        ts.append(timestamp)
        subjrframelist = []
        utt = os.path.join(utildir,str(timestamp))
        syncfile = os.path.join(utt,(timestamp+'.bpr.sync.txt'))
        print(syncfile)
        if not os.path.isfile(syncfile):
            print("can't find syncfile")
#        print(syncfile)
#        if not syncfile:
            continue
        stimfile = os.path.join(utt, 'stim.txt')
        stim = open(stimfile)
        stimtext = stim.read()
        print (stimtext)
        stimulus.append(str(stimtext))
        if any(substring in stimtext for substring in WORDS): 
            tg = os.path.join(utt,(timestamp+'.TextGrid')) # may need to change to just TextGrid depending on when data is from
            if not tg:
                tg = os.path.join(utt,(timestamp+'.bpr.ch1.TextGrid'))
            wav = os.path.join(utt,(timestamp+'.bpr.ch1.wav')) # may need to change to bpr.ch1.wav depending on when data is from
            if not wav:
                wav = os.path.join(utt,(timestamp+'.wav'))
            print(wav)
            # syncfile = os.path.join(utt,(timestamp+'.bpr.sync.txt'))
       #     if not syncfile:
       #         continue
            subjrframelist = re.compile('.*\.jpg')
            regex = re.compile('(pc[0-9])|(mean).jpg')
            bmpdir = os.path.join(subbmpdir, timestamp)
            imlist = [i for i in os.listdir(bmpdir) if (subjrframelist.search(i) and not regex.search(i))]
            print(imlist)
            try:
                im = np.array(Image.open(os.path.join(bmpdir,(imlist[0])))) #open one image to get the size
            except IndexError as e:
                myrframe = 'NA'
                continue
            if len(imlist) > 3:
#                framenos.append('NA')
#                frmtimes.append('NA')
#                continue 
                q,s = im.shape[0:2] #get the size of the images
                print(q)
                print(s)
                imnbr = len(imlist) #get the number of images
            if rpca is None:
                rpca = np.empty([len(os.listdir(bmpdir))]+list(im.shape[0:2])) * np.nan
            difflist = []
            for i in range(imnbr):
                print(i)
                print(imlist[i])
                print(os.path.join(bmpdir,imlist[i]))
  #              call(["convert", os.path.join(utt,imlist[i]), "-set", "colorspace", "RGB", os.path.join(utt,imlist[i])])
  #              call(["convert", os.path.join(utt,imlist[i]), "-set", "colorspace", "Gray", os.path.join(utt,imlist[i])])
                if i > 0:
                    diffmatrix = np.array(Image.open(os.path.join(bmpdir,imlist[i])))-np.array(Image.open(os.path.join(bmpdir,imlist[i-1])))
                    difflist.append(np.linalg.norm(diffmatrix))
         	    # identify frame where the difference in tongue shape is lowest.
            print(difflist)
            min_val, min_idx = min((val,idx) for (idx, val) in enumerate(difflist))
            print(min_val)
            print(min_idx)
	    # add frame to array
            frame_n = imlist[min_idx]
            frame_number = os.path.splitext(os.path.splitext(frame_n)[0])[1][1:] # actual frame number, not INDEX
            framenos.append(frame_number)
            openframe = np.array(Image.open(os.path.join(bmpdir,frame_n)))
            myrframe = ndimage.median_filter(openframe,5)
            print(myrframe)
            rpca[tindex,:,:] = myrframe
            print(syncfile)

		# do the acoustics: find timepoint (time) in file where that frame occurs using the syncfile
           # if int(args.subject) > 120: # change in syncfile format
            #    sm = audiolabel.LabelManager(from_file=syncfile, from_type='table',sep = '\t',fields_in_head=True)# False, fields='t1,frameidx')
           # else:
####            sm = audiolabel.LabelManager(from_file=syncfile, from_type='table',sep = '\t',fields_in_head=False, fields='t1,frameidx')
####            print(sm)
####            for v,m in sm.tier('frameidx').search(frame_number, return_match=True):
####                print(v)
####                #if int(args.subject) > 120: 
####                frmtime = v.t1
####                print(frmtime)
####                frmtimes.append(frmtime)
####		# oh, christ. Perform rformant and get F1/F2/F3 at frmtime. But let's stop here for now.
####                fbfile = (os.path.splitext(wav)[0]+'.fb') # make fbfile name that rformant will write to
####
####                try:
####                    formant_proc = subprocess.check_call(["rformant", wav])#, stdin=rdc_proc.stdout) # also remove 20
####                except subprocess.CalledProcessError as e:
####                    print(e)
####                    print(e.stdout)
####                ppl_proc = subprocess.Popen(
####                    ["pplain", fbfile],
####                    stdout=subprocess.PIPE)
####                print(fbfile)   
####				
####                lm = audiolabel.LabelManager(
####                   from_file=ppl_proc.stdout,   # read directly from pplain output
####                   from_type='table',
####                   sep=" ",
####                   fields_in_head=False,
####                   fields="f1,f2,f3,f4",
####                   t1_col=None,                 # esps output doesn't have a t1 column
####                   t1_start=0.0,          # autocreate t1 starting with this value and
####                   t1_step=0.01)             # increase by this step
####                f3 = lm.tier('f3')
####                f2 = lm.tier('f2')
####                f1 = lm.tier('f1')
####                framef3 = (f3.label_at(frmtime)).text
####	#            if float(framef3) > 2300:
####	#                framef3 = "NA"
####                framef2 = (f2.label_at(frmtime)).text
####                framef1 = (f1.label_at(frmtime)).text
####		#        row_out = '\t'.join([str(i),str(a.timestamp), str(val), str(meas_t1), framef2, framef3, str(indx), str(len(mindiffs))])
####		#        out.write(row_out+'\n')            
####                frmf3.append(framef3)
####                frmf2.append(framef2)
####                frmf1.append(framef1)
####
####            else:
####	# ok now we append a placeholder to anything that wouldn't otherwise get defined:
####	# frame no
####	# rpca doesn't need to because we used tindx
####	# frmtimes
####                framenos.append('NULL')
####                frmtimes.append('NULL')
####                frmf3.append('NULL')
####                frmf2.append('NULL')
####                frmf1.append('NULL')
####                            
print(len(rpca))


# now we are going to take all the differences.
subavg = np.linalg.norm(np.mean(rpca,axis=0))# find average along axis 0
raw_list=[]
norm_list = []
normalized_list = []

rpca_seq = list(range((len(rpca)-1)))

for double in combinations(rpca_seq,2):
    minuend = double[0]
    subtrahend = double[1]    
    rawdiff = rpca[minuend,:,:]-rpca[subtrahend,:,:]
    normdiff = np.linalg.norm(rawdiff)
    normalizeddiff = normdiff/subavg
    raw_list.append(rawdiff)
    norm_list.append(normdiff)
    normalized_list.append(normalizeddiff)

norm_avg = np.mean(norm_list)
normalized_avg = np.mean(normalized_list)
sd_norm = np.std(norm_list)
sd_normalized = np.std(normalized_list)
print(rawdiff)
print(normdiff)
print(normalizeddiff)
print(norm_list)


outdiffs = os.path.join(bmpdir,'diffs.txt')

od = open(outdiffs, 'w')
od.write('\t'.join(['subject','norm_avg','normalized_avg','sd_norm','sd_normalized'])+'\n')
od.write('\t'.join([SUB,str(norm_avg),str(normalized_avg),str(sd_norm),str(sd_normalized)])+'\n')
#f.write(str(covariance))
od.close()

if args.pca:
####################################################################################################


# remove any indices for all objects generated above where frames have NaN values (due to skipping or otherwise)
# ooooo she's squeezin'!
    rpca = np.squeeze(rpca)
    framenos = np.squeeze(np.array(framenos))
    ts = np.squeeze(np.array(ts))
    frmtimes = np.squeeze(np.array(frmtimes))
    frmf3 = np.squeeze(np.array(frmf3))
    frmf2 = np.squeeze(np.array(frmf2))
    frmf1 = np.squeeze(np.array(frmf1))
    
    
    print(rpca)
    
    
    keep_indices = np.where(~np.isnan(rpca).any(axis=(1,2)))[0]
    print(keep_indices)
    
    kept_rpca = rpca[keep_indices]
    kept_framenos = framenos[keep_indices]
    kept_ts = ts[keep_indices]
    kept_frmtimes = frmtimes[keep_indices] # when the frame occurs from sync.txt
    kept_stimulus = np.array(stimulus,str)[keep_indices] # stimulus word
    kept_f3 = frmf3[keep_indices]
    kept_f2 = frmf2[keep_indices]
    kept_f1 = frmf1[keep_indices]



    n_components = 5
    
    pca = PCA(n_components = n_components)
    
    print(kept_rpca.shape[0]) # 48ish?
    print(kept_rpca.shape[1]) # q?
    print(kept_rpca.shape[2]) # s?
    
    
    frames_reshaped = kept_rpca.reshape([kept_rpca.shape[0], kept_rpca.shape[1]*kept_rpca.shape[2]])
    
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
    
    outfile = os.path.join(bmpdir, 'pca.csv')
    d = np.row_stack((headers,np.column_stack((kept_ts,kept_framenos,kept_f1,kept_f2,kept_f3,analysis))))
    np.savetxt(outfile, d, fmt="%s", delimiter =',')
    
    expl = os.path.join(bmpdir, 'explain.txt')
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
            pcn = misc.imsave(os.path.join(bmpdir,'pc'+str(n+1)+'.png'),pic)
    


#r_mean = np.kept_rpca.mean(0)
#for i in np.arange(length(kept_rpca[0])):
#   std=(x-mu)


#standev = 
#f = open(





