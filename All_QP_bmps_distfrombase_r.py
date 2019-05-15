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
frmtimes = [] # when the frame occurs from sync.txt
stimulus = [] # stimulus word
ts = []
frmf3 = []
frmf2= []
frmf1 = []
# First, for each directory of bmps, read in the stim.txt.

rpca = None
mrpca = None


for dirs, times, files in os.walk(subbmpdir):
# When you see this a year from now, you can remember sitting in the Illini union. ;)
    for tindex, timestamp in enumerate(times):
# append some kind of placeholder for when you skip a file
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
        stimulus.append(str(stimtext))
        if stimmy in WORDS:
#            print('yay!')
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
                    rt1 = lab.t1
                    rt2 = lab.t2
                    rt_frmtime = ((rt1 + rt2)/2)
            subjrframelist = re.compile('.*\.jpg')
            regex = re.compile('(pc[0-9])|(mean).jpg')
            stimex = re.compile(stimtext) 
            bmpdir = os.path.join(subbmpdir, timestamp)
            imlist = [i for i in os.listdir(bmpdir) if (subjrframelist.search(i) and not regex.search(i) and not stimex.search(i))]
            # print(imlist)
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
                #print(q)
                #print(s)
                imnbr = len(imlist) #get the number of images
                # print(imnbr)
                if rpca is None:
                    rpca = np.empty([len(os.listdir(subbmpdir))]+list(im.shape[0:2])) * np.nan
                if mrpca is None:
                    mrpca = np.empty([len(os.listdir(subbmpdir))]+list(im.shape[0:2])) * np.nan
                difflist = []
                for i in range(imnbr):
                    #print(i)
                    # print(imlist[i])
                    # print(os.path.join(bmpdir,imlist[i]))
                    bmpnorm = np.linalg.norm(np.array(Image.open(os.path.join(bmpdir, imlist[i]))))
                    if bmpnorm > 25000:
                        continue
  #                  call(["convert", os.path.join(utt,imlist[i]), "-set", "colorspace", "RGB", os.path.join(utt,imlist[i])])
  #                  call(["convert", os.path.join(utt,imlist[i]), "-set", "colorspace", "Gray", os.path.join(utt,imlist[i])])
                    if i > 0:
                        diffmatrix = np.array(Image.open(os.path.join(bmpdir,imlist[i])))-np.array(Image.open(os.path.join(bmpdir,imlist[i-1])))
                        if (np.linalg.norm(diffmatrix)) > 0:
                            difflist.append(np.linalg.norm(diffmatrix))
                        else:
                            difflist.append('NA')
         	        # identify frame where the difference in tongue shape is lowest.
                print(difflist)
                min_val, min_idx = min((val,idx) for (idx, val) in enumerate(difflist))
                print(min_val)
                print(min_idx)
	               # add frame to array
                frame_n = imlist[min_idx]
                print(frame_n)
     
                frame_number = os.path.splitext(os.path.splitext(frame_n)[0])[1][1:] # actual frame number, not INDEX
                print(frame_number)
                framenos.append(frame_number)
                openframe = np.array(Image.open(os.path.join(bmpdir,frame_n)))
                myrframe = ndimage.median_filter(openframe,5)
                print(myrframe)
                rpca[tindex,:,:] = myrframe
                print(syncfile)
                    
                try:
                    sm = audiolabel.LabelManager(from_file = syncfile, from_type = 'table', sep = '\t', fields_in_head = False, fields = 't1,frameidx')
                    rmidfrm = int(sm.tier('frameidx').label_at(rt_frmtime).text)
                except ValueError:
                    syncfile = os.path.join(utt, (timestamp + '.bpr.sync.TextGrid'))
                    # print(syncfile)
                    sm = audiolabel.LabelManager(from_file = syncfile, from_type='praat')#,sep = '\t', fields_in_head = True, fields = 'seconds,pulse_idx,raw_data_idx')
#                    print(sm)
                    rmidfrm = int(sm.tier('pulse_idx').label_at(rt_frmtime).text)
                framename = timestamp + '.' + str(rmidfrm) + '.jpg'
                if not os.path.isfile(os.path.join(bmpdir,framename)):
                    framename = timestamp + '.' + str(rmidfrm - 1) + '.jpg'
                openframe = np.array(Image.open(os.path.join(bmpdir,framename)))
                myrframe = ndimage.median_filter(openframe,5)
                mrpca[tindex,:,:] = myrframe


		# do the acoustics: find timepoint (time) in file where that frame occurs using the syncfile
#            if int(args.subject) > 120: # change in syncfile format
#                sm = audiolabel.LabelManager(from_file=syncfile, from_type='table',sep = '\t',fields_in_head=True)# False, fields='t1,frameidx')
#            else:
#                sm = audiolabel.LabelManager(from_file=syncfile, from_type='table',sep = '\t',fields_in_head=False, fields='t1,frameidx')
#            print(sm)
# uncomment these for audio            for v,m in sm.tier('frameidx').search(frame_number, return_match=True):
# uncomment these for audio                print(v)
# uncomment these for audio                #if int(args.subject) > 120: 
# uncomment these for audio                frmtime = v.t1
# uncomment these for audio                print(frmtime)
# uncomment these for audio                frmtimes.append(frmtime)
# uncomment these for audio		# oh, christ. Perform rformant and get F1/F2/F3 at frmtime. But let's stop here for now.
# uncomment these for audio                fbfile = (os.path.splitext(wav)[0]+'.fb') # make fbfile name that rformant will write to
# uncomment these for audio
# uncomment these for audio                try:
# uncomment these for audio                    formant_proc = subprocess.check_call(["rformant", wav])#, stdin=rdc_proc.stdout) # also remove 20
# uncomment these for audio                except subprocess.CalledProcessError as e:
# uncomment these for audio                    print(e)
# uncomment these for audio                    print(e.stdout)
# uncomment these for audio                ppl_proc = subprocess.Popen(
# uncomment these for audio                    ["pplain", fbfile],
# uncomment these for audio                    stdout=subprocess.PIPE)
# uncomment these for audio                print(fbfile)   
# uncomment these for audio				
# uncomment these for audio                lm = audiolabel.LabelManager(
# uncomment these for audio                   from_file=ppl_proc.stdout,   # read directly from pplain output
# uncomment these for audio                   from_type='table',
# uncomment these for audio                   sep=" ",
# uncomment these for audio                   fields_in_head=False,
# uncomment these for audio                   fields="f1,f2,f3,f4",
# uncomment these for audio                   t1_col=None,                 # esps output doesn't have a t1 column
# uncomment these for audio                   t1_start=0.0,          # autocreate t1 starting with this value and
# uncomment these for audio                   t1_step=0.01)             # increase by this step
# uncomment these for audio                f3 = lm.tier('f3')
# uncomment these for audio                f2 = lm.tier('f2')
# uncomment these for audio                f1 = lm.tier('f1')
# uncomment these for audio                framef3 = (f3.label_at(frmtime)).text
# uncomment these for audio	#            if float(framef3) > 2300:
# uncomment these for audio	#                framef3 = "NA"
# uncomment these for audio                framef2 = (f2.label_at(frmtime)).text
# uncomment these for audio                framef1 = (f1.label_at(frmtime)).text
# uncomment these for audio		#        row_out = '\t'.join([str(i),str(a.timestamp), str(val), str(meas_t1), framef2, framef3, str(indx), str(len(mindiffs))])
# uncomment these for audio		#        out.write(row_out+'\n')            
# uncomment these for audio                frmf3.append(framef3)
# uncomment these for audio                frmf2.append(framef2)
# uncomment these for audio                frmf1.append(framef1)
# uncomment these for audio
# uncomment these for audio            else:
# uncomment these for audio	# ok now we append a placeholder to anything that wouldn't otherwise get defined:
# uncomment these for audio	# frame no
# uncomment these for audio	# rpca doesn't need to because we used tindx
# uncomment these for audio	# frmtimes
# uncomment these for audio                framenos.append('NULL')
# uncomment these for audio                frmtimes.append('NULL')
# uncomment these for audio                frmf3.append('NULL')
# uncomment these for audio                frmf2.append('NULL')
# uncomment these for audio                frmf1.append('NULL')
# uncomment these for audio                            
print(len(rpca))




# now we are going to take all the differences.
raw_list=[]
norm_list = []
normalized_list = []


rpca = np.squeeze(rpca)
framenos = np.squeeze(np.array(framenos))
ts = np.squeeze(np.array(ts))

# frmtimes = np.squeeze(np.array(frmtimes))
# frmf3 = np.squeeze(np.array(frmf3))
# frmf2 = np.squeeze(np.array(frmf2))
# frmf1 = np.squeeze(np.array(frmf1))


print(rpca)


keep_indices = np.where(~np.isnan(rpca).any(axis=(1,2)))[0]
print(keep_indices)

kept_rpca = rpca[keep_indices]
kept_ts = ts[keep_indices]

# try mask

testframe = kept_rpca[5] # 
rds = 65 # test radius

success = 0
while success == 0:
    h, w = testframe.shape[:2]
    mask = create_circular_mask(h, w, center = [w/2,h], radius = rds)
    masked_img = testframe.copy()
    masked_img[mask] = 0
    plt.imshow(masked_img, cmap = "Greys_r")
    plt.show()

    resp = input('Is the mask too big, too small, or juuuuust right? (tb, ts, gdx)')
    if resp == 'tb':
        rds = rds - 3
    elif resp == 'ts':
        rds = rds + 7
    elif resp == 'gdx':
        success = 1
    print(rds)
    plt.close()



for i in range(0,(kept_rpca.shape[0])):
    maskframe = kept_rpca[i]
    mask = create_circular_mask(h, w, center = [w/2,h], radius = rds)
    masked_img = maskframe.copy()
    masked_img[mask] = 0
    kept_rpca[i] = masked_img

testframe = kept_rpca[5]

# semi-circular mask
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


for i in range(0,(kept_rpca.shape[0])):
    maskframe = kept_rpca[i]
    mask = create_circular_mask(h, w, center = [s/2,q], radius = radius, shape = 'semi')
    masked_img = maskframe.copy()
    masked_img[mask] = 0
    if i == 1:
        plt.imshow(masked_img, cmap = "Greys_r")
        plt.show()
    kept_rpca[i] = masked_img
plt.close()
# add mask
# for i in range(0,(kept_rpca.shape[0])):
#    (kept_rpca)[i][0:50] = 0

# d = (kept_rpca)[4]
# mag = np.max(d) - np.min(d)
# d = (d-np.min(d))/mag*255
#testframe = np.flipud(e.acquisitions[0].image_converter.as_bmp(d)) # converter from any frame will work; here we use the first
# testframe = Image.fromarray(d,'RGB')
# testframe.save('lol.jpg')

#plt.title("Did I turn the bottom row to 0?") #.format?
#plt.imshow(testframe, cmap="Greys_r") 
#savepath = "a_test_frame_lol.pdf" #.format?
#plt.savefig(savepath)


# kept_framenos = framenos[keep_indices]
# kept_ts = ts[keep_indices]
# kept_frmtimes = frmtimes[keep_indices] # when the frame occurs from sync.txt
# kept_stimulus = np.array(stimulus,str)[keep_indices] # stimulus word
# kept_f3 = frmf3[keep_indices]
# kept_f2 = frmf2[keep_indices]
# kept_f1 = frmf1[keep_indices]

avgfrm = np.mean(kept_rpca,axis=0)
subavg = np.linalg.norm(np.mean(kept_rpca,axis=0))# find average along axis 0

rpca_seq = list(range((len(kept_rpca)-1)))


for i in np.arange(len(kept_rpca)):
#    print(len(kept_rpca))
    rawdiff = kept_rpca[i,:,:] - avgfrm
    raw_list.append(rawdiff)
    normdiff = np.linalg.norm(rawdiff)
    norm_list.append(normdiff)
    normalizeddiff = normdiff/subavg
    normalized_list.append(normalizeddiff)


# for double in combinations(rpca_seq,2):
#     minuend = double[0]
#     subtrahend = double[1]    
#     rawdiff = kept_rpca[minuend,:,:]-kept_rpca[subtrahend,:,:]
#     normdiff = np.linalg.norm(rawdiff)
#     normalizeddiff = normdiff/subavg
#     raw_list.append(rawdiff)
#     norm_list.append(normdiff)
#     normalized_list.append(normalizeddiff)
#     #print(double)



norm_avg = np.mean(norm_list)
normalized_avg = np.mean(normalized_list)
sd_norm = np.std(norm_list)
sd_normalized = np.std(normalized_list)
#print(rawdiff)
#print(normdiff)
#print(normalizeddiff)
#print(norm_list)


# now we are going to take all the differences using midpoint frames.
mid_raw_list=[]
mid_norm_list = []
mid_normalized_list = []


mrpca = np.squeeze(mrpca)
#mframenos = np.squeeze(np.array(mframenos))
#ts = np.squeeze(np.array(mts))
#frmtimes = np.squeeze(np.array(frmtimes))
#frmf3 = np.squeeze(np.array(frmf3))
#frmf2 = np.squeeze(np.array(frmf2))
#frmf1 = np.squeeze(np.array(frmf1))
# deal with acoustics later


print(mrpca)

mid_keep_indices = np.where(~np.isnan(mrpca).any(axis=(1,2)))[0]
print(mid_keep_indices)

kept_mrpca = mrpca[mid_keep_indices]
#kept_framenos = framenos[keep_indices]
#kept_ts = ts[keep_indices]
#kept_frmtimes = frmtimes[keep_indices] # when the frame occurs from sync.txt
#kept_stimulus = np.array(stimulus,str)[keep_indices] # stimulus word
#kept_f3 = frmf3[keep_indices]
#kept_f2 = frmf2[keep_indices]
#kept_f1 = frmf1[keep_indices]

mid_subavg = np.linalg.norm(np.mean(kept_mrpca,axis=0))# find average along axis 0
midavgfrm = np.mean(kept_mrpca,axis=0)

mrpca_seq = list(range((len(kept_mrpca)-1)))



for i in np.arange(len(kept_mrpca)):
#    print(len(kept_rpca))
    rawdiff = kept_mrpca[i,:,:] - avgfrm
    mid_raw_list.append(rawdiff)
    normdiff = np.linalg.norm(rawdiff)
    mid_norm_list.append(normdiff)
    normalizeddiff = normdiff/subavg
    mid_normalized_list.append(normalizeddiff)


# for double in combinations(mrpca_seq,2):
#    minuend = double[0]
#    subtrahend = double[1]    
#    rawdiff = kept_mrpca[minuend,:,:]-kept_mrpca[subtrahend,:,:]
#    normdiff = np.linalg.norm(rawdiff)
#    normalizeddiff = normdiff/subavg
#    mid_raw_list.append(rawdiff)
#    mid_norm_list.append(normdiff)
#    mid_normalized_list.append(normalizeddiff)
#    #print(double)
mid_norm_avg = np.mean(mid_norm_list)
mid_normalized_avg = np.mean(mid_normalized_list)
mid_sd_norm = np.std(mid_norm_list)
mid_sd_normalized = np.std(mid_normalized_list)



outfiname = TARG + '_mask_diffs_frombase.txt'
outdiffs = os.path.join(subbmpdir, outfiname)
# middiffs = os.path.join(subbmpdir,'r_middiffs.txt')

od = open(outdiffs, 'w')
od.write('\t'.join(['subject','norm_avg','normalized_avg','sd_norm','sd_normalized','mid_norm_avg','mid_normalized_avg','mid_sd_norm','mid_normalized'])+'\n')
od.write('\t'.join([args.subject,str(norm_avg),str(normalized_avg),str(sd_norm),str(sd_normalized),str(mid_norm_avg),str(mid_normalized_avg),str(mid_sd_norm),str(mid_sd_normalized)])+'\n')
#f.write(str(covariance))
od.close()

datafiname = TARG+'_data.txt'
datafi = os.path.join(subbmpdir, datafiname)
data_headers = ["timestamp","framenum","normdiff","midnormdiff"] 


h = np.row_stack((data_headers,np.column_stack((kept_ts,framenos,norm_list,mid_norm_list))))
np.savetxt(datafi, h, fmt = "%s", delimiter = ",") 

maskfi = TARG + 'maskparams.txt'
maskparams = os.path.join(subbmpdir, maskfi)
mp = open(maskparams, 'w')
mp.write('\t'.join(['posq','low_radius','high_radius'])+'\n')
mp.write('\t'.join([str(q), str(rds), str(radius)]))
mp.close()




print(framenos)

if args.pca:
####################################################################################################


# remove any indices for all objects generated above where frames have NaN values (due to skipping or otherwise)
# ooooo she's squeezin'!


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
    


#r_mean = np.kept_rpca.mean(0)
#for i in np.arange(length(kept_rpca[0])):
#   std=(x-mu)


#standev = 
#f = open(





