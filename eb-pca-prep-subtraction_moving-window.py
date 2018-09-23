"""
TODO docstring
"""

from __future__ import absolute_import, division, print_function

import os, sys, glob, re
import struct
import argparse
import audiolabel
from operator import itemgetter
import numpy as np
from scipy import ndimage
import subprocess
import pandas as pd
from hashlib import sha1
from collections import OrderedDict
import csv
from ultratils.rawreader import RawReader
from ultratils.pysonix.scanconvert import Converter
import matplotlib.pyplot as plt

## define frame-getting functions.
#def read_metadata(mfile):
#    '''Read the metadata from an .img.txt file.'''
#    mdict = {}
#    with open(mfile, 'r') as f:
#        k = f.readline().strip().split("\t")
#        v = f.readline().strip().split("\t")
#        for fld,val in zip(k, v):
#            mdict[fld] = int(val)
#    return mdict

def read_echob_metadata(rawfile):
    '''
    Gather information about a .raw file from its .img.txt file. 
    For legacy .raw data without a header; if a header exists,
    use ultratils utilities.
    Inputs: a .raw file, which is assumed to have an .img.txt file
      with the same base name.
    Outputs:
      nscanlines, the number of scan lines ("width" of unconverted img)
      npoints, the number of pixels in each scan line ("height" of img)
      not_junk, the pixel index in each scan line where junk data begins
    '''
    mfile = os.path.splitext(rawfile)[0] + ".img.txt"
    mdict = {}
    with open(mfile, 'r') as mf:
        k = mf.readline().strip().split("\t")
        v = mf.readline().strip().split("\t")
        for fld,val in zip(k, v):
            mdict[fld] = int(val)
    
    nscanlines = mdict['Height']
    npoints = mdict['Pitch']
    junk = npoints - mdict['Width']
    
    return nscanlines, npoints, junk




#
#def get_frame_size(mfile):
#    mdict = read_metadata(mfile)
#    frame_size_in_bytes = mdict['nInBufferLen']
#    return frame_size_in_bytes
#
#def get_frame(fname, framenum, frame_size, med_filter=False):
#    """Requires raw filename, desired frame number, and size of frame as arguments."""
#    framesize = frame_size
#    data_fmt = 'I' * np.int(framesize / 4) 
#    x_start = framesize * framenum
#    frame_dim_1 = 127
#    frame_dim_2 = 255
#    with open(fname, 'rb') as fh:
#        fh.seek(x_start) # define starting point at leftmost X coord.
#        packed_data = fh.read(framesize) # take a frame-sized chunk at this coord.
#        unpacked_data = struct.unpack(data_fmt, packed_data) 
#        data = np.array(unpacked_data)
#        rdata = np.flipud(data[np.arange(frame_dim_1*frame_dim_2)].reshape(frame_dim_1,frame_dim_2).T)
#        if med_filter == True:
#            rdata = ndimage.median_filter(rdata, 10)
#        return(rdata)

def read_stimfile(stimfile):
    with open(stimfile, "r") as stfile:
        stim = stfile.read().rstrip('\n')
    return stim

# a list of targets from dict.local, to be updated as required.
target_list = ['PESU', 'PO', 'AZHA']



#iz_list = ['IZ', 'BIZX', 'SIZ', 'XIZ']

vre = re.compile(
         "^(IY1|IH1|UH1|UW1|OW1|AE1|EY1|R)$" 
    )

# AH1 and AA1 not to be analyzed

# read in arguments 
parser = argparse.ArgumentParser()
parser.add_argument("ultradirectory", help="Experiment directory containing all the ultrasound")
parser.add_argument("acousticsdir", help="Experiment directory containing audio")
parser.add_argument("sub", help="Subject number")
parser.add_argument("stimulus",help="stimulus")
args = parser.parse_args()

# check for appropriate directory
try:
    expdir = args.ultradirectory
except IndexError:
    print("\tDirectory provided doesn't exist")
    ArgumentParser.print_usage
    ArgumentParser.print_help
    sys.exit(2)

# assemble a list of all *.raw files
rawfile_glob_exp = os.path.join(os.path.normpath(expdir),"*","*.raw")

print("Starting...")

# initialize rest of params for output array
data = None # array to hold ultrasound data
frame_dim_1 = 127
frame_dim_2 = 255
recs = [] # metadata store

# distance (in frames) away from intended time point that can be subbed in
threshhold = 3

# output filepaths
logfile = os.path.join(expdir,"frames_log.txt")
discard_folder = os.path.join(expdir,"discards")
frames_out = os.path.join(expdir,"moving_frames.npy")
metadata_out = os.path.join(expdir,"moving_frames_metadata.pickle")

# create the dictionary

tuplist = []
for dirs, subdir, files in os.walk(args.acousticsdir):
    for wav in files:
        if os.path.splitext(wav)[1]=='.wav':
            tup = re.split('_|.wav',wav)
            if tup[1] == 'inbuf':
#                print(tup[3])
                tuplist.append(tuple(tup))
#print(tuplist)
# sort tuples by trial number, index 3

ordtrials = []
ordtrials = sorted(tuplist, key = lambda x: int(x[3]))
#print(ordtrials)
# make list of ultrasounds

ultras = []
for dirs, subdir, files in os.walk(expdir):
    for tp in subdir:
# split up timestamp
        tptup = re.split('-',tp)
        if tptup[0] == '2017':
            ultras.append(tuple(tptup))
#        print(tptup)

# split order by 3rd element, time

ordultras = []
ordultras = sorted(ultras, key = lambda x: x[2])

# make dictionary for looking up wav given timestamp

tamildict = dict(zip(ordultras,ordtrials))
#print(tamildict)
#dictitems = tamildict.items()
#tamtam = (sorted(dictitems, key= lambda x:int(x[0][3])))
#print(tamtam)

####################
with open(logfile,"w") as header:
    header.write("acq"+"\t"+"stim"+"\t"+"phone"+"\t"+"status"+"\t"+"problem"+"\n")

for rf in glob.glob(rawfile_glob_exp):

    # get stim; skip over practice and bolus acqs
    parent = os.path.dirname(rf)
    acq = os.path.split(parent)[1]
    stim = args.stimulus
#    stimfile = os.path.join(parent,"stim.txt")
#    stim = read_stimfile(stimfile)
#    print(stim)
    if stim == "bolus" or stim == "practice":
        continue

    print("Now working on {:}...".format(acq))
    # get base name and get frame size from local .img.txt file
    barename = os.path.splitext(rf)[0]

#    # TODO proof this against breaking if first acq doesn't have img.txt
#    try:
#        myframesize = get_frame_size(barename + '.img.txt')
#        print(myframesize)
#    except ValueError:
#        print("\tUsed previous frame's size, .img.txt not generated")
#        myframesize = 129540
#        pass
    # new Matt and Ron utils
    nscanlines, npoints, junk = read_echob_metadata(rf)
    print(nscanlines, npoints, junk)

    # use dictionary linking wav files to ultrasound to look up wav file
       
    for key,value in tamildict.items():
#        print(key)
        joined_tp = ('-'.join(key))
#        print(joined_tp)
#        print(str(joined_tp))
#        print(str(barename))
        if joined_tp == os.path.basename(barename):
            print(barename)
            print('banana')
            trial = (value[3])
#            print(trial)
            soundname = ('_'.join(value))[:-1]
            wav = os.path.join(args.acousticsdir, (soundname + '.wav'))
            print(wav)
            tg = os.path.join(args.acousticsdir, (soundname + '.TextGrid'))

    # get WAV file and associated TextGrid; create label manager
#    wav = barename + '.ch1.wav'
#    tg = barename + '.ch1.TextGrid'
    pm = audiolabel.LabelManager(from_file=tg, from_type="praat")
#    sync_tg = barename + '.sync.TextGrid'
#    sync_pm = audiolabel.LabelManager(from_file=sync_tg, from_type="praat")

    # If not sync tg, create list of diffs between frames
    # get number of frames in file
    idxfile = barename + '.idx.txt' # file should just be a column of numbers
    with open(idxfile) as infile:
        idxreader = csv.reader(infile)
        d = list(idxreader)
        rows = sum(1 for row in d)
        last_idx = int((d[rows-1])[0])
#        print(last_idx)

    # loop over all target segments present in TG
    for v,m in pm.tier('phone').search(vre, return_match=True):
        targ = pm.tier('word').label_at(v.center).text

        # skip all instances of any other words not relevant
        if targ not in target_list:
            continue
        else:
            phone = v.text

            before = pm.tier('phone').prev(v).text # get preceding C for vowels
#        after = pm.tier('phone').next(v).text # get following V for consonants
# took out after line because probably not necessary?

        # debug line for acq content
#        print("\t","phone",phone,"stim",stim,"targ",targ,"before",before)#,"after",after)

        # get midpoint time and find closest ultrasound frame in sync TG
        # TODO more efficient to duplicate ultratils frame_at approach


### CHECK INDENTATION BELOW

        mid_timepoint = v.center
        t1 = v.t1
#        print(t1)
        t2 = v.t2
        vlen = t2-t1
        midperc = (mid_timepoint-t1)/vlen

#        prebuf = 12
#        postbuf = 4
        prebufs = [15, 14, 13, 12]
        postbufs = [4, 5, 6, 7]

        frameperc = int(last_idx * midperc)
        print(frameperc)

    rlframes = []
    for w in np.arange(0,4):
        frame_start = frameperc - prebufs[w]
        frame_end = frameperc + postbufs[w]
#        print(frame_start)
#        print(frame_end)


        consec_diffs = []
        for f in np.arange(frame_start,frame_end):
    # old method of getting frames 
    #        minuend = get_frame(rf,f+1, myframesize, med_filter = True)
    #        subtrahend = get_frame(rf, f, myframesize, med_filter = True)
            rdr = RawReader(rf, nscanlines=nscanlines, npoints=npoints) # or whatever the appropriate parameters actually are
            minuend = rdr.get_frame(f+1)
            subtrahend = rdr.get_frame(f)
            cdiff = minuend-subtrahend
            cdiffnorm = np.linalg.norm(cdiff)
            consec_diffs.append(cdiffnorm)
#    print(consec_diffs)
        mindiff,mindiff_idx = min((val,idx) for (idx,val) in enumerate(consec_diffs))

#    print(mindiff)
#    print(mindiff_idx)
        rl_frame_idx = frame_start + mindiff_idx
        print(rl_frame_idx)
        rlframes.append(rl_frame_idx)
    # try using the minimum of these?
    mini, minidind = min((val,idx) for (idx,val) in enumerate(rlframes))
    maxi, maxidind = max((val,idx) for (idx,val) in enumerate(rlframes))
    indsdiff = maxi-mini
    rl_frame_idx = mini
    if (indsdiff > 6) or (0 < (maxi-midperc) < 3): #(midperc-mini)):
        print('check' + trial)

        class Header(object):
            def __init__(self):
                pass

        class Probe(object):
            def __init__(self):
                pass

        header = Header()
        print(header)
# v    isualize components on (approximately) converted fan, if desired
        header.w = 127          # input image width
        header.h = 1020 #255          # input image height
        header.sf = 4000000     # for ultrasonix this is the vec-freq value
        probe = Probe()
        probe.radius = 10000    # based on '10' in transducer model number
        probe.numElements = 128 # based on '128' in transducer model number
        probe.pitch = 185 #205  # guess based on Ultrasonix C9-5/10 transducer
        c = Converter(header, probe)

        image_shape = (1020,127)#(255,127)
#        rdr = RawReader(rf, nscanlines=nscanlines, npoints = npoints)
        print(rf)
        for f in np.arange((midperc-16),(midperc+7)):
            d = rdr.get_frame(f).reshape(image_shape)
            mag = np.max(d) - np.min(d)
            d = (d-np.min(d))/mag*255
            pcn = np.flipud(c.as_bmp(np.flipud(d)))
            #pcn = c.as_bmp(d)
            plt.title("Frame{:}, Subj. {:}".format((f+1),subject))
            plt.imshow(pcn, cmap="Greys_r")
            file_ending = "subj{:}-{:}.png".format(subject, (f+1))
            ultradir = os.path.join(os.path.basename(barename),file_ending)
            savepath = os.path.join(expdir,ultradir)
            print(savepath)
        frameno = input('select a frame number for ' + sub + 'trial' + trial)             
        rl_frame_idx = frameno-1

        # get frame, and check for NaN frames
    change = 0
    discard_acq = False
    while True:
#        pre_rawdata = get_frame(rf,rl_frame_idx,myframesize,med_filter=False)
        pre_rawdata = rdr.get_frame(rl_frame_idx) 
#        print('here is prerawdata')
#        print(pre_rawdata)
#        print('Nate is the best')
        if pre_rawdata is None:
#            print('pre_rawdata is None')
            mid_frame_num -= 1
            change += 1
            if change > threshhold:
                print('change is greater than threshhold')
                with open(logfile, "a") as log:
                    log.write(acq+"\t"+stim+"\t"+phone+"\t"+"discarded"+"\t"+"passed threshhold")
                    print("Frame change threshhold passed; acq {} discarded".format(acq))
                    discard_acq = True
                break
            else:
                pass
        else:
            if change > 0:
                print('change is greater than zero')
                with open(logfile, "a") as log:
                    log.write(acq+"\t"+stim+"\t"+phone+"\t"+"changed by {:}".format(change)+"\t"+"N/A")
                print("Changed target in {:} by".format(acq), change, "frames")
            break
        # discard the acquisition if needed
        print('more poop')
    if discard_acq:
        print('poop')
        shutil.copytree(parent, os.path.join(discard_folder,acq))
        shutil.rmtree(parent)
        continue 
        # preprocessing of images
#    print('now we do the rawdata')
    rawdata = pre_rawdata.astype(np.uint8)
#    print('here is a rawdata')
#    print(rawdata)
        # generate metadata object for the current acquisition
    recs.append(
        OrderedDict([
            ('timestamp', acq),
            ('trial', trial),
            ('time', v.center),
            ('pulseidx', int(rl_frame_idx)),
            ('rawdataidx', int(rl_frame_idx)),
            ('width', frame_dim_1),
            ('height', frame_dim_2),
            ('phone', phone),
            ('stim', stim),
            ('targ', targ),
            ('before', before),
            ('sha1', sha1(rawdata.ravel()).hexdigest()),
            ('sha1_dtype', rawdata.dtype)
        ])
    )
    print('recs')
#    print(recs)
#    print(trial)
        # add frame to frames list
    if data is None:
        data = np.expand_dims(rawdata, axis=0)
    else:
        data = np.concatenate([data, np.expand_dims(rawdata, axis=0)])
#    print(data)
#########################################

md = pd.DataFrame.from_records(recs, columns=recs[0].keys())

# make sure there is one metadata row for each image frame
assert(len(md) == data.shape[0])

# compare checksums
assert(md.loc[0, 'sha1'] == sha1(data[0].ravel()).hexdigest())
assert(md.loc[len(md)-1,'sha1'] == sha1(data[-1].ravel()).hexdigest())

np.save(frames_out, data)
md.to_pickle(metadata_out)

#subprocess.call(['speech-dispatcher'])        #start speech dispatcher
#subprocess.call(['spd-say', '"this process has finished"'])

