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

# define frame-getting functions.
def read_metadata(mfile):
    '''Read the metadata from an .img.txt file.'''
    mdict = {}
    with open(mfile, 'r') as f:
        k = f.readline().strip().split("\t")
        v = f.readline().strip().split("\t")
        for fld,val in zip(k, v):
            mdict[fld] = int(val)
    return mdict

def get_frame_size(mfile):
    mdict = read_metadata(mfile)
    frame_size_in_bytes = mdict['nInBufferLen']
    return frame_size_in_bytes

def get_frame(fname, framenum, frame_size, med_filter=False):
    """Requires raw filename, desired frame number, and size of frame as arguments."""
    framesize = frame_size
    data_fmt = 'I' * np.int(framesize / 4) 
    x_start = framesize * framenum
    frame_dim_1 = 127
    frame_dim_2 = 255
    with open(fname, 'rb') as fh:
        fh.seek(x_start) # define starting point at leftmost X coord.
        packed_data = fh.read(framesize) # take a frame-sized chunk at this coord.
        unpacked_data = struct.unpack(data_fmt, packed_data) 
        data = np.array(unpacked_data)
        rdata = np.flipud(data[np.arange(frame_dim_1*frame_dim_2)].reshape(frame_dim_1,frame_dim_2).T)
        if med_filter == True:
            rdata = ndimage.median_filter(rdata, 10)
        return(rdata)

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
frames_out = os.path.join(expdir,"frames.npy")
metadata_out = os.path.join(expdir,"frames_metadata.pickle")

# create the dictionary

tuplist = []
for dirs, subdir, files in os.walk(args.acousticsdir):
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
for dirs, subdir, files in os.walk(expdir):
    for tp in subdir:
# split up timestamp
        tptup = re.split('-',tp)
        ultras.append(tuple(tptup))

# split order by 3rd element, time

ordultras = []
ordultras = sorted(ultras, key = lambda x: x[2])

# make dictionary for looking up wav given timestamp

tamildict = dict(zip(ordultras,ordtrials))
print(tamildict)
#dictitems = tamildict.items()
#tamtam = (sorted(dictitems, key= lambda x:int(x[0][3])))


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

    # TODO proof this against breaking if first acq doesn't have img.txt
    try:
        myframesize = get_frame_size(barename + '.img.txt')
        print(myframesize)
    except ValueError:
        print("\tUsed previous frame's size, .img.txt not generated")
        myframesize = 129540
        pass

    # use dictionary linking wav files to ultrasound to look up wav file
       
    for key,value in tamildict.items():
        print(key)
        joined_tp = ('-'.join(key))
        print(joined_tp)
        if joined_tp == barename:
            print(barename)
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
    idxfile = barename + 'idx.txt' # file should just be a column of numbers
    with open(idxfile) as infile:
        idxreader = csv.reader(infile)
        d = list(idxreader)
        rows = sum(1 for row in d)
        last_idx = int((d[rows-1])[0])
        print(last_idx)

    consec_diffs = []
    for f in np.arange(0,last_idx):
        minuend = get_frame(rf,f+1, framesize, med_filter = True)
        subtrahend = get_frame(rf, f, framesize, med_filter = True)
        cdiff = minuend-subtrahend
        consec_diffs.append(cdiff)
    print(consec_diffs)
    mindiff,mindiff_idx = min((val,idx) for (idx,val) in enumerate(consec_diffs))
    print(mindiff)
    print(mindiff_idx)





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
        #print("\t","phone",phone,"stim",stim,"targ",targ,"before",before,"after",after)

        # get midpoint time and find closest ultrasound frame in sync TG
        # TODO more efficient to duplicate ultratils frame_at approach
#        mid_timepoint = v.center
#        diff_list = []
#        diff2_list = []
#        for frame in sync_pm.tier('pulse_idx'):
#            diff = abs(frame.t1 - mid_timepoint)
#            diff_list.append(diff)
#        for frame in sync_pm.tier('raw_data_idx'):
#            diff2 = abs(frame.t1 - mid_timepoint)
#            diff2_list.append(diff2)
#        mid_pulse_idx_num = min(enumerate(diff_list), key=itemgetter(1))[0] 
#        mid_raw_data_idx_num = min(enumerate(diff2_list), key=itemgetter(1))[0] 

        # get frame, and check for NaN frames
#        change = 0
#        discard_acq = False
#        while True:
#            pre_rawdata = get_frame(rf,mid_pulse_idx_num,myframesize,med_filter=False)
#            if pre_rawdata is None:
#                mid_frame_num -= 1
#                change += 1
#                if change > threshhold:
#                    with open(logfile, "a") as log:
#                        log.write(acq+"\t"+stim+"\t"+phone+"\t"+"discarded"+"\t"+"passed threshhold")
#                    print("Frame change threshhold passed; acq {} discarded".format(acq))
#                    discard_acq = True
#                    break
#                else:
#                    pass
#            else:
#                if change > 0:
#                    with open(logfile, "a") as log:
#                        log.write(acq+"\t"+stim+"\t"+phone+"\t"+"changed by {:}".format(change)+"\t"+"N/A")
#                    print("Changed target in {:} by".format(acq), change, "frames")
#                break

#        # discard the acquisition if needed
#        if discard_acq:
#            shutil.copytree(parent, os.path.join(discard_folder,acq))
#            shutil.rmtree(parent)
#            continue 

        # preprocessing of images
        rawdata = pre_rawdata.astype(np.uint8)

        # generate metadata object for the current acquisition
        recs.append(
            OrderedDict([
                ('timestamp', acq),
                ('time', v.center),
                ('pulseidx', int(mid_pulse_idx_num)),
                ('rawdataidx', int(mid_raw_data_idx_num)),
                ('width', frame_dim_1),
                ('height', frame_dim_2),
                ('phone', phone),
                ('stim', stim),
                ('targ', targ),
                ('before', before),
                ('after', after),
                ('sha1', sha1(rawdata.ravel()).hexdigest()),
                ('sha1_dtype', rawdata.dtype)
            ])
        )

        # add frame to frames list
        if data is None:
            data = np.expand_dims(rawdata, axis=0)
        else:
            data = np.concatenate([data, np.expand_dims(rawdata, axis=0)])

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

