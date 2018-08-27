


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


#rawfile_glob_exp = ('218/block1/')

globlist=(glob.glob('./218/block2/*/*.raw'))

idxlist = []

for rf in globlist:
    print(rf)
    barename = os.path.splitext(rf)[0]
   # TODO proof this against breaking if first acq doesn't have img.txt
    try:
        myframesize = get_frame_size(barename + '.img.txt')
#        print(myframesize)
    except ValueError:
        print("\tUsed previous frame's size, .img.txt not generated")
        myframesize = 129540
        pass
    idxfile = barename + '.idx.txt' # file should just be a column of numbers
    with open(idxfile) as infile:
        idxreader = csv.reader(infile)
        d = list(idxreader)
        rows = sum(1 for row in d)
        last_idx = int((d[rows-1])[0])
    print(last_idx)
    consec_diffs = []
    for f in np.arange(int(last_idx/5),int(4*last_idx/5)):
        minuend = get_frame(rf,f+1, myframesize, med_filter = True)
        subtrahend = get_frame(rf, f, myframesize, med_filter = True)
        cdiff = minuend-subtrahend
        cdiffval = np.linalg.norm(cdiff)
        consec_diffs.append(cdiffval)
#    print(consec_diffs)
    mindiff,mindiff_idx = min((val,idx) for (idx,val) in enumerate(consec_diffs))
    print(mindiff)
    print(mindiff_idx+int(last_idx/5))
    idxlist.append(last_idx)
print(idxlist)
print(np.median(idxlist))

