import os, sys, re
import numpy as np
import pandas as pd
import argparse
import subprocess
from scipy import ndimage
from hashlib import sha1
from collections import OrderedDict
from ultratils.pysonix.scanconvert import Converter
import matplotlib.pyplot as plt

# for PCA business
from sklearn import decomposition
from sklearn.decomposition import PCA

parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Experiment directory containing all subjects")
parser.add_argument("n_components", help="Number of principal components to output")
parser.add_argument("-v", "--visualize", help="Produce plots of PC loadings on fan",action="store_true")
args = parser.parse_args()

# check directory; import frame data and metadata
try:
    expdir = args.directory
except IndexError:
    print("\tDirectory provided doesn't exist")
    ArgumentParser.print_usage
    ArgumentParser.print_help
    sys.exit(2)

data_in = os.path.join(expdir, 'frames.npy')
data = np.load(data_in)
metadata_in = os.path.join(expdir, 'frames_metadata.pickle')
md = pd.read_pickle(metadata_in)
print(data.shape)



# sanity checks
assert(len(md) == data.shape[0]) # make sure one md row for each frame
assert(md.loc[0, 'sha1'] == sha1(data[0].ravel()).hexdigest()) # checksums
assert(md.loc[len(md)-1,'sha1'] == sha1(data[-1].ravel()).hexdigest())

# define mask that subsets frames and frame metadata
# analysis_list = ["IY1", "AE1", "UW1", "OW1", "SH", "IZ1"]
analysis_list = ["EY1", "IY1", "IH1", "OW1", "UW1","UH1", "R"]
mask = md['phone'].isin(analysis_list)
mask = mask.as_matrix()
sub_data = data[mask]
sub_md = md[mask]

# check for appropriate number of components
n_components = int(args.n_components)
if n_components > (sub_data.shape[0] - 1):
    print("EXITING: Number of components requested definitely exceeds number to be produced")
    sys.exit(2)

# preprocess image data (apply median filter and thresholding)
def denoise_and_threshold(rawdata,in_value=90):
    value = in_value
    rawdata_new = np.where((0 - rawdata) > value, 0, rawdata - value)
    rawdata_new2 = np.where((255 - rawdata_new) < value, 255, rawdata_new + value)
    rawdata_new3 = np.where(rawdata_new2 > value, 255, 0)
    return rawdata_new3

# uncommented 64 (just below) and commented lin 66)
sub_data = ndimage.median_filter(sub_data, 1)# we had tried 3
thresh_value = 160 
#sub_data = denoise_and_threshold(sub_data,thresh_value)

# TODO do PCA
pca = PCA(n_components=n_components) 
frames_reshaped = sub_data.reshape([
                    sub_data.shape[0],
                    sub_data.shape[1] * sub_data.shape[2]
                    ])
pca.fit(frames_reshaped)
analysis = pca.transform(frames_reshaped)

# PCA text outputs
subject = re.sub("[^0-9]", "", expdir)
pc_filestring = "tamil{:}_pc_out.txt".format(subject)
pc_out = os.path.join(expdir,pc_filestring) 
pc_headers = ["pc"+str(i+1) for i in range(0,n_components)]
meta_headers = sub_md.columns[0:11].values
headers = list(meta_headers) + pc_headers # TODO this isn't working
metadata = sub_md.as_matrix(columns = sub_md.columns[0:11]) 
print("Metadata",type(metadata),metadata.shape)
print("Analysis",type(analysis),analysis.shape)
d = np.row_stack((headers,
                    np.column_stack((metadata, analysis))
                    ))
np.savetxt(pc_out, d, fmt="%s", delimiter = ',')

# visualize components on (approximately) converted fan, if desired
if args.visualize:
    # initialize a converter object
    class Header(object):
        def __init__(self):
            pass

    class Probe(object):
        def __init__(self):
            pass
    
    header = Header()
    print(header)
    header.w = 127          # input image width
    header.h = 255#1020 #255          # input image height
    header.sf = 4000000     # for ultrasonix this is the vec-freq value
    probe = Probe()
    probe.radius = 10000    # based on '10' in transducer model number
    probe.numElements = 128 # based on '128' in transducer model number
    probe.pitch = 185 #205  # guess based on Ultrasonix C9-5/10 transducer
    c = Converter(header, probe)

    image_shape = (255,127)#(1020,127)#(255,127)

    for n in range(0,n_components):
        d = pca.components_[n].reshape(image_shape)
        mag = np.max(d) - np.min(d)
        d = (d-np.min(d))/mag*255
        pcn = np.flipud(c.as_bmp(np.flipud(d)))
        #pcn = c.as_bmp(d)
        plt.title("PC{:} min/max loadings, Subj. {:}".format((n+1),subject))
        plt.imshow(pcn, cmap="Greys_r")
        file_ending = "subj{:}-pc{:}.pdf".format(subject, (n+1))
        savepath = os.path.join(expdir,file_ending)
        plt.savefig(savepath)

print("Data saved. Explained variance ratio of PCs: {:}".format(str(pca.explained_variance_ratio_)))

#subprocess.call(['speech-dispatcher'])        #start speech dispatcher
#subprocess.call(['spd-say', '"this process has finished"'])
#subprocess.call(['killall', 'speech-dispatcher'])
