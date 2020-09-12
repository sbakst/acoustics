#!/usr/bin/env python
'''\nfricative_analysis.py
fricative_analysis.py: Do spectral analysis on fricative portions of soundfiles in specified directory tree.

Usage: fricative_analysis.py dirname

Arguments:
  dirname   Root of directory tree containing soundfiles to be analyzed.

fricative_analysis.py loops over all of the .wav files in a specified directory and all of
its subdirectories. It looks in the .wav files' associated textgrids (any textgrids with the
same name as the .wav file but with a .TextGrid extension) for a tier named 'phone', performs
spectral analysis during portions with a fricative label, and prints the analysis.
'''

# Authors: Keith Johnson (keithjohnson@berkeley.edu)
# 
# Copyright (c) 2014, The Regents of the University of California
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# * Neither the name of the University of California nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os, sys
import subprocess
import audiolabel
import re
import matplotlib.pyplot as plt 

#-----------------------------
# This script uses the following functions from the ESPS library for speech processing
#   - hditem: get information from a wav or fea file header
#   - fft: spectral analysis, producing a single spectrum or a spectrogram
#   - melspec: convert an fft spectrum into a "Mel transformed" auditory spectrum
#   - pplain: print values from fea files (spectra) to plain text

# Read about these and other ESPS routines in the Berkeley Phonetics Machine 
# using the 'man' command.   For example:
#       > man fft 
# will show the manual page for the fft routine
#-------------------------


def usage():
    print sys.exit(__doc__)

### Additions for my qp

onset = ["rah", "Rome", "ream", "SEW", "SOB", "SEA"]
coda = ["bar", "bore", "beer", "BOSS", "PIECE", "DOSE"]
aset = ["rah", "bar", "BOSS", "SOB"]
oset = ["bore", "DOSE", "Rome", "SEW"]
###



head = '\t'.join(('frame phone word vowel position max maxfreq').split())
print head

# spectral computation parameter
window = 0.01 # window size in seconds

# Here is a definition of the phonetic symbols we will analyze in this script
fricatives = re.compile("S|SH|F|V|TH|DH")

# I use a "shared" directory in the Phonetics machine that I call BPM
# you can set up a shared directory in the Virtual Box settings for your machine
try:
    directory = sys.argv[1]
except IndexError:
    usage()
    sys.exit(2)

# the following loop looks up all of the sound files in the target directory, 
#   looks for a TextGrid counter-part to the sound file
#   looks for fricatives in a tier called "phone" in the TextGrid
#   computes a spectrum from the midpoint of the fricative 
#   displays the spectrum and computes something about the spectrum
#   -- comment out the spectrum plotting if you have lots of files
#   -- save the results print out by redirecting the script output to a file:
#          $ python fricative_analysis.py > myresults.txt

for root,dirs,files in os.walk(directory):  # walk the directory
    for soundfile in files:     # check each sound file
        if not '.wav' in soundfile:  # if not a .wav, go on to the next file
            continue   
            
                 
        frame = os.path.splitext(soundfile)[0]  # get the frame name
        tg = os.path.join(root,os.path.splitext(soundfile)[0]+'.TextGrid')  # expect a TextGrid file
	       
        sf=float(subprocess.check_output(["hditem", "-i","record_freq", os.path.join(root,soundfile)]))
        
        wsamps = int(round(window*sf))  # window size in samples
        
        # open the praat text grid        
        pm = audiolabel.LabelManager(from_file=(tg),from_type="praat")
         
        # loop through all of the labels on the "phone" tier that match the set of fricatives
        for f in pm.tier('phone').search(fricatives):
            # from the label get the location of the midpoint, and the text of the label
            
            midpoint = f.center
            
            phone = f.text
            
            #calculate a spectrum - store in temp.spec
            start = int(round((midpoint*sf) - (wsamps/2)))  # scoot over 1/2 a window from the center
            ret=subprocess.call(['fft','-z','-wHamming','-l{:d}'.format(wsamps),'-S{:d}'.format(wsamps),
                                    '-r{:d}:+1'.format(start),(os.path.join(root,soundfile)),'temp1.fea'])
            
            # calculate the mel frequency spectrum from the raw FFT - from 300 Hz to 1/2 the sampling freq
            nyquist = sf/2
            ret=subprocess.call(["melspec","-H300:{}".format(nyquist),"-n60","temp1.fea","temp2.fea"])
                
            # read the spectrum into this script    
            spectrum_string = subprocess.check_output(["pplain","-fre_spec_val","temp2.fea"])
            spec = spectrum_string.rstrip().split(' ')    # break the string into separate values
            freq_string = subprocess.check_output(["hditem","-ifreqs","temp2.fea"])
            freq = freq_string.rstrip().split(' ')
            spectrum = map(float,spec)           # convert array from string to floating point number
            maxx=0
            maxfreq = 0
            for i in range(10, 59):
                if spectrum[i] > maxx:
                    maxfreq = freq[i]
                    maxx = spectrum[i]
            
            
            #low=sum(spectrum[0:19])/20              # sum the amplitudes in the bottom half
            #high = sum(spectrum[40:59])/20          # and those in the top half
            #hl_diff = high-low                  # take the ratio of the amplitudes

            word_label = pm.tier('word').label_at(midpoint)
            word = word_label.text

            if word in onset:
                position = "onset"
            else:
                position = "coda"
            if word in aset:
                vowel = "a"
            else:
                if word in oset:
                    vowel = "o"
                else:
                    vowel = "i"


            # this line prints results
            print frame, phone, str(word), str (vowel), position, maxx, maxfreq

            
            # -------------------- code to show a spectrum plot ----------------
            #freq_string = subprocess.check_output(["hditem","-ifreqs","temp2.fea"])
            #freq = freq_string.rstrip().split(' ')
            #frequency = map(float,freq)
            #smax = max(spectrum)               # a useful number to have for plotting the text label

            #fig = plt.figure(1)
            #plt.plot(freq,spectrum,color="black")
            #plt.plot(freq[0:19],spectrum[0:19],color="blue")
            #plt.plot(freq[40:60],spectrum[40:60],color="red")
            #plt.xlabel('Frequency')
            #plt.ylabel('Amplitude')
            #plt.grid(True)
            #plt.text(100,smax-2,"H-L = " + str(hl_diff))
            #plt.show()         
            # ------------------- end of plotting code ----------------------
            

