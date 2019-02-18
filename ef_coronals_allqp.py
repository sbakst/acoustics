#!/usr/bin/env python

## USAGE: ef_coronals.py SOURCE DESTINATION TARGETS (TARGETWORDS) 

import sys
import os
import re
import subprocess
import audiolabel
import numpy as np
import shutil
import datetime

# DIR is the directory containing folders of images, audio, etc

PARENTDIR = sys.argv[1]
DESTDIR = sys.argv[2]
TAR = sys.argv[3]
WORDFILE = sys.argv[4]

w = open(WORDFILE)
WORDS = str.split(w.read(),"\n")



subjects = [d for d in os.listdir(PARENTDIR) if os.path.isdir(os.path.join(PARENTDIR,d))]

# create a logfile
# time = datetime.datetime.now()
# ts = str(time.year*100000000+time.month*1000000+time.day*10000+time.hour*100+time.minute)
# logfile = "./log_"+ts+".txt"
# log = open(logfile,'wb')

for s in [102] : #np.arange(120,121) : # loop over all subjects; used to say for s in range(0,len(subjects))
#for s in [104, 105, 108, 110, 112, 113, 114, 115, 116, 117, 118, 119, 120]:
    subject = str(s) # str(subjects[s]) 
    DIR = os.path.join(PARENTDIR,subject)
    subs = os.listdir(DIR)
    print (DIR)
    for i in range (len(subs)) : # loop over all acquisitions from SUBJECT
#   for i in range (10) :
        thepath=subs[i]
        print(thepath)
        if thepath.startswith('.'):
            continue
        if thepath == '2015-04-27T183707-0700':# this was done by hand, 116-130.
            continue 
      # generate names of all necessary files
        stimfile = DIR+"/"+subs[i]+"/stim.txt"
        wavfile = DIR+"/"+subs[i]+"/"+subs[i]+".bpr.ch1.wav"
        outfile = DIR+"/"+subs[i]+"/"+subs[i]+".TextGrid"
        syncfile = DIR+"/"+subs[i]+"/"+subs[i]+".bpr.sync.txt"
        
        if re.match('\d{4}-\d{2}-\d{2}T\d{6}-\d{4}',subs[i]) and (os.path.isfile(syncfile)) and (os.path.isfile(stimfile)) and (os.path.isfile(wavfile)) :
        
            stim = open(stimfile)
            stimtext = stim.read()
        
            for w in WORDS :
                w_reg = re.escape(w)+r"$"
                if w and re.findall(w_reg,stimtext) :

                    if not (os.path.isfile(outfile)) :
                        proc = subprocess.check_call(['pyalign',wavfile,stimfile,outfile])

                    fname = os.path.splitext(outfile)[0]
                    pm = audiolabel.LabelManager(from_file = outfile, from_type='praat')

                    for lab in pm.tier('phone') :
                        if (re.match(TAR,lab.text)) :
                            label = lab.text
                            t1 = lab.t1
                            t2 = lab.t2
                            #mid = (t2+t1)/2
                            if s < 128:
                                sm = audiolabel.LabelManager(from_file=syncfile, from_type='table',fields_in_head=False, fields='t1,frameidx')

    # this assumes all the bprs have already been converted to bmps.
                                fr1 = sm.tier('frameidx').label_at(t1).text # assuming that tgs have been corrected, this is where we want to start
                                fr2 = sm.tier('frameidx').label_at(t2).text # assuming that tgs have been corrected, this is where we want to stop
                            else:
                                sm = audiolabel.LabelManager(from_file=syncfile, from_type='table',fields_in_head=True)
    # this assumes all the bprs have already been converted to bmps.
                                
                                fr1 = sm.tier('pulse_idx').label_at(t1).text # assuming that tgs have been corrected, this is where we want to start
                                fr2 = sm.tier('pulse_idx').label_at(t2).text # assuming that tgs have been corrected, this is where we want to stop
                                
                                if fr1 == 'NA':
                                    fr1 = sm.tier('pulse_idx').label_at(t1 + 0.01).text # assuming that tgs have been corrected, this is where we want to start
                                if fr2 == 'NA':
                                    fr2 = sm.tier('pulse_idx').label_at(t2 + 0.01).text # assuming that tgs have been corrected, this is where we want to start
                            if not (os.path.isfile(os.path.splitext(outfile)[0] + '.' + str(fr1)+'.bmp')):
                                print(os.path.splitext(outfile)[0] + '.' + str(fr1)+'.bmp')
                                probe = '19'
                                dir1 = os.path.join(DIR,subs[i])
                             # call bpr2bmp
                            
                                bmp_proc = subprocess.check_call(['bpr2bmp','--probe',probe,'--seek',dir1])

                            newdir=DESTDIR +"/" +subject+"/"+subs[i]
                            if not os.path.exists(newdir) :
                                os.makedirs(newdir)

                            missing = ""
                            for j in range(int(fr1),int(fr2)+1) :
                                tocopy = DIR+"/"+subs[i]+"/"+subs[i]+"."+str(j)+".bmp"
                                dest = newdir+"/"+subs[i]+"."+str(j)+".bmp"
                                if (os.path.isfile(tocopy)) :
                                    shutil.copyfile(tocopy,dest)
                                else :
                                    missing = missing + tocopy + ";"

                            dest = newdir+"/stim.txt"
                            shutil.copyfile(stimfile,dest)
                            dest = newdir+"/"+subs[i]+".wav"
                            shutil.copyfile(wavfile,dest)
                            dest = newdir+"/"+subs[i]+".bpr.sync.txt"
                            shutil.copyfile(syncfile,dest)
                            dest = newdir+"/"+subs[i]+".TextGrid"
                            shutil.copyfile(outfile,dest)
    #                     log.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:0.3f}\t{}\n".format(subject, subs[i], stimtext, label, fr1, fr2, frm, mid,missing))
    #                     print w + "\t" + stimtext
