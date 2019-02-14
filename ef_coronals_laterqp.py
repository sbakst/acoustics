#!/usr/bin/env python

## USAGE: ef_coronals.py SOURCE DESTINATION TARGETS (TARGETWORDS) 
## use with a bce that has NOT BEEN UPDATED from summer 2016

# for subjects above 120

import sys
import os
import re
import subprocess
import audiolabel
import numpy as np
import shutil
import datetime
#import Image
from PIL import Image
import matplotlib.pyplot as plt
from ultratils.exp import Exp


# DIR is the directory containing folders of images, audio, etc

PARENTDIR = sys.argv[1]
#DESTDIR = sys.argv[2]
#TARGETS = sys.argv[3] # In this script, these refer to vowels--a vowel file I suppose?
WORDFILE = sys.argv[2]
SUB = sys.argv[3]

w = open(WORDFILE)
WORDS = str.split(w.read(),"\n")
tarlist = []
tarlist =['S']

#subjects = [d for d in os.listdir(PARENTDIR) if os.path.isdir(os.path.join(PARENTDIR,d))]

# create a logfile
time = datetime.datetime.now()
ts = str(time.year*100000000+time.month*1000000+time.day*10000+time.hour*100+time.minute)
logfile = "./log_"+ts+".txt"
log = open(logfile,'wb')

expdir = os.path.join(PARENTDIR,SUB)

e = Exp(expdir=expdir)
e.gather()
for a in e.acquisitions:
   print(a)
   tp = a.timestamp
   if tp == '2015-11-17T103748-0800':
      continue
#   print(tp)
   stimfile = a.abs_stim_file 
   stim = open(stimfile)
   stimtext = stim.read()
#   outfile = a.abspath+'/'+tp+'.TextGrid'
   outfile = expdir+'/'+tp+'/'+tp+'.TextGrid'
   print(outfile)
   bprtg = expdir+'/'+tp+'/'+tp+'.bpr.sync.TextGrid'
   if not os.path.isfile(bprtg):
      continue
      print('no bpr sync')
   print (outfile)
   print ('gonna give you')
   fname = os.path.splitext(outfile)[0]
   print (fname)
   pm = audiolabel.LabelManager(from_file = outfile, from_type='praat')
   for plab in pm.tier('phone') :
      if plab.text == 'S' : 
         print (plab.text)
         frt1 = plab.t1
#         print frt1
#         print 'that was frt1'
         frt2 = plab.t2

         for l in a.pulse_idx.tslice(t1=frt1, t2 = frt2):
            # l should look something like Label( t1=0.3947, t2=0.4035, text='b'21'' )
            print(l)
            try:
               d = a.frame_at(l.t1, convert = True)
               if d is None:
                  continue
                  print('skipping' + str(l.text))
            except ValueError:
               print('what is a struct.error')
               continue 
            print(d)
            print('this is l.text' + str(l.text))
#            someidx = int(a.sync_lm.tier('raw_data_idx').label_at(l.t1).text)
#            d=a.image_converter.as_bmp(d)
#         pfrm =  np.flipud(e.acquisitions[int(frmid)].image_converter.as_bmp(d)) 
#         print(pcn)


            d = np.flipud(d).astype(np.uint8)
           
#            print (a)
#            printme = a.image_converter.as_bmp(d).astype(np.uint8)
#            print(printme)

            frame = Image.fromarray(d)
            imgname = '{:}.{:}.bmp'.format(a.timestamp, l.text)

            print(frame)
            print(imgname)
            print(a.abspath)
            frame.save(os.path.join(a.abspath, imgname))
#               a.pulse_idx.tslice(t1=frt1, t2=frt2)
#               for l in a.pulse_idx.tslice(frmid): 
#               d = a.frame_at(frmid, convert = True)
#               if d is None:
#                  prevl = a.pulse_idx.prev(l)
#                  d = a.frame_at(prevl.center, convert = True)
#                  print "WARNING: duplicated frame {:} in {:}.".format(
#                     prevl.text, slice[0]
#                  )
#               d = np.flipud(d).astype(np.uint8)
#               print d
#               print a
#               frame = Image.fromarray(d)
#               imgname = '{:}.{:}.{:}.bmp'.format(a.timestamp, a.text, stimtext)
#               frame.save(os.path.join(a.abspath, imgname))

#                     mid = (t2+t1)/2
#                     sm = audiolabel.LabelManager(fromFile=syncfile, fromType='table',fieldsInHead=False, fields='t1,frameidx')
#                     frm = sm.tier('frameidx').label_at(mid).text
#                     fr1 = sm.tier('frameidx').label_at(t1).text
#                     fr2 = sm.tier('frameidx').label_at(t2).text
#                  if lab.text == 'FOR':
#                     fort = lab.t1()


#                     newdir=DESTDIR +"/" +subject+"/"+subs[i]
#                     if not os.path.exists(newdir) :
#                        os.makedirs(newdir)

#                     missing = ""
#                     for j in range(int(fr1),int(fr2)+1) :
#                        tocopy = DIR+"/"+subs[i]+"/"+subs[i]+"."+str(j)+".bmp"
#                        dest = newdir+"/"+subs[i]+"."+str(j)+".bmp"
#                        if (os.path.isfile(tocopy)) :
#                           shutil.copyfile(tocopy,dest)
#                        else :
#                           missing = missing + tocopy + ";"
#
#                     dest = newdir+"/stim.txt"
#                     shutil.copyfile(stimfile,dest)
#                     dest = newdir+"/"+subs[i]+".wav"
#                     shutil.copyfile(wavfile,dest)

#                     log.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:0.3f}\t{}\n".format(subject, subs[i], stimtext, label, fr1, fr2, frm, mid,missing))
#                     print w + "\t" + stimtext
