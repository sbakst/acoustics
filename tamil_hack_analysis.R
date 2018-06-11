# using tamil 2/3 into the file to estimate f3
# I think you used rformant to create this.
# later you will do this by ultrasound frame, hopefully.

# negative values of adaptation are following, but over 1 = probably a tracking error, or compensated more than the perturbation.

tam23rds<-read.delim('/Users/Sarah/Box Sync/Berkeley/Dissertation/India/tamil_corrected_twothirds_buffs.txt', header = TRUE, sep = '\t')

tam23rds$buf_file<-as.character(tam23rds$buf_file)
split1<-strsplit(tam23rds$buf_file, "F3|_|.wav")
tam23rds$buffer<-as.factor(sapply(split1, "[", 4))
tam23rds$sub<-as.factor(sapply(split1, "[", 5))
tam23rds$trial<-as.numeric(sapply(split1, "[", 6))

## add bark

tam23rds$f1b<-(13*atan(0.00076*tam23rds$f1))+(3.5*atan((tam23rds$f1/7500)^2))
tam23rds$f2b<-(13*atan(0.00076*tam23rds$f2))+(3.5*atan((tam23rds$f2/7500)^2))
tam23rds$f3b<-(13*atan(0.00076*tam23rds$f3))+(3.5*atan((tam23rds$f3/7500)^2))

tam23rds_in<-subset(tam23rds,buffer=='inbuf')
tam23rds_out<-subset(tam23rds,buffer=='outbuf')



tam23rds_in$group<-cut(tam23rds_in$trial, 18)
ar18_med<-aggregate(cbind(f1,f1b,f2,f2b,f3,f3b)~group+sub, data = tam23rds_in, median)

mindf<-NULL
for(i in (1:length(levels(tam23rds_in$sub)))){
  subj = levels(tam23rds_in$sub)[i]
  tempar<-NULL
  tempar2<-NULL
  minf3<-NULL
  trif3<-NULL
  mbf3<-NULL
  tempar<-subset(ar18_med, sub == subj)
  tempar2<-subset(tam23rds_in, sub == subj)
  minf3 = min(tempar$f3[4:15])
  mbf3 = median(subset(tempar2, trial<16)$f3)  
  trif3 = (tempar2$trial[which(tempar2$f3==minf3)])
{ if (length(trif3) > 1)
  trialf3 = trif3[-1]
  else
    ( trialf3 = trif3)}    
f2match = (tempar2$f2[which(tempar2$trial==trialf3)])
arf2base = median(subset(tempar2, trial<16)$f2)
#    cat("For subject", subj, "here's the minimum f3:", minf3,"at trial", trial, "where f2 was", f2match,"\n")
mindf$sub[i]<-(subj)
mindf$mintrialf3[i]<-as.numeric(trialf3)
mindf$minf3[i]<-as.numeric(minf3)
mindf$f2match[i]<-as.numeric(f2match)
mindf$arf2base[i]<-as.numeric(arf2base) # what is the basleine f2 in the f3 block
mindf$basef3[i]<-as.numeric(mbf3)
}

mindf<-data.frame(mindf)



mindf$f3max<-mindf$basef3+300

mindf$f3diff<-mindf$basef3-mindf$minf3
mindf$f2diff<-mindf$arf2base-mindf$f2match

mindf$f3perc<-(mindf$f3diff/300)

mindf$oppf2d<-mindf$f2match-mindf$arf2base








f3inbuf<-subset(tam23rds, buffer=='inbuf')
f3outbuf<-subset(tam23rds, buffer=='outbuf')






# remember, no idx because not ifcformant!

f3in_ave <- aggregate(f3~trial, data = f3inbuf, median)
f3out_ave <- aggregate(f3~trial, data = f3outbuf, median)

plot(f3out_ave, ylim = c(2100,2800))
points(f3in_ave, col = 'red')


plot(f3~trial, subset(tam23rds, sub=='207' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='212' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='213' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='215' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='218' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='218' & buffer == 'inbuf'), ylim=c(1800,2300))
plot(f3~trial, subset(tam23rds, sub=='219' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='220' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='221' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='221' & buffer == 'inbuf'), ylim=c(1800,2300))
plot(f3~trial, subset(tam23rds, sub=='221' & buffer == 'inbuf'), ylim=c(1800,2600))
plot(f3~trial, subset(tam23rds, sub=='222' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='222' & buffer == 'inbuf'), ylim=c(1800,2600))
plot(f3~trial, subset(tam23rds, sub=='224' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='224' & buffer == 'inbuf'), ylim=c(1800,2600))
plot(f3~trial, subset(tam23rds, sub=='227' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='230' & buffer == 'inbuf'))
points(f3~trial, subset(tam23rds, sub=='230' & buffer == 'outbuf'), col = 'red')
plot(f3~trial, subset(tam23rds, sub=='231' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='232' & buffer == 'inbuf'))
points(f3~trial, subset(tam23rds, sub=='232' & buffer == 'outbuf'), col = 'red')
plot(f3~trial, subset(tam23rds, sub=='233' & buffer == 'inbuf'))
abline(v=15)
plot(f3~trial, subset(tam23rds, sub=='234' & buffer == 'inbuf'))
points(f3~trial, subset(tam23rds, sub=='234' & buffer == 'outbuf'), col = 'red')
plot(f3~trial, subset(tam23rds, sub=='235' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='236' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='236' & buffer == 'inbuf'), ylim=c(1800,2600))
points(f3~trial, subset(tam23rds, sub=='236' & buffer == 'outbuf'), col = 'red')
plot(f3~trial, subset(tam23rds, sub=='237' & buffer == 'inbuf'))
plot(f3~trial, subset(tam23rds, sub=='237' & buffer == 'inbuf'), ylim=c(2200,2700))



