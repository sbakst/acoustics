setwd('~/Desktop/vmshare/QP/raw_qp_data/sib_bmps/')


allalphas <- read.delim('../../101-128data_2.csv', header = TRUE, sep = ',') # changed Subject to subject
# allalphas<-na.omit(allalphas[,c(1,2)])
# allalphas <- na.omit(allalphas)

#s_diffs101 <- read.delim('101/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs103 <- read.delim('103/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs104 <- read.delim('104/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs105 <- read.delim('105/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
#s_diffs106 <- read.delim('106/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs108 <- read.delim('108/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs109 <- read.delim('109/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs110 <- read.delim('110/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs112 <- read.delim('112/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs113 <- read.delim('113/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs114 <- read.delim('114/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs115 <- read.delim('115/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs116 <- read.delim('116/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs117 <- read.delim('117/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs118 <- read.delim('118/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs119 <- read.delim('119/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs120 <- read.delim('120/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs121 <- read.delim('../121/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs122 <- read.delim('../122/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs123 <- read.delim('../123/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs124 <- read.delim('../124/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs125 <- read.delim('../125/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs126 <- read.delim('../126/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs127 <- read.delim('../127/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
s_diffs128 <- read.delim('../128/S_mask_diffs_frombase.txt', header = TRUE, sep = '\t')

allsdiffs <- rbind(#s_diffs101, 
                  s_diffs103,
                  s_diffs104, 
                  s_diffs105, 
                 # s_diffs106,
                  s_diffs108,
                  s_diffs109,
                  s_diffs110, 
                  s_diffs112, 
                  s_diffs113, 
                  s_diffs114, 
                  s_diffs115, 
                  s_diffs116, 
                  s_diffs117, 
                  s_diffs118, 
                  s_diffs119, 
                  s_diffs120, 
                  s_diffs121,
                  s_diffs122,
                  s_diffs123, 
                  s_diffs124, 
                  s_diffs125, 
                  s_diffs126, 
                  s_diffs127, 
                  s_diffs128)
allsdata <- merge(allsdiffs, allalphas, by = 'subject')
allsdata$maxfreqsd<-as.numeric(levels(allsdata$maxfreqsd))[allsdata$maxfreqsd]
plot(allsdata$alpha,allsdata$normalized_avg)