setwd('~/Desktop/vmshare/QP/raw_qp_data/sib_bmps/')


allalphas <- read.delim('../../101-128data_2.csv', header = TRUE, sep = ',') # changed Subject to subject
# allalphas<-na.omit(allalphas[,c(1,2)])
allalphas <- na.omit(allalphas)

s_diffs101 <- read.delim('101/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs104 <- read.delim('104/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs105 <- read.delim('105/midframe_s_diffs.txt', header = TRUE, sep = '\t')
#s_diffs106 <- read.delim('106/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs108 <- read.delim('108/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs110 <- read.delim('110/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs112 <- read.delim('112/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs113 <- read.delim('113/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs114 <- read.delim('114/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs115 <- read.delim('115/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs116 <- read.delim('116/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs117 <- read.delim('117/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs118 <- read.delim('118/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs119 <- read.delim('119/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs120 <- read.delim('120/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs121 <- read.delim('../../../../121/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs122 <- read.delim('../../../../122/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs123 <- read.delim('123/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs124 <- read.delim('124/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs125 <- read.delim('125/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs126 <- read.delim('126/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs127 <- read.delim('127/midframe_s_diffs.txt', header = TRUE, sep = '\t')
s_diffs128 <- read.delim('128/midframe_s_diffs.txt', header = TRUE, sep = '\t')

middiffs <- rbind(#s_diffs101, 
  s_diffs104, 
  s_diffs105, 
  #s_diffs106, 
  s_diffs108, 
  s_diffs110, 
  s_diffs112, 
  #s_diffs113, 
  s_diffs114, 
  s_diffs115, 
  s_diffs116, 
  #s_diffs117, 
  #s_diffs118, 
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
midata <- merge(middiffs, allalphas, by = 'subject')
midata$maxfreqsd<-as.numeric(levels(midata$maxfreqsd))[midata$maxfreqsd]
