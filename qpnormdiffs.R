setwd('~/Desktop/vmshare/QP/raw_qp_data/sub_bmps/')


allalphas <- read.delim('../../101-128data_2.csv', header = TRUE, sep = ',') # changed Subject to subject
# allalphas<-na.omit(allalphas[,c(1,2)])
# allalphas <- na.omit(allalphas)

#r_diffs101 <- read.delim('101/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs103 <- read.delim('103/bar_diffs.txt', header = TRUE, sep = '\t')
r_diffs104 <- read.delim('104/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs105 <- read.delim('105/r_diffs.txt', header = TRUE, sep = '\t')
#r_diffs106 <- read.delim('106/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs108 <- read.delim('108/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs109 <- read.delim('109/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs110 <- read.delim('110/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs112 <- read.delim('112/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs113 <- read.delim('113/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs114 <- read.delim('114/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs115 <- read.delim('115/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs116 <- read.delim('116/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs117 <- read.delim('117/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs118 <- read.delim('118/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs119 <- read.delim('119/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs120 <- read.delim('120/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs121 <- read.delim('../121/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs122 <- read.delim('../122/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs123 <- read.delim('../123/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs124 <- read.delim('../124/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs125 <- read.delim('../125/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs126 <- read.delim('../126/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs127 <- read.delim('../127/r_diffs.txt', header = TRUE, sep = '\t')
r_diffs128 <- read.delim('../128/r_diffs.txt', header = TRUE, sep = '\t')

allr_diffs <- rbind(#r_diffs101, 
                  r_diffs103,
                  r_diffs104, 
                  r_diffs105, 
                  #r_diffs106, 
                  r_diffs108, 
                  r_diffs109, 
                  r_diffs110, 
                  r_diffs112, 
                  r_diffs113, 
                  r_diffs114, 
                  r_diffs115, 
                  r_diffs116, 
                  r_diffs117, 
                  r_diffs118, 
                  r_diffs119, 
                  r_diffs120, 
                  r_diffs121,
                  r_diffs122,
                  r_diffs123, 
                  r_diffs124, 
                  r_diffs125, 
                  r_diffs126, 
                  r_diffs127, 
                  r_diffs128)
alldata <- merge(allr_diffs, allalphas, by = 'subject')
alldata$F3sd<-as.numeric(levels(alldata$F3sd))[alldata$F3sd]

# test for heteroscedasticity
m <- lm(normalized_avg~alpha, data = alldata)
mresid <- m$residuals^2
alphasq <- alldata$alpha^2
alpha <- alldata$alpha
bp <- lm(mresid~alpha+alphasq)
summary(bp)