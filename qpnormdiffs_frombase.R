setwd('~/Desktop/vmshare/QP/raw_qp_data/sub_bmps/')


allalphas <- read.delim('../../101-128data_2.csv', header = TRUE, sep = ',') # changed Subject to subject
# allalphas<-na.omit(allalphas[,c(1,2)])
# allalphas <- na.omit(allalphas)

#r_diffs_frombase101 <- read.delim('101/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase102 <- read.delim('102/R_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase103 <- read.delim('103/r_mask_diffs_frombase_retroflex.txt', header = TRUE, sep = '\t')
r_diffs_frombase104 <- read.delim('104/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase105 <- read.delim('105/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase106 <- read.delim('106/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase108 <- read.delim('108/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase109 <- read.delim('109/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase110 <- read.delim('110/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase112 <- read.delim('112/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase113 <- read.delim('113/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase114 <- read.delim('114/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase115 <- read.delim('115/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase116 <- read.delim('116/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase117 <- read.delim('117/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase118 <- read.delim('118/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase119 <- read.delim('119/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase120 <- read.delim('120/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase121 <- read.delim('../121/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase122 <- read.delim('../122/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase123 <- read.delim('../123/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase124 <- read.delim('../124/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase125 <- read.delim('../125/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase126 <- read.delim('../126/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase127 <- read.delim('../127/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')
r_diffs_frombase128 <- read.delim('../128/r_mask_diffs_frombase.txt', header = TRUE, sep = '\t')

allr_diffs_frombase <- rbind(r_diffs_frombase102, 
                  r_diffs_frombase103,
                  r_diffs_frombase104, 
                  r_diffs_frombase105, 
                  r_diffs_frombase106, 
                  r_diffs_frombase108, 
                  r_diffs_frombase109, 
                  r_diffs_frombase110, 
                  r_diffs_frombase112, 
                  r_diffs_frombase113, 
                  r_diffs_frombase114, 
                  r_diffs_frombase115, 
                  r_diffs_frombase116, 
                  r_diffs_frombase117, 
                  r_diffs_frombase118, 
                  r_diffs_frombase119, 
                  r_diffs_frombase120, 
                  r_diffs_frombase121,
                  r_diffs_frombase122,
                  r_diffs_frombase123, 
                  r_diffs_frombase124, 
                  r_diffs_frombase125, 
                  r_diffs_frombase126, 
                  r_diffs_frombase127, 
                  r_diffs_frombase128)
alldata <- merge(allr_diffs_frombase, allalphas, by = 'subject')
alldata$F3sd<-as.numeric(levels(alldata$F3sd))[alldata$F3sd]

# test for heteroscedasticity
m <- lm(normalized_avg~alpha, data = alldata)
mresid <- m$residuals^2
alphasq <- alldata$alpha^2
alpha <- alldata$alpha
bp <- lm(mresid~alpha+alphasq)
summary(bp)

plot(alldata$alpha, alldata$normalized_avg)

cor.test(alldata$alpha, alldata$normalized_avg)