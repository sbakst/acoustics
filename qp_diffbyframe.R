# framenum = framenum-1 IF subj < 121
# run acoustics first
setwd('~/Desktop/vmshare/QP/raw_qp_data/sub_bmps/')


allalphas <- read.delim('../../101-128data_2.csv', header = TRUE, sep = ',') # changed Subject to subject
# allalphas<-na.omit(allalphas[,c(1,2)])
# allalphas <- na.omit(allalphas)

#r_diffs_frombase101 <- read.delim('101/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase102 <- read.delim('102/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase103 <- read.delim('103/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase104 <- read.delim('104/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase105 <- read.delim('105/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
# r_diffs_frombase106 <- read.delim('106/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase108 <- read.delim('108/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase109 <- read.delim('109/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase110 <- read.delim('110/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase112 <- read.delim('112/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase113 <- read.delim('113/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase114 <- read.delim('114/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase115 <- read.delim('115/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase116 <- read.delim('116/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase117 <- read.delim('117/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase118 <- read.delim('118/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase119 <- read.delim('119/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase120 <- read.delim('120/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase121 <- read.delim('../121/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase122 <- read.delim('../122/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase123 <- read.delim('../123/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase124 <- read.delim('../124/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase125 <- read.delim('../125/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase126 <- read.delim('../126/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase127 <- read.delim('../127/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')
r_diffs_frombase128 <- read.delim('../128/R_mask_diffs_frmtime.txt', header = TRUE, sep = '\t')


allr_diffs_frombase <- rbind(r_diffs_frombase102, 
                             #r_diffs_frombase103,
                             r_diffs_frombase104, 
                             r_diffs_frombase105, 
                             # r_diffs_frombase106, 
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



#r_diffs_frombase101 <- read.delim('101/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all102 <- read.delim('102/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all103 <- read.delim('103/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all104 <- read.delim('104/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all105 <- read.delim('105/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
# r_diffs_all106 <- read.delim('106/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all108 <- read.delim('108/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all109 <- read.delim('109/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all110 <- read.delim('110/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all112 <- read.delim('112/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all113 <- read.delim('113/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all114 <- read.delim('114/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all115 <- read.delim('115/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all116 <- read.delim('116/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all117 <- read.delim('117/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all118 <- read.delim('118/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all119 <- read.delim('119/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all120 <- read.delim('120/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')



r_diffs_all102$alpha<-1.900000
#r_diffs_all103$alpha<-
r_diffs_all104$alpha<-2.300000
r_diffs_all105$alpha<-1.800000
r_diffs_all108$alpha<-2.000000
r_diffs_all109$alpha<-2.700000
r_diffs_all110$alpha<-2.3
r_diffs_all112$alpha<-2.2
r_diffs_all113$alpha<-2.2
r_diffs_all114$alpha<-2.0
r_diffs_all115$alpha<-1.7
r_diffs_all116$alpha<-2.0
r_diffs_all117$alpha<-1.9
r_diffs_all118$alpha<-1.8
r_diffs_all119$alpha<-2.0
r_diffs_all120$alpha<-1.9









rdiffs120<-rbind(#r_diffs_all101,
  r_diffs_all102,
#  r_diffs_all103,
  r_diffs_all104,
  r_diffs_all105,
  # r_diffs_all1
  r_diffs_all108,
  r_diffs_all109,
  r_diffs_all110,
  r_diffs_all112,
  r_diffs_all113,
  r_diffs_all114,
  r_diffs_all115,
  r_diffs_all116,
  r_diffs_all117,
  r_diffs_all118,
  r_diffs_all119,
  r_diffs_all120)

rdiffs120$framenum <-rdiffs120$framenum-1 # i think

r_diffs_all121 <- read.delim('../121/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all122 <- read.delim('../122/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all123 <- read.delim('../123/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all124 <- read.delim('../124/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all125 <- read.delim('../125/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all126 <- read.delim('../126/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all127 <- read.delim('../127/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')
r_diffs_all128 <- read.delim('../128/R_mask_diffs_byts_frmtime.txt', header = TRUE, sep = ',')

r_diffs_all121$alpha<-1.789223
r_diffs_all122$alpha<-2.063465 
r_diffs_all123$alpha<-1.951767
r_diffs_all124$alpha<-1.821079
r_diffs_all125$alpha<-1.982110
r_diffs_all126$alpha<-2.457076
r_diffs_all127$alpha<-2.280988
r_diffs_all128$alpha<-1.910943





rdiffs128<-rbind(r_diffs_all121, 
r_diffs_all122,
r_diffs_all123,
r_diffs_all124,
r_diffs_all125,
r_diffs_all126,
r_diffs_all127,
r_diffs_all128)

rdiffs <-rbind(rdiffs120,rdiffs128)








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