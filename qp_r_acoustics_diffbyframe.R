# R acoustics

setwd('~/Desktop/vmshare/QP/raw_qp_data/sub_bmps/')

#r_acoustics101 <- read.delim('101/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics102 <- read.delim('102/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics102$subject <- 102
r_acoustics103 <- read.delim('103/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics103$subject <- 103
r_acoustics104 <- read.delim('104/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics104$subject <- 104
r_acoustics105 <- read.delim('105/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics105$subject <- 105
# r_acoustics106 <- read.delim('106/Racoustic_data.txt', header = TRUE, sep = ',')
#r_acoustics102$subject <- 106
r_acoustics108 <- read.delim('108/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics108$subject <- 108
r_acoustics109 <- read.delim('109/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics109$subject <- 109
r_acoustics110 <- read.delim('110/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics110$subject <- 110
r_acoustics112 <- read.delim('112/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics112$subject <- 112
r_acoustics113 <- read.delim('113/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics113$subject <- 103
r_acoustics114 <- read.delim('114/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics114$subject <- 114
r_acoustics115 <- read.delim('115/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics115$subject <- 115
r_acoustics116 <- read.delim('116/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics116$subject <- 116
r_acoustics117 <- read.delim('117/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics117$subject <- 117
r_acoustics118 <- read.delim('118/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics118$subject <- 118
r_acoustics119 <- read.delim('119/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics119$subject <- 119
r_acoustics120 <- read.delim('120/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics120$subject <- 120
r_acoustics121 <- read.delim('../121/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics121$subject <- 121
r_acoustics122 <- read.delim('../122/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics122$subject <- 122
r_acoustics123 <- read.delim('../123/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics123$subject <- 123
r_acoustics124 <- read.delim('../124/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics124$subject <- 124
r_acoustics125 <- read.delim('../125/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics125$subject <- 125
r_acoustics126 <- read.delim('../126/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics126$subject <- 126
r_acoustics127 <- read.delim('../127/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics127$subject <- 127
r_acoustics128 <- read.delim('../128/Racoustic_data.txt', header = TRUE, sep = ',')
r_acoustics128$subject <- 128
rac<-rbind(#r_acoustics10
  r_acoustics102,
  r_acoustics103,
  r_acoustics104,
  r_acoustics105,
  # r_acoustics1
  r_acoustics108,
  r_acoustics109,
  r_acoustics110,
  r_acoustics112,
  r_acoustics113,
  r_acoustics114,
  r_acoustics115,
  r_acoustics116,
  r_acoustics117,
  r_acoustics118,
  r_acoustics119,
  r_acoustics120,
  r_acoustics121,
  r_acoustics122,
  r_acoustics123,
  r_acoustics124,
  r_acoustics125,
  r_acoustics126,
  r_acoustics127,
  r_acoustics128)



#f3sd101   <- sd(r_acoustics101$frmf3)
f3sd102    <- sd(r_acoustics102$frmf3)
f3sd103   <- sd(r_acoustics103$frmf3)
f3sd104    <- sd(r_acoustics104$frmf3)
f3sd105    <- sd(r_acoustics105$frmf3)
# f3sd106  <- sd(r_acoustics106$frmf3)
f3sd108    <- sd(r_acoustics108$frmf3)
f3sd109    <- sd(r_acoustics109$frmf3)
f3sd110    <- sd(r_acoustics110$frmf3)
f3sd112    <- sd(r_acoustics112$frmf3)
f3sd113    <- sd(r_acoustics113$frmf3)
f3sd114    <- sd(r_acoustics114$frmf3)
f3sd115    <- sd(r_acoustics115$frmf3)
f3sd116    <- sd(r_acoustics116$frmf3)
f3sd117    <- sd(r_acoustics117$frmf3)
f3sd118    <- sd(r_acoustics118$frmf3)
f3sd119    <- sd(r_acoustics119$frmf3)
f3sd120    <- sd(r_acoustics120$frmf3)
f3sd121    <- sd(r_acoustics121$frmf3)
f3sd122    <- sd(r_acoustics122$frmf3)
f3sd123    <- sd(r_acoustics123$frmf3)
f3sd124    <- sd(r_acoustics124$frmf3)
f3sd125    <- sd(r_acoustics125$frmf3)
f3sd126    <- sd(r_acoustics126$frmf3)
f3sd127    <- sd(r_acoustics127$frmf3)
f3sd128    <- sd(r_acoustics128$frmf3)


subject<-c(#101, 
        102, 
        #103,
        104, 
        105, 
        #106,
        108, 
        109, 
        110, 
        112, 
        113, 
        114, 
        115, 
        116, 
        117, 
        118, 
        119, 
        120, 
        121, 
        122, 
        123, 
        124, 
        125, 
        126, 
        127, 
        128)     

f3sd2019<-rbind(#f3sd101,
      f3sd102, 
      #f3sd103,
      f3sd104, 
      f3sd105, 
      #f3sd106,
      f3sd108, 
      f3sd109, 
      f3sd110, 
      f3sd112, 
      f3sd113, 
      f3sd114, 
      f3sd115, 
      f3sd116, 
      f3sd117, 
      f3sd118, 
      f3sd119, 
      f3sd120, 
      f3sd121, 
      f3sd122, 
      f3sd123, 
      f3sd124, 
      f3sd125, 
      f3sd126, 
      f3sd127, 
      f3sd128 )

racoustics = data.frame(subject, f3sd2019)
