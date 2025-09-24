#!/usr/bin/env Rscript
library(tidyverse)

# Hfib muscle
codes = c('1_1', '1_2', '2_1', '2_2', '3_1', '3_2', '4_1', '4_2', '5_1', '5_2', '6_1', 
          '6_2', '7_1', '8_1', '8_2', '9_1', '9_2', '10_1', '10_2', '11_1', '11_2', '12_1', '12_2', 
          '13_1', '13_2', '14_1', '14_2', '15_1', '15_2', '16_1', '16_2', '17_1', '18_1', '18_2')

colors = c(rep("navyblue",15), rep("seagreen",19))

pdf(paste0('alignment_fig_hfib_muscle.pdf'),height=12, width=15)

h1 = 2
h2 = 4
minY = 18
maxY = (h1+h1+h2) + (length(codes)-2) * (h1+h2) + minY
t=.8

plot(NULL, xlim = c(-500, 8000), ylim = c(-2,220), axes=F, ylab = '', xlab = '')

text(0, maxY+4, "H-Fibroin", pos=3)
text(3000, -2, "Length of Allele in Amino Acids", pos=4)
axis(1, at = seq(0,7500, by=500), labels = seq(0,7500, by=500), line=-5)

for (i in 1:(length(codes)-1)) {
  
  code1 = codes[i]
  code2 = codes[i+1]
  
  coords <- read.csv(paste0('hfib_muscle/', code1, 'vs', code2, 'coords.csv'), header=TRUE)
  
  rect(0,maxY-h1, tail(coords, n=1)$X2, maxY, col = colors[i], border='lightgray', lwd=.1)
  rect(0,maxY-(h1+h1+h2), tail(coords, n=1)$X4, maxY-(h1+h2), col = colors[i+1], border='lightgray', lwd=.1)
  
  text(-275, maxY+h1+t, code1, pos = 1, cex = t)
  text(-275, maxY-h2+t, code2, pos = 1, cex = t)
  
  for(i in 1:nrow(coords)){
    
    polygon(x = c(coords$X1[i], coords$X2[i], coords$X4[i], coords$X3[i]),
            y = c(maxY-h1,maxY-h1,maxY-(h1+h2),maxY-(h1+h2)),
            col = adjustcolor('black', alpha.f = 0.35), border = 'lightgray',lwd = .1)
  }
  maxY = maxY-(h1+h2)
}

dev.off()

# Hfib mafft
codes = c('1_1', '1_2', '2_1', '2_2', '3_1', '3_2', '4_1', '4_2', '5_1', '5_2', '6_1', 
          '6_2', '7_1', '8_1', '8_2', '9_1', '9_2', '10_1', '10_2', '11_1', '11_2', '12_1', '12_2', 
          '13_1', '13_2', '14_1', '14_2', '15_1', '15_2', '16_1', '16_2', '17_1', '18_1', '18_2')

colors = c(rep("navyblue",15), rep("seagreen",19))

pdf(paste0('alignment_fig_hfib_mafft.pdf'),height=12, width=15)

h1 = 2
h2 = 4
minY = 18
maxY = (h1+h1+h2) + (length(codes)-2) * (h1+h2) + minY
t=.8

plot(NULL, xlim = c(-500, 8000), ylim = c(-2,220), axes=F, ylab = '', xlab = '')

text(0, maxY+4, "H-Fibroin", pos=3)
text(3000, -2, "Length of Allele in Amino Acids", pos=4)
axis(1, at = seq(0,7500, by=500), labels = seq(0,7500, by=500), line=-5)

for (i in 1:(length(codes)-1)) {
  
  code1 = codes[i]
  code2 = codes[i+1]
  
  coords <- read.csv(paste0('hfib_mafft/', code1, 'vs', code2, 'coords.csv'), header=TRUE)
  
  rect(0,maxY-h1, tail(coords, n=1)$X2, maxY, col = colors[i], border='lightgray', lwd=.1)
  rect(0,maxY-(h1+h1+h2), tail(coords, n=1)$X4, maxY-(h1+h2), col = colors[i+1], border='lightgray', lwd=.1)
  
  text(-275, maxY+h1+t, code1, pos = 1, cex = t)
  text(-275, maxY-h2+t, code2, pos = 1, cex = t)
  
  for(i in 1:nrow(coords)){
    
    polygon(x = c(coords$X1[i], coords$X2[i], coords$X4[i], coords$X3[i]),
            y = c(maxY-h1,maxY-h1,maxY-(h1+h2),maxY-(h1+h2)),
            col = adjustcolor('black', alpha.f = 0.35), border = 'lightgray',lwd = .1)
  }
  maxY = maxY-(h1+h2)
}

dev.off()

# Hfib prank
codes = c('1_1', '1_2', '2_1', '2_2', '3_1', '3_2', '4_1', '4_2', '5_1', '5_2', '6_1', 
          '6_2', '7_1', '8_1', '8_2', '9_1', '9_2', '10_1', '10_2', '11_1', '11_2', '12_1', '12_2', 
          '13_1', '13_2', '14_1', '14_2', '15_1', '15_2', '16_1', '16_2', '17_1', '18_1', '18_2')

colors = c(rep("navyblue",15), rep("seagreen",19))

pdf(paste0('alignment_fig_hfib_prank.pdf'),height=12, width=15)

h1 = 2
h2 = 4
minY = 18
maxY = (h1+h1+h2) + (length(codes)-2) * (h1+h2) + minY
t=.8

plot(NULL, xlim = c(-500, 8000), ylim = c(-2,220), axes=F, ylab = '', xlab = '')

text(0, maxY+4, "H-Fibroin", pos=3)
text(3000, -2, "Length of Allele in Amino Acids", pos=4)
axis(1, at = seq(0,7500, by=500), labels = seq(0,7500, by=500), line=-5)

for (i in 1:(length(codes)-1)) {
  
  code1 = codes[i]
  code2 = codes[i+1]
  
  coords <- read.csv(paste0('hfib_prank/', code1, 'vs', code2, 'coords.csv'), header=TRUE)
  
  rect(0,maxY-h1, tail(coords, n=1)$X2, maxY, col = colors[i], border='lightgray', lwd=.1)
  rect(0,maxY-(h1+h1+h2), tail(coords, n=1)$X4, maxY-(h1+h2), col = colors[i+1], border='lightgray', lwd=.1)
  
  text(-275, maxY+h1+t, code1, pos = 1, cex = t)
  text(-275, maxY-h2+t, code2, pos = 1, cex = t)
  
  for(i in 1:nrow(coords)){
    
    polygon(x = c(coords$X1[i], coords$X2[i], coords$X4[i], coords$X3[i]),
            y = c(maxY-h1,maxY-h1,maxY-(h1+h2),maxY-(h1+h2)),
            col = adjustcolor('black', alpha.f = 0.35), border = 'lightgray',lwd = .1)
  }
  maxY = maxY-(h1+h2)
}

dev.off()

# Lfib
codes = c('1_1', '1_2', '2_1', '2_2', '3_1', '3_2', '4_1', '4_2', '5_1', '5_2', '6_2', 
          '7_1', '7_2', '8_1', '8_2', '9_1', '9_2', '10_1', '10_2', '11_1', '11_2', '12_1', '12_2', 
          '13_1', '13_2', '14_2', '15_1', '15_2', '16_1', '16_2', '17_1', '17_2', '18_1', '18_2')

colors = c(rep("navyblue",15), rep("seagreen",19))

pdf(paste0('alignment_fig_lfib.pdf'),height=12, width=15)

h1 = 2
h2 = 4
minY = 18
maxY = (h1+h1+h2) + (length(codes)-2) * (h1+h2) + minY
t=.8

plot(NULL, xlim = c(-20, 295), ylim = c(-2,220), axes=F, ylab = '', xlab = '')

text(0, maxY+4, "L-Fibroin", pos=3)
text(100, -2, "Length of Allele in Amino Acids", pos=4)
axis(1, at = seq(0,275, by=25), labels = seq(0,275, by=25), line=-5)

for (i in 1:(length(codes)-1)) {
  
  code1 = codes[i]
  code2 = codes[i+1]
  
  coords <- read.csv(paste0('lfib/', code1, 'vs', code2, 'coords.csv'), header=TRUE)
  
  rect(0,maxY-h1, tail(coords, n=1)$X2, maxY, col = colors[i], border='lightgray', lwd=.1)
  rect(0,maxY-(h1+h1+h2), tail(coords, n=1)$X4, maxY-(h1+h2), col = colors[i+1], border='lightgray', lwd=.1)
  
  text(-11, maxY+h1+t, code1, pos = 1, cex = t)
  text(-11, maxY-h2+t, code2, pos = 1, cex = t)
  
  for(i in 1:nrow(coords)){
    
    polygon(x = c(coords$X1[i], coords$X2[i], coords$X4[i], coords$X3[i]),
            y = c(maxY-h1,maxY-h1,maxY-(h1+h2),maxY-(h1+h2)),
            col = adjustcolor('black', alpha.f = 0.35), border = 'lightgray',lwd = .1)
  }
  maxY = maxY-(h1+h2)
}

dev.off()

# PEVK
codes = c('1_1', '1_2', '2_1', '2_2', '3_1', '3_2', '4_1', '4_2', '5_1', '5_2', '6_1', '6_2', '7_1', 
          '7_2', '8_1', '8_2', '9_1', '9_2', '10_1', '10_2', '11_1', '11_2', '12_1', '12_2', '13_1', 
          '13_2', '14_1', '14_2', '15_1', '15_2', '16_1', '16_2', '17_1', '17_2', '18_1', '18_2')

colors = c(rep("navyblue",16), rep("seagreen",20))

pdf(paste0('alignment_fig_pevk.pdf'),height=12, width=15)

h1 = 2
h2 = 4
minY = 20
maxY = (h1+h1+h2) + (length(codes)-2) * (h1+h2) + minY
t=.8

plot(NULL, xlim = c(-30, 530), ylim = c(-2,235), axes=F, ylab = '', xlab = '')

text(0, maxY+4, "PEVK", pos=3)
text(190, -2, "Length of Allele in Amino Acids", pos=4)
axis(1, at = seq(0,500, by=50), labels = seq(0,500, by=50), line=-5)

for (i in 1:(length(codes)-1)) {
  
  code1 = codes[i]
  code2 = codes[i+1]
  
  coords <- read.csv(paste0('pevk/', code1, 'vs', code2, 'coords.csv'), header=TRUE)
  
  rect(0,maxY-h1, tail(coords, n=1)$X2, maxY, col = colors[i], border='lightgray', lwd=.1)
  rect(0,maxY-(h1+h1+h2), tail(coords, n=1)$X4, maxY-(h1+h2), col = colors[i+1], border='lightgray', lwd=.1)
  
  text(-20, maxY+h1+t, code1, pos = 1, cex = t)
  text(-20, maxY-h2+t, code2, pos = 1, cex = t)
  
  for(i in 1:nrow(coords)){
    
    polygon(x = c(coords$X1[i], coords$X2[i], coords$X4[i], coords$X3[i]),
            y = c(maxY-h1,maxY-h1,maxY-(h1+h2),maxY-(h1+h2)),
            col = adjustcolor('black', alpha.f = 0.35), border = 'lightgray',lwd = .1)
  }
  maxY = maxY-(h1+h2)
}

dev.off()

# Resilin
codes = c('1_1', '1_2', '2_1', '2_2', '3_1', '3_2', '4_1', '4_2', '5_1', '5_2', '6_1', '6_2', '7_1', 
          '7_2', '8_1', '8_2', '9_1', '9_2', '10_1', '10_2', '11_1', '11_2', '12_1', '12_2', '13_1', 
          '13_2', '14_1', '14_2', '15_1', '15_2', '16_1', '16_2', '17_1', '17_2', '18_1', '18_2')

colors = c(rep("navyblue",16), rep("seagreen",20))

for (j in c("A","B")){
  pdf(paste0('alignment_fig_resilin', j, '.pdf'), height=12, width=15)
  
  h1 = 2
  h2 = 4
  minY = 20
  maxY = (h1+h1+h2) + (length(codes)-2) * (h1+h2) + minY
  t=.8
  
  plot(NULL, xlim = c(-50, 750), ylim = c(-2,235), axes=F, ylab = '', xlab = '')
  
  text(0, maxY+4, paste0("Resilin ",j), pos=3)
  text(275, -2, "Length of Allele in Amino Acids", pos=2)
  axis(1, at = seq(0,700, by=50), labels = seq(0,700, by=50), line=-5)
  
  for (i in 1:(length(codes)-1)) {
    
    code1 = paste0(codes[i],j)
    code2 = paste0(codes[i+1],j)
    
    coords <- read.csv(paste0('resilin/', code1, 'vs', code2, 'coords.csv'), header=TRUE)
    
    rect(0,maxY-h1, tail(coords, n=1)$X2, maxY, col = colors[i], border='lightgray', lwd=.1)
    rect(0,maxY-(h1+h1+h2), tail(coords, n=1)$X4, maxY-(h1+h2), col = colors[i+1], border='lightgray', lwd=.1)
    
    text(-20, maxY+h1+t, code1, pos = 1, cex = t)
    text(-20, maxY-h2+t, code2, pos = 1, cex = t)
    
    for(i in 1:nrow(coords)){
      
      polygon(x = c(coords$X1[i], coords$X2[i], coords$X4[i], coords$X3[i]),
              y = c(maxY-h1,maxY-h1,maxY-(h1+h2),maxY-(h1+h2)),
              col = adjustcolor('black', alpha.f = 0.35), border = 'lightgray',lwd = .1)
    }
    maxY = maxY-(h1+h2)
  }
  dev.off()
}
