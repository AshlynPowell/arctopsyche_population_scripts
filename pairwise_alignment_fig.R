library(tidyverse)

codes = c('ArcGran', 'LepLine', 'AtoDavi', 'HimTibe', 'EubRegi', 'GlyPell', 'HesMag',
          'LimLuna', 'LimMarm', 'LimRhom', 'AthCine', 'CerDiss', 'MysLong', 'OdoAlbi')

# code = codes[1]

for (code in codes) {
  species_lookup = read.csv('species_lookup.csv', header = TRUE)
  species = pull(filter(species_lookup, Code == code), Name)
  
  coords <- read.csv(paste0(code, '_both_alleles.csv'), header = TRUE) 
  
  pdf(paste0(code, '_both_alleles.pdf'), height = 4, width = 9) 
  
  plot(NULL, xlim = c(0, 10000), ylim = c(-3,3), axes=F, ylab = '', xlab = '')
  
  rect(0,1.4,tail(coords, n=1)$X2,1.7, col = 'navyblue', border = 'lightgray',lwd = .1)
  rect(0,.2,tail(coords, n=1)$X4,.5, col = 'seagreen', border = 'lightgray',lwd = .1)
  
  axis(1, at = seq(0,10000, by=2000), labels = seq(0,10000, by=2000), line = -4, cex.axis = 2)
  
  text(5000, -2.7, "Allele Length (Amino Acids)", cex = 2, adj = 0.5)
  text(100, 2.5, species, cex = 2, adj = 0)
  
  for(i in 1:nrow(coords)){
    
    polygon(x = c(coords$X1[i], coords$X2[i], coords$X4[i], coords$X3[i]), 
            y = c(1.4,1.4,.5,.5), 
            col = adjustcolor('black', alpha.f = 0.35), border = 'lightgray',lwd = .1)
  }
  
  dev.off()
  
}

# ?axis
