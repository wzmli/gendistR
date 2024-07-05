library(dplyr)
library(data.table)

consensus <- fread("./data/nextclade_consensus_2024-05-17.tsv")
# consensus sequences taken from https://github.com/corneliusroemer/pango-sequences/tree/main on May 17th 2024.

refs <- consensus %>%
  filter(Nextclade_pango %in% c("B.1.1.7", "B.1.351", "P.1", "B.1.617.2", "B.1.427", "B.1.525", "B.1.526", "B.1.617.1", "C.37", "B.1.621", "BA.1", "BA.2", "BA.2.75", "BA.4", "BA.5", "BQ.1", "XBB", "XBB.1.5", "EG.5", "BA.2.86", "JN.1.9.2", "JN.1.11"))

#B.1.1.7 alpha
#B.1.351 beta
#P.1 gamma
#B.1.617.2 delta # note delta is a combination of many things, this is the first
#B.1.427 #note epsilon also includes B.1.429
#B.1.525 # eta
#B.1.526 # iota
#B.1.617.1 # kappa
#C.37 #lambda
#B.1.621 # mu
#BA.2.75
#EG.5 is XBB.1.9

#no consensus for JN.1.9.2 yet.

unique(refs$qc.overallStatus)# all have good qc

save(refs, file = "./data/ref_consensus.RData")
