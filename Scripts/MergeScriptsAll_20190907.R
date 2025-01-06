##########################
##### Merge Tables  #######
###########################


library(dplyr)
library(tidyverse)

# A <- read.csv("data/suba4_withoutnames.csv", header = T, stringsAsFactors = F)
# 
# suba4_Iso_ED <- A %>% distinct(AGI,location, .keep_all = TRUE)
# 
# write.csv(suba4_Iso_ED, "suba4_Iso_ED.csv") # Then removed acending number column

# Merging the files > Final

c <- read.csv("data/MapBIN_ED.csv", header = T, stringsAsFactors = F)
d <- read.csv("data/suba4_Iso_ED.csv", header = T, stringsAsFactors = F)
e <- read.csv("data/APR18_0.05_LFC1.5_190907_ED.csv", header = T, stringsAsFactors = F)


MapBIN_SUBA_LFQ <- left_join(d, e, by = "AGI")

MapBIN_SUBA_LFQ$AGI_ED <- str_sub(MapBIN_SUBA_LFQ$AGI, end = -3)
MapBIN_SUBA_LFQ <- left_join(MapBIN_SUBA_LFQ, c, by = "AGI_ED")
MapBIN_SUBA_LFQ <- filter(MapBIN_SUBA_LFQ, !is.na(ID))


write.csv(MapBIN_SUBA_LFQ, "18APRDEP_Results_190907.csv")

##########################
## Merge with Araport11 ##
##########################

d <- read.csv("data/suba4_Iso_ED.csv", header = T, stringsAsFactors = F)
e <- read.csv("data/APR19_0.05_LFC1.5_190903_ED.csv", header = T, stringsAsFactors = F)
f <- read.csv("data/Araport11_ED.csv", header = T, stringsAsFactors = F)

Araport11_SUBA_LFQ <- left_join(e, f, by = "AGI") %>% left_join(d, by = "AGI")


write.csv(Araport11_SUBA_LFQ, "19APR_Araport11_SUBA_LFQ_20190903.csv")





