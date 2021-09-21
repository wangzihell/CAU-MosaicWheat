#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-g", "--geoDB"), dest = "geoDB", default = "/data2/rawdata2/hapmap/script/sample_location_tetraintro_WEC_201221.txt",
              help = "geographic data file, accession location lon lat [default: sample_location.txt]"),
  make_option(c("-i", "--sampleF"), dest = "sampleF", default = "clade_map_input.txt",
              help = "sample list, one accessions per line. [default: id_group.txt]"),
  make_option(c("-m", "--mergeD"), dest = "mergeD", default = 3,
              help = "merge samples within distance [5]"),
  make_option(c("-r", "--removeLoc"), dest = "removeLoc", default = 3,
              help = "merge samples within distance [5]"),  
  make_option(c("-o", "--outfile"), dest = "outfile", default = "test.txt",
              help = "output prefix [default: output]")
)

parser <- OptionParser(usage = "Rscript treeWholeNJ.R [options] Tip group files",
                       description = "change the location of each accession to mergeed location.",
                       option_list = option_list)

arguments <- parse_args(parser)

geoDB <- arguments$geoDB
sampleF <- arguments$sampleF
mergeD <- as.numeric(arguments$mergeD)
removeLoc <- as.numeric(arguments$removeLoc)
outfile <- arguments$outfile

# merge locs between which dist less than disthr
mergeloc <- function(locdf, disthr = 0.1){
  
  # merge same location, ggmap use location as group
  if(disthr < 0.1){disthr = 0.1}
  
  newlocdf <- locdf[FALSE,]
  
  for (i in seq_len(nrow(locdf))){
    for (j in seq_len(nrow(locdf))){
      newlocdf[i, j] <- as.numeric(dist(rbind(locdf[i,c(3,4)], locdf[j,c(3,4)])))
    }
  }
  dis_matrix <- cutree(hclust(as.dist(newlocdf)), h = disthr)
      
  return(dis_matrix)
}

# sample geographical location
geoDB <- read.table(geoDB, header = TRUE, stringsAsFactors = F, sep="\t")
geoDB$lon <- as.numeric(geoDB$lon)
geoDB$lat <- as.numeric(geoDB$lat)

# genotype data of one site eg. chr1A.1:1000 (a snp site)
DF1 <- read.table(sampleF, header = T, stringsAsFactors = F)
DF1 <- DF1[,1, drop=F]

DF1 <- inner_join(geoDB, DF1, by = 'Accession')
if (!all(complete.cases(DF1))) {
  print("incomplete cases droped")
  DF1 <- DF1[complete.cases(DF1),]
}

# disthr should change with size of map to avoid points overlaping
DF2 <- mergeloc(DF1, disthr = mergeD)
DF1[,5] <- DF2

DF3 <- DF1 %>% group_by(V5) %>% mutate(lon1 = mean(lon), lat1 = mean(lat)) %>% ungroup()
DF3 <- DF3[,c(1,2,6,7,5)]

DF3 <- DF3[order(DF3$V5),]
if(is.null(removeLoc)){
  for (i in unique(DF3$V5)){
    tmp_DF3 <- DF3[DF3$V5 == i,]
    DF3[DF3$V5 == i,2] <- rep(paste0(tmp_DF3[1,2, drop=T],"_",i), nrow(tmp_DF3))
  }
} else {
  row_remove <- c()
  for (i in unique(DF3$V5)){
    tmp_DF3 <- DF3[DF3$V5 == i,]
    if (nrow(tmp_DF3) <= removeLoc) row_remove <- c(row_remove, i)
    DF3[DF3$V5 == i,2] <- rep(paste0(tmp_DF3[1,2, drop=T],"_",i), nrow(tmp_DF3))
  }
  DF3 <- DF3[! DF3$V5 %in% row_remove,]
}

colnames(DF3) <- c("Accession", "location", "lon", "lat", "Group")
write.table(DF3, outfile, quote = F, sep = "\t", col.names = T, row.names = F)
