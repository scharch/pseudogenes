#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


library(dplyr)
library(stringr)

rawdata <- read.delim(args, header = F)


## IGKL ##
names(rawdata)[5] <- "Bowtie_count"
names(rawdata)[8] <- "V_gene_assignment"
names(rawdata)[48] <- "V_J_frame"
rawdata$V_gene_assignment <- (str_sub(rawdata$V_gene_assignment, 1, str_length(rawdata$V_gene_assignment)-3))

## IGH ##
#names(rawdata)[5] <- "Bowtie_count"
#names(rawdata)[8] <- "V_gene_assignment"
#names(rawdata)[49] <- "V_J_frame"
#rawdata$V_gene_assignment <- (str_sub(rawdata$V_gene_assignment, 1, str_length(rawdata$V_gene_assignment)-3))

## Pseudogene Data ##
pseudogenes <- read.delim("/hpcdata/vrc/vrc1_data/douek_lab/projects/BCRSeq/2018213_gDNA/baldr/allpseudogenes.list", header = FALSE)

## Remove duplicates ##
rawdata <- distinct(rawdata, rawdata$V51, rawdata$V52, rawdata$V53, .keep_all = TRUE)


## Remove low counts ##
highestvalue <- max(rawdata$Bowtie_count)
threshold <- highestvalue/5
rawdata <- rawdata[(rawdata$Bowtie_count >threshold),]


## Filter by pseudogene or functional (prod or nonprod) ##
for (i in rownames(rawdata)){
  if (rawdata[i,"V_gene_assignment"] %in% pseudogenes$V1){
    rawdata[i, "Rearrangement_Type"] <- "Pseudogene"
  } else if (rawdata[i,"V_J_frame"] == "In-frame"){
    rawdata[i, "Rearrangement_Type"] <- "Functional_Productive"
  } else if (rawdata[i,"V_J_frame"] == "Out-of-frame"){
    rawdata[i, "Rearrangement_Type"] <- "Functional_Nonproductive"
  }
}

numberofProductive <- length(which(rawdata$Rearrangement_Type == "Functional_Productive"))
total <- length(rawdata$Rearrangement_Type)

if ("Pseudogene" %in% rawdata$Rearrangement_Type) {
  write(args, file = "IGKL_cellswithpseudo.txt", append = TRUE)
} else if ("Functional_Nonproductive" %in% rawdata$Rearrangement_Type && numberofProductive <=1){
  write(args, file = "IGKL_cellswithfunctional_nonproductive.txt", append = TRUE)
} else if ("Functional_Productive" %in% rawdata$Rearrangement_Type && numberofProductive > 1){
  write(args, file = "IGKL_cellstoignore.txt", append = TRUE)
} else if ("Functional_Productive" %in% rawdata$Rearrangement_Type && numberofProductive <= 1 && total == 1){
  write(args, file = "IGKL_cellswithfunctional_productive.txt", append = TRUE)
}
