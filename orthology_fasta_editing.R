library(tidyverse)

###--- FUNCTIONS

## NAMES to species
adeq <- function (fasta) {
  
  raw_fasta <- readLines(fasta)
  index <- str_which(raw_fasta, ">")
  
  raw_fasta [index] <- 
    str_remove(raw_fasta [index], "PREDICTED: ")
  raw_fasta [index] <- 
    paste0 (">",
           map(str_split (raw_fasta [index], " "), 2), 
           "_",
           map(str_split (raw_fasta [index], " "), 3))
  
  return (raw_fasta)
  
}

## COMPARE of species names
fasta_compare <- function (fastas) {
  len = length(fastas)
  
  i <- 1
  sp_fasta <- list ()
  while (i <= len) {
    raw_fasta <- adeq (fastas[i])
    sp_fasta [[i]] <- raw_fasta[str_which(raw_fasta, "_")]
    i <- i + 1
  }
  
  i <- 1
  sp <- sp_fasta [[i]]
  while (i < len) {
    sp <- intersect (sp, sp_fasta[[i+1]])
    i <- i + 1
  }
  
  return (sp)

}

## SELECT OF GROUP/TAXA
sp_to_groups <- function (sp, grp, group_name) {
  
  groups <- readLines (grp)
  groups <- str_split(groups, "=")
  groups <- data.frame(sp = unlist(map (groups, 1)), 
                       group = unlist(map (groups, 2)))
  groups$sp <- paste0 (">", groups$sp)
  
  groups <- groups[groups$group == group_name,]
  
  sp <- intersect (sp, groups$sp)
  return (sp)
  
}

## ALL FILTRATIONS
fasta_filtration <- function (fastas, output, grp = "NA", group_name = "NA") {
  len = length(fastas)
  
  sp <- fasta_compare(fastas)
  
  if (group_name != "NA") {
    sp <- sp_to_groups (sp, grp, group_name)
  }
  
  i <- 1
  while (i <= len) {
    raw_fasta <- adeq (fastas[i])
    
    index <- str_which(raw_fasta, "_")
    sp_fasta <- raw_fasta[index]
    compare <- is.element (sp_fasta, sp)
    #delete of filtered data
    index_first <- index[!compare]
    index_second <- index[c(F, !compare[-length(compare)])]
    #deduplication
    indexF1 <- index_first[!(index_first %in% index_second)]
    indexF2 <- index_second[!(index_second %in% index_first)] - 1
    sorter <- sort(c(indexF1, indexF2))
    
    j <- 1
    MEAN <- T
    trusher <- c()
    while (j < length (raw_fasta)) {
      if (j %in% sorter) {
        MEAN <- !MEAN
      }
      trusher <- c (trusher, MEAN)
      j <- j + 1
    }
    
    raw_fasta <- raw_fasta [trusher]
    
    if (group_name != "NA") {
      filename <- paste0(group_name, "_", output[i], ".fasta")
    } else {
      filename <- paste0(output[i], ".fasta")
    }
    
    
    
    writeLines(raw_fasta, filename)
    
    print (paste0(output[i], ".fasta was written..."))
    
    i <- i + 1
  }
  print ("Done! Petrov loh)))")
}



###---
setwd ("~/CCR5_speedrun")
fastas <- c ("CCR5_refseq_transcript.fasta",
            "CD4_refseq_transcript.fasta",
            "CXCR4_refseq_transcript.fasta")
grp <- "species.grp"
group_name <- "Primates"
output <- c ("CCR5", "CD4", "CXCR4")

fasta_filtration (fastas, output)
fasta_filtration (fastas, output, grp, group_name)
