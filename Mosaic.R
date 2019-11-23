## I used the terminal in BioLinux to change the extensions from .fna to .fasta
#       for f in *.fna
#       do
#       [ -f "$f" ] && mv "$f" "${f%fna}fasta"
#       done

library("seqinr")
library("msa")
library("tidyverse")

### Selection for genes ###################################################################################################################################
# If on my own computer:
setwd("C:/Users/Kim/Documents/School/2019_3Fall/Biology_498/Mosaic_Mixta/")
fastaFiles <- as.data.frame(list.files(path = "C:/Users/Kim/Documents/School/2019_3Fall/Biology_498/Mosaic_Mixta/Genes", pattern = ".fasta"))
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("C:/Users/Kim/Documents/School/2019_3Fall/Biology_498/Mosaic_Mixta/Genes", 
                              fastaFiles$File_name, sep = "/")            # Creates a file pathway for each gene

# If in the office:
setwd("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/")
fastaFiles <- as.data.frame(list.files(path = "3Homologues_10/"))              # Makes a dataframe where the first column is a list of fasta gene files
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/3Homologues_10", 
                              fastaFiles$File_name, sep = "/")            # Creates a file pathway for each gene

####
ten_seq <- function(gene) {                                               # Function that checks if there are 10 sequences in the gene file
  if (nrow(gene) == 20) {                                                 # Checks if there's 20 rows because odd rows are names and even rows are seq.
    print("You may continue")
    twenty <- "Yes"                                                       # Confirmation variable for when the function is called
  } else {
    print("Not this gene")
  }
}
#
gene_length <- function(gene_file) {                                      # Function that checks for similar gene lengths
  count = 0
  for(row in 1:nrow(gene_file)) {
    if(nchar(gene_file$sequences)[row] >= max(nchar(gene_file$sequences), 
                                              na.rm = TRUE) * 0.9) {      # All genes must have a length of at least 90% of the longest gene
      count <- count + 1                                                  # Confirmation variable for when the function is called, must be 10
    }
  }
  count
}

for(row in 1:nrow(fastaFiles)) {
  path <- fastaFiles$Path_name[row]
  gene_file <- read.table(file = path, header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE)                       # Reads in the gene file according to the pathway, separation is newline
  
  twenty <- ten_seq(gene_file)
  if(twenty == "Yes") {                                                   # If 10 sequences, then continue
    print("    Yes")
    gene_file <- data.frame(sequences = gene_file$V1[1:10 * 2], 
                            species = gene_file$V1[1:10 * 2 - 1], 
                            stringsAsFactors = FALSE)                     # Dataframe where first column are sequences and second are corresponing names
    
    gene_file$species <- c("ID:IHMAAOCK_01692__Tatumella saanichensis__NML_06-3099", "ID:JADOELLM_02794__Citrobacter freundii__NCTC_9750",
                           "ID:OJPMODEH_01576__Enterobacter cloacae_subsp_cloacae__ATCC 13047", "ID:DGOHCAMD_01524__Erwinia amylovora__CFBP_1232",
                           "ID:LMPCFIAF_01940__Erwinia tasmaniensis__ET1-99", "ID:MLHGDOMH_01877__Mixta calida__DSM_22759",
                           "ID:MHMNNPCM_01846__Mixta gaviniae__DSM_22758", "ID:MEFHALAL_00346__Pantoea agglomerans__NBRC_102470",
                           "ID:KOGPCAHI_00156__Pantoea septica__LMG_5345", "ID:IBJKPIAN_01672__Tatumella ptyseos__NCTC11468")
    # Renames the species names so they aren't ridiculously long
    
    count <- gene_length(gene_file)
    if(count == 10) {                                                     # If all 10 are of similar lengths, then continue
      print("        YES!!")
      write.fasta(sequences = as.list(gene_file$sequences),
                  names = gene_file$species,
                  file.out = paste("4Organize/", fastaFiles$File_name[row],
                                   sep = ''),
                  open = "w", nbchar = 10000, as.string = TRUE)           # Creates a fasta file for each gene, will continue on to alignment
    } else {
      print("        Nope")
    }
  }
}
rm(gene_file, count, path, row, twenty)

### Aligning genes ########################################################################################################################################
# Using ClustalW in R
fastaFilesOrg <- as.data.frame(list.files(path = "C:/Users/Kim/Documents/School/2019_3Fall/Biology_498/Mosaic_Mixta/4Organize/", pattern = ".fasta"))
colnames(fastaFilesOrg) <- "File_name"                                       # Changes the column name
fastaFilesOrg$Path_name <- paste("C:/Users/Kim/Documents/School/2019_3Fall/Biology_498/Mosaic_Mixta/4Organize", 
                                 fastaFilesOrg$File_name, sep = "/")            # Creates a file pathway for each gene

for(row in 1:nrow(fastaFilesOrg)) {
  path <- fastaFilesOrg$Path_name[row]
  
  gene_file <- readDNAStringSet(path)
  
  align <- msa::msaClustalW(inputSeqs = gene_file, maxiters = 100, type = "dna", order = "input")
  alignConv <- msaConvert(align, type = "seqinr::alignment")
  
  alignFA <- base::as.data.frame(matrix(ncol = 2, nrow = 10))
  colnames(alignFA) <- c("sequences", "species")
  alignFA$sequences <- alignConv$seq
  alignFA$species <- alignConv$nam
  
  write.fasta(sequences = as.list(alignFA$sequences),
              names = alignFA$species,
              file.out = paste("5Aligned/", fastaFilesOrg$File_name[row],
                               sep = ''),
              open = "w", nbchar = 10000, as.string = TRUE)
}
beep(8)
rm(align, alignConv, alignFA, gene_file, path, row)

### Distance matrices #####################################################################################################################################






#
### reading in the distance matrix #########
library("tidyverse")
library("readxl")
test <- read.csv(file = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/3Homologues_10/Test_csv", header = FALSE)
