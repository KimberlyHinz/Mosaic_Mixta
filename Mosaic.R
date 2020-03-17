### Information ###########################################################################################################################################
# Author: Kim Hinz
# Date of study: 2019-09-05 -- 2020-04-01
# Purpose: Phylogenetic analyses of Mixta genes to determine origination.
# Name of study: Mosaic Mixta


# The following code is for the phylogenetic study of two Mixta species (bacteria). Depending on the genes, model, and statistical method chosen for the 
# phylogenetic tree, the Mixta species group with different genera (see Palmer et al. 2018 and Rezzonico et al. 2016 for an example). Primarily, Mixta 
# appears to be a close relative to Pantoea with some leaning towards Erwinia. The purpose of my research is to determine why this might be the case by 
# performing distance matrix analyses for the homologous genes between two Mixta species, two Pantoea species, two Erwinia species, two Tatumella species, 
# one Citrobacter species, and one Enterobacter species. The Citrobacter and Enterobacter species form the outgroup.


# Publicly available genomes of the type strains of the species included in this study were retrieved from NCBI.
# # Mixta calida DSM_22759 ----------------------------> complete genome
# # Mixta gaviniae DSM_22758 --------------------------> complete genome
# # Pantoea agglomerans NBRC_102470 -------------------> contigs
# # Pantoea septica LMG_5345 --------------------------> contigs
# # Erwinia amylovora CFBP_1232 -----------------------> contigs
# # Erwinia tasmaniensis ET1/99 -----------------------> complete genome
# # Tatumella ptyseos NCTC_11468 ----------------------> complete genome
# # Tatumella saanichensis NML_06-3099 ----------------> contigs
# # Citrobacter freundii NCTC_9750 --------------------> complete genome
# # Enterobacter cloacae subsp cloacae ATCC_13047 -----> complete genome


# These genomes were annotated using PROKKA version 1.14.1 and core genes were extracted using the GET_HOMOLOGUES software package with the bidirectional
# best-hit search algorithm using default parameters.

# The output file type of GET_HOMOLOGUES is .fna. Therefore, I used the terminal in BioLinux to change the extensions from .fna to .fasta
#       for f in *.fna
#       do
#       [ -f "$f" ] && mv "$f" "${f%fna}fasta"
#       done

### Packages ##############################################################################################################################################
library("seqinr")

library("plyr")
library("msa")
library("beepr")
library("dplyr")

library("ape")
library("adegenet")
library("Rfast")

library("plyr")
library("tidyr")
library("ggplot2")
theme_set(theme_bw())

setwd("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/")
setwd("D:/")
#
### Selection for genes ###################################################################################################################################
# This portion of the code acts as filter; it passes forward files that have full sequences (genes aren't split) and that don't have truncated sequences
# (genes must be at least 90% of the length of the longest gene in each file). For example, if the longest sequence is 1000 bp, then the rest of the 
# sequences in that file must be at least 900 bp. If at least one is shorter than 900 bp, the whole file is excluded.

fastaFiles <- as.data.frame(list.files(path = "3Homologues_10/",
                                       pattern = ".fasta"))               # Makes a dataframe where the first column is a list of fasta gene files
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("3Homologues_10/", 
                              fastaFiles$File_name, 
                              sep = "")                                   # Creates a file pathway for each gene

ten_seq <- function(gene) {                                               # Function that checks if there are 10 sequences in the gene file
  if(nrow(gene) == 20) {                                                  # Checks if there's 20 rows because odd rows are names and even rows are seq.
    print("You may continue")
    twenty <- "Yes"                                                       # Confirmation variable for when the function is called
  } else {
    print("Not this gene")
    twenty <- "No"
  }
}

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

for(row in 1:nrow(fastaFiles)) {                                          # Lets pass genes that meet requirements
  path <- fastaFiles$Path_name[row]
  gene_file <- read.table(file = path, 
                          header = FALSE, 
                          sep = "\n", 
                          stringsAsFactors = FALSE)                       # Reads in the gene file according to the pathway, separation is newline
  
  twenty <- ten_seq(gene_file)
  if(twenty == "Yes") {                                                   # If 10 sequences, then continue
    print("    Yes")
    gene_file <- data.frame(sequences = gene_file$V1[1:10 * 2], 
                            species = gene_file$V1[1:10 * 2 - 1], 
                            stringsAsFactors = FALSE)                     # Dataframe where first column are sequences and second are corresponing names
    
    gene_file$species <- c("Tatumella saanichensis__NML_06-3099", "Citrobacter freundii__NCTC_9750", "Enterobacter cloacae_subsp_cloacae__ATCC 13047", 
                           "Erwinia amylovora__CFBP_1232", "Erwinia tasmaniensis__ET1-99", "Mixta calida__DSM_22759", "Mixta gaviniae__DSM_22758", 
                           "Pantoea agglomerans__NBRC_102470", "Pantoea septica__LMG_5345", 
                           "Tatumella ptyseos__NCTC_11468")               # Renames the species names so they aren't ridiculously long
    
    count <- gene_length(gene_file)
    if(count == 10) {                                                     # If all 10 are of similar lengths, then continue
      print("        YES!!")
      write.fasta(sequences = as.list(gene_file$sequences),
                  names = gene_file$species,
                  file.out = paste("4Organize/", 
                                   fastaFiles$File_name[row], sep = ''),
                  open = "w", nbchar = 10000, as.string = TRUE)           # Creates a fasta file for each gene, will continue on to alignment
    } else {
      print("        Nope")
    }
  }
}
rm(gene_file, count, path, row, twenty)

### Aligning genes ########################################################################################################################################
# This section aligns the genes that passed the filter using ClustalW through the R package msa. The parameters are 100 maximum iterations (default is 16)
# and default parameters. Then, the genes are written into a new fasta file.

fastaFilesOrg <- as.data.frame(list.files(path = "4Organize/", 
                                          pattern = ".fasta"))            # Makes a dataframe listing the fasta files in the folder
colnames(fastaFilesOrg) <- "File_name"                                    # Changes the column name
fastaFilesOrg$Path_name <- paste("4Organize/", 
                                 fastaFilesOrg$File_name, 
                                 sep = "")                                # Creates a file pathway for each gene

align_gene <- function(gene_file) {                                       # Aligns the sequences in the file and converts to a dataframe
  align <- msa::msaClustalW(inputSeqs = gene_file, 
                            maxiters = 100, 
                            type = "dna",
                            order = "input")                              # Aligns the genes
  
  alignConv <- msaConvert(align, type = "seqinr::alignment")              # Converts the aligned genes into a readable form
  
  alignFA <- base::as.data.frame(matrix(ncol = 0, nrow = 10))             # Create the necessary format for write.fasta()
  alignFA <- mutate(alignFA, 
                    sequences = alignConv$seq,
                    species = alignConv$nam)                              # Copies two columns from alignConv to alignFA
}

for(row in 1:nrow(fastaFilesOrg)) {                                       # Aligns genes that passed the filter
  path <- fastaFilesOrg$Path_name[row]
  gene_file <- readDNAStringSet(path)                                     # Reads in the genes (this way ensures the names are kept)
  
  alignFA <- align_gene(gene_file)                                        # Aligns the sequences and converts to a readable dataframe
  
  write.fasta(sequences = as.list(alignFA$sequences),
              names = alignFA$species,
              file.out = paste("5Aligned/", fastaFilesOrg$File_name[row],
                               sep = ""),
              open = "w", 
              nbchar = 10000, 
              as.string = TRUE)                                           # Creates a fasta file for each gene, will continue on to distance matrices
}
beep(8)
rm(alignFA, gene_file, path, row)

### Best model for genes ##################################################################################################################################
# This section extracts the best (available) model for distance matrices for the genes. GTR and HKY are not options for distance matrices, so they are
# removed from the order of best models. Accounting for invariant sites is also not an available option, so that parameter is ignored. 
# The second half of this section creates a txt file for each model with all the file pathways to the genes requiring that model.

fastaFilesModel <- as.data.frame(list.files(path = "6Model/",
                                            pattern = ".csv"))            # Make a dataframe listing the csv files in the folder
colnames(fastaFilesModel) <- "File_name"                                  # Changes the column name
fastaFilesModel$Path_name <- paste("6Model/", 
                                   fastaFilesModel$File_name, 
                                   sep = "")                              # Creates a file pathway for each gene

model_code <- function(model1) {
  best_model <- case_when(
    model1 %in% c("JC", "JC+I") ~ "JC",
    model1 %in% c("JC+G", "JC+G+I") ~ "JC_G",
    model1 %in% c("K2", "K2+I") ~ "K2",
    model1 %in% c("K2+G", "K2+G+I") ~ "K2_G",
    model1 %in% c("T92", "T92+I") ~ "T92",
    model1 %in% c("T92+G", "T92+G+I") ~ "T92_G",
    model1 %in% c("TN93", "TN93+I") ~ "TN93",
    model1 %in% c("TN93+G", "TN93+G+I") ~ "TN93_G",
  )
  return(best_model)
}

best_model <- as.data.frame(matrix(ncol = 4, nrow = 0))                   # Dataframe for each gene's best model
for(row in 1:nrow(fastaFilesModel)) {                                     # Organize by model since invariant sites are not an option for distance matrices
  path <- fastaFilesModel$Path_name[row]                                  # Path to model for each gene
  gene_model_test <- read.csv(file = path)                                # Read in model test csv
  
  gene_model_test$Dist_Matr <- case_when(                                 # Sets Dist_Matr column has FALSE if that model does not exist for distance 
    gene_model_test$Model %in% c("GTR+G+I", "GTR+G", "GTR+I",             #   matrices and TRUE if they do
                                 "GTR", "HKY+G+I", "HKY+G", 
                                 "HKY+I", "HKY") ~ FALSE,
    TRUE ~ TRUE
  )
  
  gene_model_test <- subset(gene_model_test, Dist_Matr == "TRUE")         # Subsets for available models
  
  model <- as.data.frame(gene_model_test$Model[1])                        # Assign the name of the best model to model
  colnames(model) <- "Model1"                                             # Rename column
  model$File_name <- fastaFilesModel$File_name[row]                       # File name for model testing csv file
  
  model$ModelCode <- model_code(model$Model1)                             # Get the model codes
  
  best_model <- rbind(best_model, model)                                  # Combine all best models to one dataset
}
rm(gene_model_test, model, path, row)                                     # Remove unneeded variables from the for loop

best_model$Gene <- gsub(pattern = "-4212.csv", replacement = "",          # Creates a column with just the gene name
                        x = best_model$File_name)
best_model$Path_Name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/5Aligned/",
                              best_model$Gene, ".fasta", sep = "")        # New pathway to aligned fasta file

write.csv(x = best_model, file = "8Results/best_model_full.csv", row.names = FALSE)
best_model <- read.csv(file = "8Results/best_model_full.csv")
best_model <- subset(x = best_model, select = c("Gene", "Model1", "ModelCode"))
write.csv(x = best_model, file = "8Results/Best_Model.csv", row.names = FALSE)

Uniq_mods <- as.data.frame(unique(best_model$ModelCode))                  # All unique models for genes
colnames(Uniq_mods) <- "Model_Name"

for(row in 1:nrow(Uniq_mods)) {                                           # Writes a txt file with all pathways for genes of each model (for MEGAX)
  Name <- as.character(Uniq_mods$Model_Name[row])                         # Takes each model name in turn
  
  datframe <- subset(best_model, ModelCode == Name)                       # Subsets best_model according to model name
  datframe <- as.data.frame(datframe$Path_Name)                           # Keep only the pathway to fasta files
  
  write.table(datframe, file = paste(Name, ".txt", sep = ""), sep = "\n", # Creates a txt file listing the gene pathways
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
rm(datframe, Name, row)

### Closest relative ######################################################################################################################################
# This section uses the distance matrices to extract Mixta calida's and Mixta gaviniae's relatives in order. To do this, the code must first read in the
# .meg files and transform them into a more readable format. This is what the first portion of the for loop does. The second portion names the relatives
# in order (more recent the evolutionary divide, the lower the distance number).

megFiles <- as.data.frame(list.files(path = "7Distance/",
                                     pattern = ".meg"))                   # Makes a list of all .meg file in this diretory
colnames(megFiles) <- "File_name"                                         # Changes the column name
megFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/7Distance/",
                            megFiles$File_name, sep = "")                 # Adds the path name for each gene
megFiles$Gene <- best_model$Gene                                          # Adds the gene name (no extension)

close_relative <- function(gen, spcs) {                                   # Returns the row number of the three closest relatives in the matrices
  min1 <- Rfast::nth(x = spcs, k = 1, descending = FALSE, 
                     index.return = TRUE)
  min2 <- Rfast::nth(x = spcs, k = 2, descending = FALSE, 
                     index.return = TRUE)
  min3 <- Rfast::nth(x = spcs, k = 3, descending = FALSE, 
                     index.return = TRUE)
  min4 <- Rfast::nth(x = spcs, k = 4, descending = FALSE, 
                     index.return = TRUE)
  min5 <- Rfast::nth(x = spcs, k = 5, descending = FALSE, 
                     index.return = TRUE)
  min6 <- Rfast::nth(x = spcs, k = 6, descending = FALSE, 
                     index.return = TRUE)
  min7 <- Rfast::nth(x = spcs, k = 7, descending = FALSE, 
                     index.return = TRUE)
  min8 <- Rfast::nth(x = spcs, k = 8, descending = FALSE, 
                     index.return = TRUE)
  min9 <- Rfast::nth(x = spcs, k = 9, descending = FALSE, 
                     index.return = TRUE)
  min10 <- Rfast::nth(x = spcs, k = 10, descending = FALSE, 
                      index.return = TRUE)
  
  rela <- data.frame(Gene = gen,
                     One = relative(min1),
                     Two = relative(min2),
                     Three = relative(min3),
                     Four = relative(min4),
                     Five = relative(min5),
                     Six = relative(min6),
                     Seven = relative(min7),
                     Eight = relative(min8),
                     Nine = relative(min9),
                     Ten = relative(min10))
}

relative <- function(number) {                                            # Returns the relative name given row number
  rltv <- case_when(
    number == 1 ~ "Tatumella_saanichensis",
    number == 2 ~ "Citrobacter_freundii",
    number == 3 ~ "Enterobacter_cloacae",
    number == 4 ~ "Erwinia_amylovora",
    number == 5 ~ "Erwinia_tasmaniensis",
    number == 6 ~ "Mixta_calida",
    number == 7 ~ "Mixta_gaviniae",
    number == 8 ~ "Pantoea_agglomerans",
    number == 9 ~ "Pantoea_septica",
    number == 10 ~ "Tatumella_ptyseos")
}

M_cal_rel <- as.data.frame(matrix(ncol = 11, nrow = 0))                   # Dataframe for M. calida's closest relatives
colnames(M_cal_rel) <- c("Gene", "Itself_check", "First_rel", "Second_rel", "Third_rel", "Fourth_rel", "Fifth_rel", "Sixth_rel", "Seventh_rel",
                         "Eighth_rel", "Ninth_rel")                        # Changes the column names

M_gav_rel <- as.data.frame(matrix(ncol = 11, nrow = 0))                   # Dataframe for M. gaviniae's closest relatives
colnames(M_gav_rel) <- c("Gene", "Itself_check", "First_rel", "Second_rel", "Third_rel", "Fourth_rel", "Fifth_rel", "Sixth_rel", "Seventh_rel",
                         "Eighth_rel", "Ninth_rel")                      # Changes the column names

for(row in 1:nrow(megFiles)) {                                            # Finds the two closest relatives to Mixta species
  path <- megFiles$Path_name[row]                                         # Path name
  gene <- megFiles$Gene[row]                                              # Gene name
  
  mega <- case_when(
    gene %in% c("37818_hypothetical_protein", "38262_ygbE", "38956_hypothetical_protein", "39709_yciH", "39916_eamA") 
    ~ read.table(file = path, stringsAsFactors = FALSE, skip = 37,        # These five genes had to be run manually (therefore different format)
                 fill = TRUE),
    TRUE ~ read.table(file = path, stringsAsFactors = FALSE, skip = 45,   # For the rest
                      fill = TRUE)
  )
  
  mega2 <- as.data.frame(matrix(ncol = 1, nrow = 10))
  for(i in 1:length(mega)) {                                              # Removes the square brackets
    hel <- as.character(mega[[i]])
    
    for(j in 1:length(hel)) {
      hel[j] <- gsub(pattern = "\\[|\\]", replacement = "", x = hel[j])
    }
    mega2 <- cbind(mega2, hel, stringsAsFactors = FALSE) 
  }
  rm(hel, i, j)
  colnames(mega2) <- paste("V", 1:13, sep = "")                           # Changes the column names
  
  dist <- subset(mega2, select = V3:V12)                                  # Subsets mega2, keeping only the important columns
  colnames(dist) <- paste("V", 1:10, sep = "")                            # Changes the column names
  
  for(row in 2:(nrow(dist) - 1)) {                                        # Moves things over so that the SE's are separate
    for(i in 1:(row - 1)) {
      dist[row, i] <- dist[row, i + 1]
    }
  }
  rm(i, row)
  diag(dist) <- 0                                                         # Adds zeros down the diagonal since each species' gene is closest to itself
  
  M_cal <- as.numeric(rbind(t(dist[6, 1:6]), dist[7, 6], dist[8, 6],      # Grabs the distances for each species relative to M. calida
                            dist[9, 6], dist[10, 6]))
  M_gav <- as.numeric(rbind(t(dist[7, 1:7]), dist[8, 7], dist[9, 7],      # Grabs the distances for each species relative to M. gaviniae
                            dist[10, 7]))
  
  MCclorel <- close_relative(gene, M_cal)                                 # Calls the close_relative function
  MGclorel <- close_relative(gene, M_gav)
  
  M_cal_rel <- rbind(M_cal_rel, MCclorel)                                 # Combines everything together
  M_gav_rel <- rbind(M_gav_rel, MGclorel)
}
beep(8)
rm(dist, MCclorel, mega, mega2, MGclorel, gene, M_cal, M_gav, path)

write.csv(x = M_cal_rel, file = "8Results/M_calida_Relatives.csv", row.names = FALSE)
write.csv(x = M_gav_rel, file = "8Results/M_gaviniae_Relatives.csv", row.names = FALSE)

### Categorizing closest relative #########################################################################################################################
# This section simply identifies the closest Mixta relative to both Mixta species (besides itself) and the closest non-Mixta relative. Identifying the 
# closest Mixta relative allows me to identify any genes wherein the M. calida and M. gaviniae copies are perfectly identical (distance = 0.0000).
# Identifying the closest non-Mixta relative allows me to see with which species Mixta has the most recent ancestor.

M_cal_rel <- read.csv(file = "8Results/M_calida_Relatives.csv",           # Reads in the M. calida results
                      stringsAsFactors = FALSE)

M_cal_rel$Results_Mixta <- case_when(                                     # Grabs the closest Mixta relative that is not itself, used as a check
  M_cal_rel$One == "Mixta_calida" ~ M_cal_rel$Two,
  M_cal_rel$One == "Mixta_gaviniae" ~ M_cal_rel$One
)

M_cal_rel$Results_Other <- case_when(
  M_cal_rel$Two %in% c("Mixta_calida", "Mixta_gaviniae") ~ M_cal_rel$Three,
  TRUE ~ M_cal_rel$Two
)

write.csv(x = M_cal_rel, file = "8Results/M_calida_Relatives.csv", 
          row.names = FALSE)                                              # Write the results to a csv file

####
M_gav_rel <- read.csv(file = "8Results/M_gaviniae_Relatives.csv",         # Reads in the M. gaviniae results   
                      stringsAsFactors = FALSE)

for(row in 1:nrow(M_gav_rel)) {                                           # To check if the first closest relative is the other Mixta species
  M_gav_rel$Results_Mixta[row] <- case_when(
    M_gav_rel$One[row] == "Mixta_gaviniae" ~ M_gav_rel$Two[row],
    M_gav_rel$One[row] == "Mixta_calida" ~ M_gav_rel$One[row]
  )
}

for(row in 1:nrow(M_gav_rel)) {                                          # Check the next (non-Mixta) relative
  M_gav_rel$Results_Other[row] <- case_when(
    M_gav_rel$Two[row] %in% c("Mixta_calida", "Mixta_gaviniae") ~ 
      M_gav_rel$Three[row],
    TRUE ~ M_gav_rel$Two[row]
  )
}

write.csv(x = M_gav_rel, file = "8Results/M_gaviniae_Relatives.csv", 
          row.names = FALSE)                                              # Write the results to a csv file



length(which(M_cal_rel$Results_Other == "Pantoea_septica"))               # 828
length(which(M_cal_rel$Results_Other == "Pantoea_agglomerans"))           # 81
length(which(M_cal_rel$Results_Other == "Erwinia_amylovora"))             # 49
length(which(M_cal_rel$Results_Other == "Erwinia_tasmaniensis"))          # 53
length(which(M_cal_rel$Results_Other == "Tatumella_ptyseos"))             # 4
length(which(M_cal_rel$Results_Other == "Tatumella_saanichensis"))        # 2
length(which(M_cal_rel$Results_Other == "Citrobacter_freundii"))          # 8
length(which(M_cal_rel$Results_Other == "Enterobacter_cloacae"))          # 10

length(which(M_gav_rel$Results_Other == "Pantoea_septica"))               # 838
length(which(M_gav_rel$Results_Other == "Pantoea_agglomerans"))           # 68
length(which(M_gav_rel$Results_Other == "Erwinia_amylovora"))             # 50
length(which(M_gav_rel$Results_Other == "Erwinia_tasmaniensis"))          # 55
length(which(M_gav_rel$Results_Other == "Tatumella_ptyseos"))             # 3
length(which(M_gav_rel$Results_Other == "Tatumella_saanichensis"))        # 3
length(which(M_gav_rel$Results_Other == "Citrobacter_freundii"))          # 7
length(which(M_gav_rel$Results_Other == "Enterobacter_cloacae"))          # 11

### Distances and standard errors #########################################################################################################################
# This section grabs the distances and standard errors for each species in comparison to the Mixta species (NOT the tidy version).

megFiles <- as.data.frame(list.files(path = "7Distance/",
                                     pattern = ".meg"))                   # Makes a list of all .meg file in this diretory
colnames(megFiles) <- "File_name"                                         # Changes the column name
megFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/7Distance/",
                            megFiles$File_name, sep = "")                 # Adds the path name for each gene
megFiles$Gene <- best_model$Gene                                          # Adds the gene name (no extension)


M_calida_dist <- as.data.frame(matrix(ncol = 21, nrow = 0))               # Dataframe for M. calida's closest relatives

M_gaviniae_dist <- as.data.frame(matrix(ncol = 21, nrow = 0))             # Dataframe for M. gaviniae's closest relatives

for(row in 1:nrow(megFiles)) {                                            # Finds the two closest relatives to Mixta species
  path <- megFiles$Path_name[row]                                         # Path name
  gene <- megFiles$Gene[row]                                              # Gene name
  
  if(gene %in% c("37818_hypothetical_protein", "38262_ygbE", "38956_hypothetical_protein", "39709_yciH", 
                 "39916_eamA")) {                                         # These five genes had to be run manually (therefore different format)
    mega <- read.table(file = path, stringsAsFactors = FALSE, skip = 37, 
                       fill = TRUE)
  } else {                                                                # For the rest
    mega <- read.table(file = path, stringsAsFactors = FALSE, skip = 45, 
                       fill = TRUE)
  }
  
  mega2 <- as.data.frame(matrix(ncol = 1, nrow = 10))
  for(i in 1:length(mega)) {                                              # Removes the square brackets
    hel <- as.character(mega[[i]])
    
    for(j in 1:length(hel)) {
      hel[j] <- gsub(pattern = "\\[|\\]", replacement = "", x = hel[j])
    }
    mega2 <- cbind(mega2, hel, stringsAsFactors = FALSE) 
  }
  rm(hel, i, j)
  colnames(mega2) <- paste("V", 1:13, sep = "")                           # Changes the column names
  
  dist <- subset(mega2, select = V3:V12)                                  # Subsets mega2, keeping only the important columns
  colnames(dist) <- paste("V", 1:10, sep = "")                            # Changes the column names
  
  for(row in 2:(nrow(dist) - 1)) {                                        # Moves things over so that the SE's are separate
    for(i in 1:(row - 1)) {
      dist[row, i] <- dist[row, i + 1]
    }
  }
  rm(i, row)
  diag(dist) <- 0                                                         # Adds zeros down the diagonal since each species' gene is closest to itself
  
  M_cal <- as.data.frame(t(as.numeric(cbind(dist[6, 1], dist[1, 6],       # Grabs the distances and standard errors for each species relative to M. calida
                                            dist[6, 2], dist[2, 6], 
                                            dist[6, 3], dist[3, 6],
                                            dist[6, 4], dist[4, 6],
                                            dist[6, 5], dist[5, 6],
                                            dist[6, 6], dist[6, 6],
                                            dist[7, 6], dist[6, 7],
                                            dist[8, 6], dist[6, 8],
                                            dist[9, 6], dist[6, 9],
                                            dist[10, 6], dist[6, 10]))))
  M_cal <- cbind(gene, M_cal)
  
  M_gav <- as.data.frame(t(as.numeric(cbind(dist[7, 1], dist[1, 7],       # Grabs the distances and standard errorsfor each species relative to M. gaviniae
                                            dist[7, 2], dist[2, 7], 
                                            dist[7, 3], dist[3, 7],
                                            dist[7, 4], dist[4, 7],
                                            dist[7, 5], dist[5, 7],
                                            dist[7, 6], dist[6, 7],
                                            dist[7, 7], dist[7, 7],
                                            dist[8, 7], dist[7, 8],
                                            dist[9, 7], dist[7, 9],
                                            dist[10, 7], dist[7, 10]))))
  M_gav <- cbind(gene, M_gav)
  
  M_calida_dist <- rbind(M_calida_dist, M_cal)
  M_gaviniae_dist <- rbind(M_gaviniae_dist, M_gav)
}
beep(8)
rm(dist, M_cal, M_gav, mega, mega2, gene, path)

colnames(M_calida_dist) <- c("Gene", "Tatumella_saanichensis", "TS_Error", "Citrobacter_freundii", "CF_Error", "Enterobacter_cloacae", "EC_Error",
                             "Erwinia_amylovora", "EA_Error", "Erwinia_tasmaniensis", "ET_Error", "Mixta_calida", "MC_Error", "Mixta_gaviniae", 
                             "MG_Error", "Pantoea_agglomerans", "PA_Error", "Pantoea_septica", "PS_Error", 
                             "Tatumella_ptyseos", "TP_Error")             # Changes the column names

colnames(M_gaviniae_dist) <- c("Gene", "Tatumella_saanichensis", "TS_Error", "Citrobacter_freundii", "CF_Error", "Enterobacter_cloacae", "EC_Error",
                               "Erwinia_amylovora", "EA_Error", "Erwinia_tasmaniensis", "ET_Error", "Mixta_calida", "MC_Error", "Mixta_gaviniae", 
                               "MG_Error", "Pantoea_agglomerans", "PA_Error", "Pantoea_septica", "PS_Error", 
                               "Tatumella_ptyseos", "TP_Error")           # Changes the column names

write.csv(x = M_calida_dist, file = "8Results/M_calida_Distances.csv", row.names = FALSE)
write.csv(x = M_gaviniae_dist, file = "8Results/M_gaviniae_Distances.csv", row.names = FALSE)

### Retrieve gene IDs #####################################################################################################################################
# This section retrieves the gene IDs for the Mixta species.

M_calida_dist <- read.csv(file = "8Results/M_calida_Distances.csv", 
                          stringsAsFactors = FALSE)
M_gaviniae_dist <- read.csv(file = "8Results/M_gaviniae_Distances.csv", 
                            stringsAsFactors = FALSE)

fastaFiles <- as.data.frame(list.files(path = "3Homologues_10/",
                                       pattern = ".fasta"))               # Makes a dataframe where the first column is a list of fasta gene files
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/3Homologues_10/", 
                              fastaFiles$File_name, sep = "")             # Creates a file pathway for each gene

ID_calida <- as.data.frame(matrix(ncol = 2, nrow = 0))
ID_gaviniae <- as.data.frame(matrix(ncol = 2, nrow = 0))
for(row in 1:nrow(fastaFiles)) {                                          # Lets pass genes that meet requirements
  path <- fastaFiles$Path_name[row]
  file <- as.character(fastaFiles$File_name[row])
  gene_file <- read.table(file = path, header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE)                       # Reads in the gene file according to the pathway, separation is newline
  
  calida <- gene_file[11, ]                                               # Get M. calida info row
  calida <- gsub(pattern = "\\|.*", replacement = "", x = calida)         # Remove everything after the "|" symbol
  calida <- gsub(pattern = ">ID:MLHGDOMH_", replacement = "", x = calida) # Remove the ID code
  calida <- cbind(file, calida)                                           # Combine file name with ID
  ID_calida <- rbind(ID_calida, calida)                                   # Combine everything
  
  gaviniae <- gene_file[13, ]                                             # Get M. gaviniae info row
  gaviniae <- gsub(pattern = "\\|.*", replacement = "", x = gaviniae)     # Remove everything after the "|" symbol
  gaviniae <- gsub(pattern = ">ID:MHMNNPCM_", replacement = "", 
                   x = gaviniae)                                          # Remove the ID code
  gaviniae <- cbind(file, gaviniae)                                       # Combine file name with ID
  ID_gaviniae <- rbind(ID_gaviniae, gaviniae)                             # Combine everything
}
rm(calida, gaviniae, gene_file, file, path, row)

colnames(ID_calida) <- c("Gene", "ID")                                    # Change the column names
colnames(ID_gaviniae) <- c("Gene", "ID")
ID_calida$Gene <- as.character(ID_calida$Gene)                            # Change structure of columns
ID_calida$ID <- as.character(ID_calida$ID)
ID_gaviniae$Gene <- as.character(ID_gaviniae$Gene)
ID_gaviniae$ID <- as.character(ID_gaviniae$ID)

ID_calida$Gene <- gsub(pattern = ".fasta", replacement = "", 
                       x = ID_calida$Gene)                                # Remove ".fasta" from file names so only gene names
ID_gaviniae$Gene <- gsub(pattern = ".fasta", replacement = "", 
                         x = ID_gaviniae$Gene)

ID_calida$Check <- case_when(
  ID_calida$Gene %in% M_calida_dist$Gene ~ TRUE,
  TRUE ~ FALSE
)

ID_gaviniae$Check <- case_when(
  ID_gaviniae$Gene %in% M_gaviniae_dist$Gene ~ TRUE,
  TRUE ~ FALSE
)

ID_calida <- subset(x = ID_calida, Check == TRUE)                         # Subset for only the 1035 genes
ID_gaviniae <- subset(x = ID_gaviniae, Check == TRUE)

M_calida_dist <- cbind(ID_calida$Gene, ID_calida$ID, M_calida_dist)       # Combine IDs with distance and standard errors dataframe
M_gaviniae_dist <- cbind(ID_gaviniae$Gene, ID_gaviniae$ID, 
                         M_gaviniae_dist)

colnames(M_calida_dist)[1:2] <- c("Gene_Check", "ID")                     # Change column names for first two columns
colnames(M_gaviniae_dist)[1:2] <- c("Gene_Check", "ID")

M_calida_dist <- subset(M_calida_dist, select = -Gene_Check)
M_gaviniae_dist <- subset(M_gaviniae_dist, select = -Gene_Check)

write.csv(x = M_calida_dist, file = "8Results/M_calida_Distances.csv",    # Save distances, standard errors, and gene IDs
          row.names = FALSE)
write.csv(x = M_gaviniae_dist, file = "8Results/M_gaviniae_Distances.csv", 
          row.names = FALSE)

### Sort distances ########################################################################################################################################
# This section tidies up the distances and standard errors in a more readable format for when making plots

M_calida_dist <- read.csv(file = "8Results/M_calida_Distances.csv", 
                          stringsAsFactors = FALSE)                       # Read in M. calida's distances
M_gaviniae_dist <- read.csv(file = "8Results/M_gaviniae_Distances.csv", 
                            stringsAsFactors = FALSE)                     # Read in M. gaviniae's distances

ID_Species_names <- c("ID", "Gene", "Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", 
                      "Erwinia_amylovora", "Erwinia_tasmaniensis", "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", 
                      "Pantoea_septica", "Tatumella_ptyseos")             # Column names for subsetting for the distances
ID_Species_errors <- c("ID", "Gene", "TS_Error", "CF_Error", "EC_Error", "EA_Error", "ET_Error", "MC_Error", "MG_Error", 
                       "PA_Error", "PS_Error", "TP_Error")                # Column names for subsetting for the standard errors

tidy_dist_C <- subset(M_calida_dist, select = ID_Species_names)           # Subset for only the distances
tidy_dist_C <- tidy_dist_C %>%                                            # Gather distances into one column according to species and by ID
  pivot_longer(cols = Tatumella_saanichensis:Tatumella_ptyseos,
               names_to = "Species", values_to = "Distance")

tidy_ster_C <- subset(M_calida_dist, select = ID_Species_errors)          # Subset for only the standard errors
tidy_ster_C <- tidy_ster_C %>%                                            # Gather standard errors into one column according to species and by ID
  pivot_longer(cols = TS_Error:TP_Error,
               names_to = "Species_se", values_to = "Std_Errors")
colnames(tidy_ster_C) <- c("ID_se", "Gene_se", "Species_se",
                           "Std_Errors")                                  # Change the column names

tidy_calida <- cbind(tidy_dist_C, tidy_ster_C)                            # Combine distances and standard errors dataframes
tidy_calida <- subset(tidy_calida, select = c("ID", "Gene", "Species", "Distance", 
                                              "Std_Errors"))              # Removes the columns that were needed only to double check things were inline

M_calida_sort <- ddply(tidy_calida, c("ID", "Gene"))                      # Sort by gene ID and then Gene name

write.csv(x = M_calida_sort, file = "8Results/M_calida_Sort_Dist.csv",    # Write this to a csv file
          row.names = FALSE)

# Mixta gaviniae #
tidy_dist_G <- subset(M_gaviniae_dist, select = ID_Species_names)         # Subset for only the distances
tidy_dist_G <- tidy_dist_G %>%                                            # Gather distances into one column according to species and by ID
  pivot_longer(cols = Tatumella_saanichensis:Tatumella_ptyseos,
               names_to = "Species", values_to = "Distance")

tidy_ster_G <- subset(M_gaviniae_dist, select = ID_Species_errors)        # Subset for only the standard errors
tidy_ster_G <- tidy_ster_G %>%                                            # Gather standard errors into one column according to species and by ID
  pivot_longer(cols = TS_Error:TP_Error,
               names_to = "Species_se", values_to = "Std_Errors")
colnames(tidy_ster_G) <- c("ID_se", "Gene_se", "Species_se", 
                           "Std_Errors")                                  # Change the column names

tidy_gaviniae <- cbind(tidy_dist_G, tidy_ster_G)                          # Combine distances and standard errors dataframes
tidy_gaviniae <- subset(tidy_gaviniae, select = c("ID", "Gene", "Species", "Distance", 
                                                  "Std_Errors"))          # Removes the columns that were needed only to double check things were inline

M_gaviniae_sort <- ddply(tidy_gaviniae, c("ID", "Gene"))                  # Sort by gene ID and then Gene name

write.csv(x = M_gaviniae_sort, file = "8Results/M_gaviniae_Sort_Dist.csv", 
          row.names = FALSE)                                              # Write this to a csv file


### M. calida plots #######################################################################################################################################
# This section creates some plots using Mixta calida data

M_calida_sort <- read.csv(file = "8Results/M_calida_Sort_Dist.csv", 
                          stringsAsFactors = FALSE)                       # Read in the tidied and sorted distances for M. calida

M_cal1 <- subset(M_calida_sort, Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                               "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                               "Tatumella_ptyseos"))      # Removes the Mixta species since they are most likely ~ 0

M_cal1$DistanceN <- M_cal1$Distance * -1                                  # Creates a column with negative distances (so 0 will be at top of plot)

extra_genes <- as.data.frame(matrix(ncol = 0, nrow = 42))                 # Ensures data set has 1:number of genes M. calida actually has
extra_genes <- mutate(extra_genes,
                      ID = 4043:4084, Gene = NA, 
                      Species = "Tatumella_saanichensis", Distance = NA,
                      Std_Errors = NA, DistanceN = NA)

M_cal1 <- rbind(M_cal1, extra_genes)

png("9_1Plots_calida/MC_Full_dist.png", width = 2000, height = 1200)
ggplot(data = M_cal1, aes(x = ID, y = DistanceN)) +                       # Full plot
  geom_point(aes(colour = M_cal1$Species), 
             size = 2, alpha = 0.75) +
  geom_line(aes(colour = M_cal1$Species), 
            linetype = "dotted") +
  scale_colour_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                      labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                                 "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
  theme(legend.position = "bottom", text = element_text(size = 36), 
        legend.text = element_text(face = "italic")) +
  labs(x = expression(paste(italic("M. calida"), " Gene ID")), 
       y = expression(paste("Negative Distance from ", italic("M. calida"))), 
       colour = "Species") +
  scale_x_continuous(limits = c(-50, 4134), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4.0, 0.1), expand = c(0, 0))
dev.off()

# Separated into 8 groups
dist_plot <- function(dtst, beg, end) {                                   # Function for the 8 distance plots (8 segments of the full plot)
  ggplot(data = dtst, aes(x = ID, y = DistanceN)) +
    geom_point(aes(colour = dtst$Species), 
               size = 1.5, alpha = 0.75) +
    geom_line(aes(colour = dtst$Species), 
              linetype = "dotted") +
    scale_colour_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                        labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                                   "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
    theme(legend.position = "bottom", text = element_text(size = 9), 
          legend.text = element_text(face = "italic")) +
    labs(x = expression(paste(italic("M. calida"), " Gene ID")), 
         y = expression(paste("Negative Distance from ", italic("M. calida"))), 
         colour = "Species") +
    scale_x_continuous(breaks = round(seq(min(dtst$ID), max(dtst$ID), by = 100), -2),
                       limits = c(beg - 10, end + 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-5.11, 0.1), expand = c(0, 0)) +
    geom_errorbar(aes(ymin = DistanceN - Std_Errors, 
                      ymax = DistanceN + Std_Errors, colour = Species), 
                  width = 0.2, position = position_dodge(0.05))
}

gene_num <- as.data.frame(matrix(ncol = 0, nrow = 8))                     # To separate the dataset into 8 segments of equal length on the genome
gene_num <- mutate(gene_num,
                   beg = c(1, 511, 1022, 1532, 2043, 2553, 3064, 3574),
                   end = c(510, 1021, 1531, 2042, 2552, 3063, 3573, 4084))

for(row in 1:nrow(gene_num)) {                                            # Creates the plots and saves them
  plot <- dist_plot(subset(M_cal1, ID %in% gene_num$beg[row]:gene_num$end[row]), gene_num$beg[row], gene_num$end[row])
  ggsave(plot, file = paste("9_1Plots_calida/MC_dist_", row, "_8.png", sep = ""), 
         width = 16.51, height = 12.38, units = "cm")
}

# Separated into 8 groups and distances no greater than 1
M_cal2 <- subset(M_cal1, DistanceN > -1)
M_cal2 <- rbind(M_cal2, extra_genes)

dist_plot_one <- function(dtst, beg, end) {                               # Function for the next 8 plots that show only distances < 1
  ggplot(data = dtst, aes(x = ID, y = DistanceN)) +
    geom_point(aes(colour = dtst$Species), 
               size = 1.5, alpha = 0.75) +
    scale_colour_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                        labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                                   "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
    theme(legend.position = "bottom", text = element_text(size = 9), 
          legend.text = element_text(face = "italic")) +
    labs(x = expression(paste(italic("M. calida"), " Gene ID")),
         y = expression(paste("Negative Distance from ", italic("M. calida"))), 
         colour = "Species") +
    scale_x_continuous(breaks = round(seq(min(dtst$ID), max(dtst$ID), by = 100), -2),
                       limits = c(beg - 10, end + 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.25, 0.05), expand = c(0, 0)) +
    geom_errorbar(aes(ymin = DistanceN - Std_Errors, 
                      ymax = DistanceN + Std_Errors, colour = Species), 
                  width = 0.2, position = position_dodge(0.05))
}

for(row in 1:nrow(gene_num)) {                                            # Creates and saves the eight plots
  plot <- dist_plot_one(subset(M_cal2, ID %in% gene_num$beg[row]:gene_num$end[row]), gene_num$beg[row], gene_num$end[row])
  ggsave(plot, file = paste("9_1Plots_calida/MC_distone_", row, "_8.png", sep = ""), 
         width = 16.51, height = 12.38, units = "cm")
}

### M. gaviniae plots #####################################################################################################################################
# This section creates some plots using Mixta gaviniae data

M_gaviniae_sort <- read.csv(file = "8Results/M_gaviniae_Sort_Dist.csv", 
                            stringsAsFactors = FALSE)                     # Read in the tidied and sorted distances for M. gaviniae

M_gav1 <- subset(M_gaviniae_sort, Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                                 "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                                 "Tatumella_ptyseos"))    # Removes the Mixta species since they are most likely ~ 0
M_gav1$DistanceN <- M_gav1$Distance * -1                                  # Creates a column with negative distances (so 0 will be at top of plot)

extra_genes <- as.data.frame(matrix(ncol = 0, nrow = 75))                 # Ensures data set has 1:number of genes M. calida actually has
extra_genes <- mutate(extra_genes,
                      ID = c(1, 4182:4255), Gene = NA, 
                      Species = "Tatumella_saanichensis", Distance = NA,
                      Std_Errors = NA, DistanceN = NA)

M_gav1 <- rbind(M_gav1, extra_genes)

png("9_2Plots_gaviniae/MG_Full_dist.png", width = 2000, height = 1200)
ggplot(data = M_gav1, aes(x = ID, y = DistanceN)) +                       # Full plot
  geom_point(aes(colour = M_gav1$Species), 
             size = 2, alpha = 0.75) +
  geom_line(aes(colour = M_gav1$Species), 
            linetype = "dotted") +
  scale_colour_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                      labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                                 "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
  theme(legend.position = "bottom", text = element_text(size = 36),
        legend.text = element_text(face = "italic")) +
  labs(x = expression(paste(italic("M. gaviniae"), " Gene ID")),
       y = expression(paste("Negative Distance from ", italic("M. gaviniae"))), 
       colour = "Species") +
  scale_x_continuous(limits = c(-50, 4305), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4.0, 0.1), expand = c(0, 0))
dev.off()

# Separated into 8 groups
dist_plot <- function(dtst, beg, end) {                                   # Function for the 8 distance plots (8 segments of the full plot)
  ggplot(data = dtst, aes(x = ID, y = DistanceN)) +
    geom_point(aes(colour = dtst$Species), 
               size = 1.5, alpha = 0.75) +
    geom_line(aes(colour = dtst$Species), 
              linetype = "dotted") +
    scale_colour_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                        labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                                   "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
    theme(legend.position = "bottom", text = element_text(size = 9),
          legend.text = element_text(face = "italic")) +
    labs(x = expression(paste(italic("M. gaviniae"), " Gene ID")),
         y = expression(paste("Negative Distance from ", italic("M. gaviniae"))), 
         colour = "Species") +
    scale_x_continuous(breaks = round(seq(min(dtst$ID), max(dtst$ID), by = 100), -2),
                       limits = c(beg - 10, end + 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-5.11, 0.1), expand = c(0, 0)) +
    geom_errorbar(aes(ymin = DistanceN - Std_Errors, 
                      ymax = DistanceN + Std_Errors, colour = Species), 
                  width = 0.2, position = position_dodge(0.05))
}

gene_num <- as.data.frame(matrix(ncol = 0, nrow = 8))                     # To separate the dataset into 8 segments of equal length on the genome
gene_num <- mutate(gene_num,
                   beg = c(1, 532, 1064, 1596, 2128, 2660, 3192, 3724),
                   end = c(531, 1063, 1595, 2127, 2659, 3191, 3723, 4255))

for(row in 1:nrow(gene_num)) {                                            # Creates the plots and saves them
  plot <- dist_plot(subset(M_gav1, ID %in% gene_num$beg[row]:gene_num$end[row]), gene_num$beg[row], gene_num$end[row])
  ggsave(plot, file = paste("9_2Plots_gaviniae/MG_dist_", row, "_8.png", sep = ""), 
         width = 16.51, height = 12.38, units = "cm")
}

# Separated into 8 groups and distances no greater than 1
dist_plot_one <- function(dtst, beg, end) {                               # Function for the next 8 plots that show only distances < 1
  ggplot(data = dtst, aes(x = ID, y = DistanceN)) +
    geom_point(aes(colour = dtst$Species), size = 1.5, alpha = 0.75) +
    scale_colour_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                        labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                                   "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
    theme(legend.position = "bottom", text = element_text(size = 9),
          legend.text = element_text(face = "italic")) +
    labs(x = expression(paste(italic("M. gaviniae"), " Gene ID")),
         y = expression(paste("Negative Distance from ", italic("M. gaviniae"))), 
         colour = "Species") +
    scale_x_continuous(breaks = round(seq(min(dtst$ID), max(dtst$ID), by = 100), -2),
                       limits = c(beg - 10, end + 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.25, 0.05), expand = c(0, 0)) +
    geom_errorbar(aes(ymin = DistanceN - Std_Errors, 
                      ymax = DistanceN + Std_Errors, colour = Species), 
                  width = 0.2, position = position_dodge(0.05))
}

M_gav2 <- subset(M_gav1, DistanceN > -1)
M_gav2 <- rbind(M_gav2, extra_genes)

for(row in 1:nrow(gene_num)) {                                            # Creates the plots and saves them
  plot <- dist_plot_one(subset(M_gav2, ID %in% gene_num$beg[row]:gene_num$end[row]), gene_num$beg[row], gene_num$end[row])
  ggsave(plot, file = paste("9_2Plots_gaviniae/MG_distone_", row, "_8.png", sep = ""), 
         width = 16.51, height = 12.38, units = "cm")
}

### Circular plots ########################################################################################################################################
# This section creates a circular plot showing the clostest non-Mixta relatives around the Mixta genome.

# Mixta calida #
M_calida <- read.csv(file = "8Results/M_calida_Relatives.csv",            # Reads in the M. calida relative results
                     stringsAsFactors = FALSE)

M_calida_dist <- read.csv(file = "8Results/M_calida_Distances.csv",
                          stringsAsFactors = FALSE)

M_calida$ID <- M_calida_dist$ID                                           # Copies the ID column into M_calida

M_calida$Results_Number <- case_when(                                     # Replaces the results names with a number
  M_calida$Results_Other == "Pantoea_agglomerans" ~ 4,
  M_calida$Results_Other == "Pantoea_septica" ~ 3,
  M_calida$Results_Other == "Erwinia_amylovora" ~ 6,
  M_calida$Results_Other == "Erwinia_tasmaniensis" ~ 5,
  M_calida$Results_Other == "Tatumella_ptyseos" ~ 7,
  M_calida$Results_Other == "Tatumella_saanichensis" ~ 8,
  M_calida$Results_Other == "Citrobacter_freundii" ~ 9,
  M_calida$Results_Other == "Enterobacter_cloacae" ~ 10
)

extra_genes <- as.data.frame(matrix(ncol = 0, nrow = 42))                 # Ensures data set has 1:number of genes M. calida actually has
extra_genes <- mutate(extra_genes,
                      Gene = NA, One = NA, Two = NA, Three = NA, 
                      Four = NA, Five = NA, Six = NA, Seven = NA, 
                      Eight = NA, Nine = NA, Ten = NA, 
                      Results_Mixta = NA, 
                      Results_Other = "Tatumella_saanichensis",           # Arbitrary Results_Other name
                      ID = 4043:4084, Results_Number = NA)

M_calida <- rbind(M_calida, extra_genes)                                  # Combines M_calida with the extra genes

# Mixta gaviniae #
M_gaviniae <- read.csv(file = "8Results/M_gaviniae_Relatives.csv",        # Reads in the M. gaviniae results
                       stringsAsFactors = FALSE)

M_gaviniae_dist <- read.csv(file = "8Results/M_gaviniae_Distances.csv",
                            stringsAsFactors = FALSE)

M_gaviniae$ID <- M_gaviniae_dist$ID                                       # Copies the ID column into M_gaviniae

M_gaviniae$Results_Number <- case_when(                                   # Replaces the results names with a number
  M_gaviniae$Results_Other == "Pantoea_agglomerans" ~ 4,
  M_gaviniae$Results_Other == "Pantoea_septica" ~ 3,
  M_gaviniae$Results_Other == "Erwinia_amylovora" ~ 6,
  M_gaviniae$Results_Other == "Erwinia_tasmaniensis" ~ 5,
  M_gaviniae$Results_Other == "Tatumella_ptyseos" ~ 7,
  M_gaviniae$Results_Other == "Tatumella_saanichensis" ~ 8,
  M_gaviniae$Results_Other == "Citrobacter_freundii" ~ 9,
  M_gaviniae$Results_Other == "Enterobacter_cloacae" ~ 10
)

extra_genes <- as.data.frame(matrix(ncol = 0, nrow = 75))                 # Ensures data set has 1:number of genes M. calida actually has
extra_genes <- mutate(extra_genes,
                      Gene = NA, One = NA, Two = NA, Three = NA, 
                      Four = NA, Five = NA, Six = NA, Seven = NA, 
                      Eight = NA, Nine = NA, Ten = NA, 
                      Results_Mixta = NA, 
                      Results_Other = "Tatumella_saanichensis", 
                      ID = c(1, 4182:4255), Results_Number = NA)

M_gaviniae <- rbind(M_gaviniae, extra_genes)                              # Combines M_calida with the extra genes

# M. calida plot #
png("9_1Plots_calida/MC_categ_results.png", width = 1000, height = 725)
ggplot(data = M_calida, aes(xmin = ID - 5, xmax = ID, ymin = 0, ymax = Results_Number, fill = Results_Other)) +
  geom_rect() +
  scale_fill_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                    labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                               "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
  coord_polar() +
  scale_y_continuous(limits = c(0, 10)) +
  theme(legend.position = "right", text = element_text(size = 20), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.text = element_text(face = "italic")) +
  labs(fill = "Species")
dev.off()

# M. gaviniae plot #
png("9_2Plots_gaviniae/MG_categ_results.png", width = 1000, height = 725)
ggplot(data = M_gaviniae, aes(xmin = ID - 5, xmax = ID, ymin = 0, ymax = Results_Number, fill = Results_Other)) +
  geom_rect() +
  scale_fill_manual(values = alpha(c("red", "orange", "green3", "darkgreen", "dodgerblue2", "blue3", "darkorchid", "violetred1")),
                    labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                               "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
  coord_polar() +
  scale_y_continuous(limits = c(0, 10)) +
  theme(legend.position = "right", text = element_text(size = 20), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.text = element_text(face = "italic")) +
  labs(fill = "Species")
dev.off()

#
### Gene Length ###########################################################################################################################################
# This section retrieves the gene lengths for each gene in each species

# Min, Mean, and Max #
M_calida <- read.csv(file = "8Results/M_calida_Relatives.csv",            # Reads in the M. calida results
                     stringsAsFactors = FALSE)
M_gaviniae <- read.csv(file = "8Results/M_gaviniae_Relatives.csv",
                       stringsAsFactors = FALSE)

fastaFiles <- as.data.frame(list.files(path = "3Homologues_10/",
                                       pattern = ".fasta"))               # Makes a dataframe where the first column is a list of fasta gene files
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("3Homologues_10/", 
                              fastaFiles$File_name, 
                              sep = "")                                   # Creates a file pathway for each gene
fastaFiles$Gene <- gsub(pattern = ".fasta", replacement = "",             # Creates a column with just the gene name
                        x = fastaFiles$File_name)

fastaFiles$Check <- case_when(                                            # Sets common genes to TRUE
  fastaFiles$Gene %in% M_calida$Gene ~ TRUE,
  TRUE ~ FALSE
)

fastaFiles <- subset(fastaFiles, Check == TRUE, select = File_name:Gene)  # Subsets based on TRUE so file contains only the 1035 genes

gene_length_mmm <- function(gene_file) {                                  # Finds the average, min, and max gene lengths in the file
  leng <- as.data.frame(matrix(ncol = 0, nrow = 0))
  for(row in 1:nrow(gene_file)) {
    len <- nchar(gene_file$sequences)[row]
    leng <- rbind(leng, len)
  }
  colnames(leng) <- "gene_mmm"
  mean_length <- mean(leng$gene_mmm, na.rm = TRUE)
  lengths <- as.data.frame(cbind(min(nchar(gene_file$sequences)), 
                                 mean_length, 
                                 max(nchar(gene_file$sequences))))
}

gene_mmm <- as.data.frame(matrix(ncol = 0, nrow = 0))
for(row in 1:nrow(fastaFiles)) {                                          # Average, min, and max gene lengths for each of 1035 gene files
  path <- fastaFiles$Path_name[row]
  gene_file <- read.table(file = path, 
                          header = FALSE, 
                          sep = "\n", 
                          stringsAsFactors = FALSE)                       # Reads in the gene file according to the pathway, separation is newline
  
  gene_file <- data.frame(sequences = gene_file$V1[1:10 * 2], 
                          species = gene_file$V1[1:10 * 2 - 1], 
                          stringsAsFactors = FALSE)                       # Dataframe where first column are sequences and second are corresponing names
  
  gene_file$species <- c("Tatumella saanichensis__NML_06-3099", "Citrobacter freundii__NCTC_9750", "Enterobacter cloacae_subsp_cloacae__ATCC 13047", 
                         "Erwinia amylovora__CFBP_1232", "Erwinia tasmaniensis__ET1-99", "Mixta calida__DSM_22759", "Mixta gaviniae__DSM_22758", 
                         "Pantoea agglomerans__NBRC_102470", "Pantoea septica__LMG_5345", 
                         "Tatumella ptyseos__NCTC_11468")                 # Renames the species names so they aren't ridiculously long
  
  lengths <- gene_length_mmm(gene_file)                                   # Calls the function
  lengths <- cbind(lengths, fastaFiles$Gene[row])                         # Adds the gene name, this will be used to check the gene order for cbind later
  colnames(lengths) <- c("Min", "Mean", "Max", "Gene_Name")
  
  gene_mmm <- rbind(gene_mmm, lengths)                                    # Creates a dataframe with all of the average, min, and max gene lengths
}
rm(gene_file, lengths, path, row)                                         # Removes unnecessary variables

unique(fastaFiles$Gene == gene_mmm$Gene_Name)                             # Checks that gene order is the same between both files. Should return just TRUE

M_calida <- cbind(M_calida, gene_mmm)                                     # Combines the two dataframes
M_calida <- subset(M_calida, select = Gene:Max)                           # Removes the Gene_name column

M_calida$Rela_Pattern <- paste(M_calida$One, M_calida$Two, M_calida$Three, 
                               M_calida$Four, M_calida$Five, M_calida$Six, 
                               M_calida$Seven, M_calida$Eight, 
                               M_calida$Nine, M_calida$Ten, sep = "_")

species <- as.data.frame(matrix(ncol = 0, nrow = 10))
species <- mutate(species,
                  Spp = c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", "Erwinia_tasmaniensis", 
                          "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", "Pantoea_septica", "Tatumella_ptyseos"),
                  Acro = c("TS", "CF", "EC", "EA", "ET", "MC", "MG", "PA", "PS", "TP"))

for(row in 1:nrow(species)) {
  M_calida$Rela_Pattern <- gsub(pattern = species$Spp[row], replacement = species$Acro[row], x = M_calida$Rela_Pattern)
}
rm(row)

write.csv(x = M_calida, file = "8Results/M_calida_Relatives_Length.csv", 
          row.names = FALSE)

M_gaviniae <- cbind(M_gaviniae, M_calida$Min, M_calida$Mean, 
                    M_calida$Max)

M_gaviniae$Rela_Pattern <- paste(M_gaviniae$One, M_gaviniae$Two, 
                                 M_gaviniae$Three, M_gaviniae$Four, 
                                 M_gaviniae$Five, M_gaviniae$Six, 
                                 M_gaviniae$Seven, M_gaviniae$Eight, 
                                 M_gaviniae$Nine, M_gaviniae$Ten, 
                                 sep = "_")

for(row in 1:nrow(species)) {
  M_gaviniae$Rela_Pattern <- gsub(pattern = species$Spp[row], replacement = species$Acro[row], x = M_gaviniae$Rela_Pattern)
}
rm(row)

colnames(M_gaviniae)[14:16] <- c("Min", "Mean", "Max")

write.csv(x = M_gaviniae, file = "8Results/M_gaviniae_Relatives_Length.csv", 
          row.names = FALSE)

### Individual gene lengths ###
M_calida_sdist <- read.csv(file = "8Results/M_calida_Sort_Dist.csv", 
                           stringsAsFactors = FALSE)
M_calida_sdist_L <- M_calida_sdist[order(M_calida_sdist$Gene),]           # Reorders the dataset by gene name

gene_length_all <- function(gene_file) {                                  # Finds the individual gene lengths for each species in the file
  leng <- as.data.frame(matrix(ncol = 0, nrow = 0))
  for(row in 1:nrow(gene_file)) {
    len <- nchar(gene_file$sequences[row])
    leng <- rbind(leng, len)
  }
  colnames(leng) <- "dist_all"
  return(leng)
}

gene_leng <- as.data.frame(matrix(ncol = 0, nrow = 0))
for(row in 1:nrow(fastaFiles)) {                                          # Gene lengths for each species for each of the 1035 gene files
  path <- fastaFiles$Path_name[row]
  gene_file <- read.table(file = path, 
                          header = FALSE, 
                          sep = "\n", 
                          stringsAsFactors = FALSE)                       # Reads in the gene file according to the pathway, separation is newline
  
  gene_file <- data.frame(sequences = gene_file$V1[1:10 * 2], 
                          species = gene_file$V1[1:10 * 2 - 1], 
                          stringsAsFactors = FALSE)                       # Dataframe where first column are sequences and second are corresponing names
  
  gene_file$species <- c("Tatumella saanichensis__NML_06-3099", "Citrobacter freundii__NCTC_9750", "Enterobacter cloacae_subsp_cloacae__ATCC 13047", 
                         "Erwinia amylovora__CFBP_1232", "Erwinia tasmaniensis__ET1-99", "Mixta calida__DSM_22759", "Mixta gaviniae__DSM_22758", 
                         "Pantoea agglomerans__NBRC_102470", "Pantoea septica__LMG_5345", 
                         "Tatumella ptyseos__NCTC_11468")                 # Renames the species names so they aren't ridiculously long
  
  lengths <- gene_length_all(gene_file)                                   # Calls the function
  lengths <- cbind(lengths, fastaFiles$Gene[row])                         # Adds the gene name, this will be used to check the gene order for cbind later
  colnames(lengths) <- c("Gene_Lengths", "Gene_Name")
  
  gene_leng <- rbind(gene_leng, lengths)                                  # Creates a dataframe with all of the distances
}
rm(gene_file, lengths, path, row)                                         # Removes unnecessary variables

unique(gene_leng$Gene_Name == M_calida_sdist_L$Gene)                      # The check, should return just TRUE

M_calida_sdist_L$Gene_Length <- gene_leng$Gene_Lengths                    # Copies gene lengths into M_calida_sdist_L
M_calida_sdist_L <- M_calida_sdist_L[order(M_calida_sdist_L$ID),]         # Reorders the dataset by ID

write.csv(x = M_calida_sdist_L, file = "8Results/M_calida_Sort_Dist_L.csv", 
          row.names = FALSE)

###
M_gaviniae_sdist <- read.csv(file = "8Results/M_gaviniae_Sort_Dist.csv", 
                             stringsAsFactors = FALSE)
M_gaviniae_sdist_L <- M_gaviniae_sdist[order(M_gaviniae_sdist$Gene),]     # Reorders the dataset by gene name
gene_leng
unique(gene_leng$Gene_Name == M_gaviniae_sdist_L$Gene)                    # The check, should return just TRUE

M_gaviniae_sdist_L$Gene_Length <- gene_leng$Gene_Lengths                  # Copies gene lengths into M_calida_sdist_L
M_gaviniae_sdist_L <- M_gaviniae_sdist_L[order(M_gaviniae_sdist_L$ID),]   # Reorders the dataset by ID

write.csv(x = M_gaviniae_sdist_L, file = "8Results/M_gaviniae_Sort_Dist_L.csv", 
          row.names = FALSE)
#
### Gene Length Models ####################################################################################################################################
# This section contains models that aim to identify if there is a significant relationship between (1) closest first relative and gene length and (2) 
# distances and gene length.

library("mgcv")
library("gratia")

M_calida_L <- read.csv(file = "8Results/M_calida_Relatives_Length.csv", 
                       stringsAsFactors = FALSE)
M_gaviniae_L <- read.csv(file = "8Results/M_gaviniae_Relatives_Length.csv",
                         stringsAsFactors = FALSE)
M_calida_sdist_L <- read.csv(file = "8Results/M_calida_Sort_Dist_L.csv",
                             stringsAsFactors = FALSE)
M_gaviniae_sdist_L <- read.csv(file = "8Results/M_gaviniae_Sort_Dist_L.csv",
                               stringsAsFactors = FALSE)

ctrl <- gam.control(nthreads = 3, trace = TRUE)                           # Control for GAM

M_calida_L <- mutate(M_calida_L,
                     Gene = as.factor(Gene),
                     Rela_Pattern = as.factor(Rela_Pattern),
                     Results_Other = as.factor(Results_Other),
                     ROther_Num = as.integer(Results_Other) - 1,
                     RP_Num = as.integer(Rela_Pattern) - 1)

M_gaviniae_L <- mutate(M_gaviniae_L,
                       Gene = as.factor(Gene),
                       Rela_Pattern = as.factor(Rela_Pattern),
                       Results_Other = as.factor(Results_Other),
                       ROther_Num = as.integer(Results_Other) - 1,
                       RP_Num = as.integer(Rela_Pattern) - 1)

### Closest relative vs gene length ###
CR_GLc <- readRDS("10Models/MC_Clos_Rel.rds")
CR_GLg <- readRDS("10Models/MG_Clos_Rel.rds")

CR_GLc <- gam(list(ROther_Num ~ s(Mean), ~ s(Mean), ~ s(Mean), ~ s(Mean), ~ s(Mean), ~ s(Mean), ~ s(Mean)),
              data = M_calida_L, family = multinom(K=7))                  # Model to determine if mean gene length has an affect on closest relative result

saveRDS(CR_GLc, file = "10Models/MC_Clos_Rel.rds")                        # Saves the model

layout(matrix(1:4, ncol = 2, byrow = TRUE))
gam.check(CR_GLc)
draw(CR_GLc)

CR_GLg <- gam(list(ROther_Num ~ s(Mean), ~ s(Mean), ~ s(Mean), ~ s(Mean), ~ s(Mean), ~ s(Mean), ~ s(Mean)),
              data = M_gaviniae_L, family = multinom(K=7))

saveRDS(CR_GLg, file = "10Models/MG_Clos_Rel.rds")                        # Saves the model

layout(matrix(1:4, ncol = 2, byrow = TRUE))
gam.check(CR_GLg)
draw(CR_GLg)

## Predictions ##
M_cal_pred <- data.frame(Mean = 156:4221)

M_cal_pred <- cbind(M_cal_pred, predict(CR_GLc, newdata = M_cal_pred, se = TRUE, type = "response"))
colnames(M_cal_pred)[2:17] <- c("Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", "Erwinia_tasmaniensis", "Pantoea_agglomerans", 
                                "Pantoea_septica", "Tatumella_ptyseos", "Tatumella_saanichensis", 
                                "CF_SE", "EC_SE", "EA_SE", "ET_SE", "PA_SE", "PS_SE", "TP_SE", "TS_SE")

tidy_MC <- subset(M_cal_pred, select = Mean:Tatumella_saanichensis)
tidy_MC <- tidy_MC %>%
  pivot_longer(cols = Citrobacter_freundii:Tatumella_saanichensis,
               names_to = "Species", values_to = "Fit")

tidy_MC_se <- subset(M_cal_pred, select = CF_SE:TS_SE)
tidy_MC_se$Mean <- M_cal_pred$Mean
tidy_MC_se <- tidy_MC_se %>%
  pivot_longer(cols = CF_SE:TS_SE,
               names_to = "Species", values_to = "Fit_SE")

unique(tidy_MC$Mean == tidy_MC_se$Mean)

M_cal_pred <- cbind(tidy_MC, subset(tidy_MC_se, select = Fit_SE))
rm(tidy_MC, tidy_MC_se)

M_cal_pred <- mutate(M_cal_pred,
                     upper = Fit + (Fit_SE * 1.96),
                     lower = Fit - (Fit_SE * 1.96))

MC_pred <- ggplot(M_cal_pred, aes(x = Mean, y = Fit, colour = Species)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = Mean), inherit.aes = FALSE, 
              data = M_cal_pred, alpha = 0.2) +
  geom_line() +
  scale_colour_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                      labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                                 "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
  theme(legend.position = "bottom", text = element_text(size = 9),
        legend.text = element_text(face = "italic")) +
  labs(x = "Gene Length Mean", y = "Probability", 
       colour = "Species")

ggsave(MC_pred, file = "9_1Plots_calida/MC_CR_pred.png", 
       width = 16.51, height = 12.38, units = "cm")


M_gav_pred <- data.frame(Mean = 156:4221)

M_gav_pred <- cbind(M_gav_pred, predict(CR_GLg, newdata = M_gav_pred, se = TRUE, type = "response"))
colnames(M_gav_pred)[2:17] <- c("Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", "Erwinia_tasmaniensis", "Pantoea_agglomerans", 
                                "Pantoea_septica", "Tatumella_ptyseos", "Tatumella_saanichensis", 
                                "CF_SE", "EC_SE", "EA_SE", "ET_SE", "PA_SE", "PS_SE", "TP_SE", "TS_SE")

tidy_MG <- subset(M_gav_pred, select = Mean:Tatumella_saanichensis)
tidy_MG <- tidy_MG %>%
  pivot_longer(cols = Citrobacter_freundii:Tatumella_saanichensis,
               names_to = "Species", values_to = "Fit")

tidy_MG_se <- subset(M_gav_pred, select = CF_SE:TS_SE)
tidy_MG_se$Mean <- M_gav_pred$Mean
tidy_MG_se <- tidy_MG_se %>%
  pivot_longer(cols = CF_SE:TS_SE,
               names_to = "Species", values_to = "Fit_SE")

unique(tidy_MG$Mean == tidy_MG_se$Mean)

M_gav_pred <- cbind(tidy_MG, subset(tidy_MG_se, select = Fit_SE))
rm(tidy_MG, tidy_MG_se)

M_gav_pred <- mutate(M_gav_pred,
                     upper = Fit + (Fit_SE * 1.96),
                     lower = Fit - (Fit_SE * 1.96))

MG_pred <- ggplot(M_gav_pred, aes(x = Mean, y = Fit, colour = Species)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = Mean), inherit.aes = FALSE, 
              data = M_gav_pred, alpha = 0.2) +
  geom_line() +
  scale_colour_manual(values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
                      labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis", 
                                 "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
  theme(legend.position = "bottom", text = element_text(size = 9),
        legend.text = element_text(face = "italic")) +
  labs(x = "Gene Length Mean", y = "Probability", 
       colour = "Species")

ggsave(MG_pred, file = "9_2Plots_gaviniae/MG_CR_pred.png", 
       width = 16.51, height = 12.38, units = "cm")

### Distances vs gene length ###
D_GLc <- readRDS("10Models/MC_Distance.rds")
D_GLg <- readRDS("10Models/MG_Distance.rds")

M_cal_DGL <- subset(M_calida_sdist_L, Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                                     "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                                     "Tatumella_ptyseos"))
D_GLc <- bam(Distance ~ s(Gene_Length, k = 42),
             data = M_cal_DGL, method = "fREML", control = ctrl, family = tw(link = "log"), discrete = TRUE)

saveRDS(D_GLc, file = "10Models/MC_Distance.rds")                         # Saves the model

layout(matrix(1:4, ncol = 2, byrow = TRUE))
gam.check(D_GLc)
draw(D_GLc)

M_gav_DGL <- subset(M_gaviniae_sdist_L, Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                                       "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                                       "Tatumella_ptyseos"))

D_GLg <- bam(Distance ~ s(Gene_Length, k = 42),
             data = M_gav_DGL, method = "fREML", control = ctrl, family = tw(link = "log"), discrete = TRUE)

saveRDS(D_GLg, file = "10Models/MG_Distance.rds")                         # Saves the model

layout(matrix(1:4, ncol = 2, byrow = TRUE))
gam.check(D_GLg)
draw(D_GLg)

## Predictions ##
M_cal_pred <- data.frame(Gene_Length = 150:4224)

M_cal_pred <- cbind(M_cal_pred, predict(D_GLc, newdata = M_cal_pred, se = TRUE, type = "response"))

M_cal_pred <- mutate(M_cal_pred,
                     upper = fit + (se.fit * 1.96),
                     lower = fit - (se.fit * 1.96))

MC_pred <- ggplot(M_cal_pred, aes(x = Gene_Length, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = Gene_Length), inherit.aes = FALSE, 
              data = M_cal_pred, alpha = 0.2) +
  geom_line() +
  theme(legend.position = "bottom", text = element_text(size = 9)) +
  labs(x = "Gene Length", y = "Distance")

ggsave(MC_pred, file = "9_1Plots_calida/MC_DGL_pred.png", 
       width = 16.51, height = 12.38, units = "cm")


M_gav_pred <- data.frame(Gene_Length = 150:4224)

M_gav_pred <- cbind(M_gav_pred, predict(D_GLg, newdata = M_gav_pred, se = TRUE, type = "response"))

M_gav_pred <- mutate(M_gav_pred,
                     upper = fit + (se.fit * 1.96),
                     lower = fit - (se.fit * 1.96))

MG_pred <- ggplot(M_gav_pred, aes(x = Gene_Length, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = Gene_Length), inherit.aes = FALSE, 
              data = M_gav_pred, alpha = 0.2) +
  geom_line() +
  theme(legend.position = "bottom", text = element_text(size = 9)) +
  labs(x = "Gene Length", y = "Distance")

ggsave(MG_pred, file = "9_2Plots_gaviniae/MG_DGL_pred.png", 
       width = 16.51, height = 12.38, units = "cm")

#
### Significant First Relative ############################################################################################################################
# This section does three things. (1) It identifies genes wherein the first non-Mixta relative is significantly different from the other species. (2) It
# identifies genes wherein - given that the first two non-Mixta species are of the same genus - the first genus is significant different from the next
# species. (3) It identifies genes wherein the first first two relatives are of the same genus.
Species_noM <- c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                 "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                 "Tatumella_ptyseos")

M_calida_sdist <- read.csv(file = "8Results/M_calida_Sort_Dist.csv",
                           stringsAsFactors = FALSE)

M_calida_sdist <- subset(M_calida_sdist, Species %in% Species_noM)        # Remove Mixta species

M_calida_sdist <- mutate(M_calida_sdist,
                         lower_conf = M_calida_sdist$Distance - (M_calida_sdist$Std_Errors * 1.96),
                         upper_conf = M_calida_sdist$Distance + (M_calida_sdist$Std_Errors * 1.96))


M_gaviniae_sdist <- read.csv(file = "8Results/M_gaviniae_Sort_Dist.csv",
                             stringsAsFactors = FALSE)

M_gaviniae_sdist <- subset(M_gaviniae_sdist, Species %in% Species_noM)    # Remove Mixta species

M_gaviniae_sdist <- mutate(M_gaviniae_sdist,
                           lower_conf = M_gaviniae_sdist$Distance - (M_gaviniae_sdist$Std_Errors * 1.96),
                           upper_conf = M_gaviniae_sdist$Distance + (M_gaviniae_sdist$Std_Errors * 1.96))


uniq_gene <- as.data.frame(unique(M_calida_sdist$Gene))                   # Dataframe with the 1035 gene names
colnames(uniq_gene) <- "Gene_name"

### Significant genes by species and genus ###
# M. calida #
sig_genes <- as.data.frame(matrix(ncol = 0, nrow = 0))
for(row in 1:nrow(uniq_gene)) {
  gene <- subset(M_calida_sdist, Gene == uniq_gene$Gene_name[row])
  
  x <- Rfast::nth(x = gene$Distance, k = 1, descending = FALSE, index.return = TRUE)
  y <- Rfast::nth(x = gene$Distance, k = 2, descending = FALSE, index.return = TRUE)
  z <- Rfast::nth(x = gene$Distance, k = 3, descending = FALSE, index.return = TRUE)
  
  genus_check <- case_when(
    x %in% c(1, 8) & y %in% c(1, 8) ~ TRUE,                               # Tatumella
    x %in% c(6, 7) & y %in% c(6, 7) ~ TRUE,                               # Pantoea
    x %in% c(4, 5) & y %in% c(4, 5) ~ TRUE,                               # Erwinia
    TRUE ~ FALSE
  )
  
  sig <- as.data.frame(matrix(ncol = 0, nrow = 1))
  sig <- mutate(sig,
                Gene_name = uniq_gene$Gene_name[row],
                First_Species = case_when(
                  gene$upper_conf[x] < gene$lower_conf[y] ~ TRUE,
                  TRUE ~ FALSE
                ),
                Two_Same_Genus = genus_check,
                First_Genus = case_when(
                  genus_check == TRUE ~ case_when(                        # If first two are same genus, compare second and third relatives
                    gene$upper_conf[y] < gene$lower_conf[z] ~ TRUE,       # Distances are significantly different
                    TRUE ~ FALSE                                          # Distances are not significantly different
                  ),
                  genus_check == FALSE ~ FALSE                            # Don't bother comparing since the first two relatives aren't the same genus
                ))
  
  sig_genes <- rbind(sig_genes, sig)
}
rm(gene, sig, genus_check, row, x, y, z)

write.csv(x = sig_genes, file = "8Results/M_calida_Sig_Rel.csv",
          row.names = FALSE)


# M.gaviniae #
sig_genes <- as.data.frame(matrix(ncol = 0, nrow = 0))
for(row in 1:nrow(uniq_gene)) {
  gene <- subset(M_gaviniae_sdist, Gene == uniq_gene$Gene_name[row])
  
  x <- Rfast::nth(x = gene$Distance, k = 1, descending = FALSE, index.return = TRUE)
  y <- Rfast::nth(x = gene$Distance, k = 2, descending = FALSE, index.return = TRUE)
  z <- Rfast::nth(x = gene$Distance, k = 3, descending = FALSE, index.return = TRUE)
  
  genus_check <- case_when(
    x %in% c(1, 8) & y %in% c(1, 8) ~ TRUE,                               # Tatumella
    x %in% c(6, 7) & y %in% c(6, 7) ~ TRUE,                               # Pantoea
    x %in% c(4, 5) & y %in% c(4, 5) ~ TRUE,                               # Erwinia
    TRUE ~ FALSE
  )
  
  sig <- as.data.frame(matrix(ncol = 0, nrow = 1))
  sig <- mutate(sig,
                Gene_name = uniq_gene$Gene_name[row],
                First_Species = case_when(
                  gene$upper_conf[x] < gene$lower_conf[y] ~ TRUE,
                  TRUE ~ FALSE
                ),
                Two_Same_Genus = genus_check,
                First_Genus = case_when(
                  genus_check == TRUE ~ case_when(                        # If first two are same genus, compare second and third relatives
                    gene$upper_conf[y] < gene$lower_conf[z] ~ TRUE,       # Distances are significantly different
                    TRUE ~ FALSE                                          # Distances are not significantly different
                  ),
                  genus_check == FALSE ~ FALSE                            # Don't bother comparing since the first two relatives aren't the same genus
                ))
  
  sig_genes <- rbind(sig_genes, sig)
}
rm(gene, sig, genus_check, row, x, y, z)

write.csv(x = sig_genes, file = "8Results/M_gaviniae_Sig_Rel.csv",
          row.names = FALSE)

### Usual Relative Order ##################################################################################################################################
M_calida_L <- read.csv(file = "8Results/M_calida_Relatives_Length.csv", 
                       stringsAsFactors = FALSE)
M_gaviniae_L <- read.csv(file = "8Results/M_gaviniae_Relatives_Length.csv",
                         stringsAsFactors = FALSE)

patt_MC <- as.data.frame(unique(M_calida_L$Rela_Pattern))
colnames(patt_MC) <- "Rel_Pattern"

M_cal_patterns <- as.data.frame(matrix(nrow = 0, ncol = 0))
for(row in 1:nrow(patt_MC)) {
  pat <- as.data.frame(matrix(nrow = 1, ncol = 0))
  pat <- mutate(pat,
                Relative_Pattern = patt_MC$Rel_Pattern[row],
                Number = length(which(M_calida_L$Rela_Pattern == patt_MC$Rel_Pattern[row])))
  
  M_cal_patterns <- rbind(M_cal_patterns, pat)
}
write.csv(x = M_cal_patterns, file = "8Results/M_calida_Relative_Pattern.csv", row.names = FALSE)

patt_MG <- as.data.frame(unique(M_gaviniae_L$Rela_Pattern))
colnames(patt_MG) <- "Rel_Pattern"

M_gav_patterns <- as.data.frame(matrix(nrow = 0, ncol = 0))
for(row in 1:nrow(patt_MG)) {
  pat <- as.data.frame(matrix(nrow = 1, ncol = 0))
  pat <- mutate(pat,
                Relative_Pattern = patt_MG$Rel_Pattern[row],
                Number = length(which(M_gaviniae_L$Rela_Pattern == patt_MG$Rel_Pattern[row])))
  
  M_gav_patterns <- rbind(M_gav_patterns, pat)
}
write.csv(x = M_gav_patterns, file = "8Results/M_gaviniae_Relative_Pattern.csv", row.names = FALSE)
#
### Nucleotide numbers ####################################################################################################################################
### Mixta calida ###
M_cal <- read.csv(file = "8Results/M_calida_Distances.csv", stringsAsFactors = FALSE)
M_cal <- subset(M_cal, select = ID:Gene)

M_cal$Gene_name <- substring(M_cal$Gene, 7)

fastaFiles <- as.data.frame(list.files(path = "3Homologues_10/",
                                       pattern = ".fasta"))               # Makes a dataframe where the first column is a list of fasta gene files
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/3Homologues_10/", 
                              fastaFiles$File_name, sep = "")             # Creates a file pathway for each gene

fastaFiles <- mutate(fastaFiles,
                     Gene = gsub(pattern = ".fasta", replacement = "", x = fastaFiles$File_name),
                     Check = case_when(
                       Gene %in% M_cal$Gene ~ TRUE,
                       TRUE ~ FALSE
                     ))

fastaFiles <- subset(fastaFiles, Check == TRUE, select = File_name:Gene)

nucl_calida <- as.data.frame(matrix(ncol = 0, nrow = 0))
for(row in 1:nrow(fastaFiles)) {
  path <- fastaFiles$Path_name[row]
  file <- as.character(fastaFiles$File_name[row])
  
  gene_file <- read.table(file = path, header = FALSE, sep = "\n", stringsAsFactors = FALSE)
  
  calida <- as.data.frame(gene_file[11, ])
  colnames(calida) <- "nucleotides"
  
  calida$gene <- fastaFiles$Gene[row]
  
  nucl_calida <- rbind(nucl_calida, calida)
}
rm(gene_file, calida, file, path, row)

nucl_calida$ID <- M_cal$ID
nucl_calida$nucleotides <- as.character(nucl_calida$nucleotides)

nucl_calida <- separate(data = nucl_calida, col = nucleotides, into = paste("V", 1:8, sep = ""), sep = ":", remove = FALSE, extra = "merge")
nucl_calida <- subset(nucl_calida, select = c(nucleotides, V3, gene, ID))

nucl_calida <- separate(data = nucl_calida, col = V3, into = c("Beg", "End"), sep = "-", remove = TRUE)
nucl_calida <- mutate(nucl_calida, 
                      Beg = as.numeric(Beg),
                      End = as.numeric(End))

nucl_calida <- nucl_calida[order(nucl_calida$ID),]

nucl_calida <- subset(nucl_calida, select = Beg:ID)

write.csv(x = nucl_calida, file = "8Results/M_calida_Nucleotide.csv", row.names = FALSE)


### Mixta gaviniae ###
M_gav <- read.csv(file = "8Results/M_gaviniae_Distances.csv", stringsAsFactors = FALSE)
M_gav <- subset(M_gav, select = ID:Gene)

M_gav$Gene_name <- substring(M_gav$Gene, 7)

fastaFiles <- as.data.frame(list.files(path = "3Homologues_10/",
                                       pattern = ".fasta"))               # Makes a dataframe where the first column is a list of fasta gene files
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/3Homologues_10/", 
                              fastaFiles$File_name, sep = "")             # Creates a file pathway for each gene

fastaFiles <- mutate(fastaFiles,
                     Gene = gsub(pattern = ".fasta", replacement = "", x = fastaFiles$File_name),
                     Check = case_when(
                       Gene %in% M_gav$Gene ~ TRUE,
                       TRUE ~ FALSE
                     ))

fastaFiles <- subset(fastaFiles, Check == TRUE, select = File_name:Gene)

nucl_gaviniae <- as.data.frame(matrix(ncol = 0, nrow = 0))
for(row in 1:nrow(fastaFiles)) {
  path <- fastaFiles$Path_name[row]
  file <- as.character(fastaFiles$File_name[row])
  
  gene_file <- read.table(file = path, header = FALSE, sep = "\n", stringsAsFactors = FALSE)
  
  gaviniae <- as.data.frame(gene_file[13, ])
  colnames(gaviniae) <- "nucleotides"
  
  gaviniae$gene <- fastaFiles$Gene[row]
  
  nucl_gaviniae <- rbind(nucl_gaviniae, gaviniae)
}
rm(gene_file, gaviniae, file, path, row)

nucl_gaviniae$ID <- M_gav$ID
nucl_gaviniae$nucleotides <- as.character(nucl_gaviniae$nucleotides)

nucl_gaviniae <- separate(data = nucl_gaviniae, col = nucleotides, into = paste("V", 1:8, sep = ""), sep = ":", remove = FALSE, extra = "merge")
nucl_gaviniae <- subset(nucl_gaviniae, select = c(nucleotides, V3, gene, ID))

nucl_gaviniae <- separate(data = nucl_gaviniae, col = V3, into = c("Beg", "End"), sep = "-", remove = TRUE)
nucl_gaviniae <- mutate(nucl_gaviniae, 
                        Beg = as.numeric(Beg),
                        End = as.numeric(End))

nucl_gaviniae <- nucl_gaviniae[order(nucl_gaviniae$ID),]

nucl_gaviniae <- subset(nucl_gaviniae, select = Beg:ID)

write.csv(x = nucl_gaviniae, file = "8Results/M_gaviniae_Nucleotide.csv", row.names = FALSE)
#
### Genes for MEGAX analysis ##############################################################################################################################
BM <- read.csv(file = "8Results/Best_Model.csv", stringsAsFactors = FALSE)

M_cal <- read.csv(file = "8Results/M_calida_Distances.csv", stringsAsFactors = FALSE)
M_cal <- subset(M_cal, select = ID:Gene)

M_gav <- read.csv(file = "8Results/M_gaviniae_Distances.csv", stringsAsFactors = FALSE)
M_gav <- subset(M_gav, select = ID:Gene)

M_cal_RL <- read.csv(file = "8Results/M_calida_Relatives_Length.csv", stringsAsFactors = FALSE)
M_cal_RL <- subset(M_cal_RL, select = c(Gene, Mean, Rela_Pattern))

M_gav_RL <- read.csv(file = "8Results/M_gaviniae_Relatives_Length.csv", stringsAsFactors = FALSE)
M_gav_RL <- subset(M_gav_RL, select = c(Gene, Mean, Rela_Pattern))

unique(M_cal$Gene == c(M_cal_RL$Gene, M_gav$Gene, M_gav_RL$Gene)) # If TRUE, then continue

rlvnt_genes <- data.frame(matrix(ncol = 0, nrow = 1035))
rlvnt_genes <- mutate(rlvnt_genes,
                      Gene = M_cal$Gene,
                      cal_ID = M_cal$ID,
                      gav_ID = M_gav$ID,
                      Mean_gene_length = M_cal_RL$Mean,
                      cal_rel = M_cal_RL$Rela_Pattern, 
                      gav_rel = M_gav_RL$Rela_Pattern)

rm(M_cal, M_cal_RL, M_gav, M_gav_RL)

rlvnt_genes$cal_rel <- substring(rlvnt_genes$cal_rel, 7)
rlvnt_genes$gav_rel <- substring(rlvnt_genes$gav_rel, 7)

rlvnt_genes$Same_rel_patt <- case_when(
  rlvnt_genes$cal_rel == rlvnt_genes$gav_rel ~ TRUE,
  TRUE ~ FALSE
)

write.csv(x = rlvnt_genes, file = "8Results/MEGAX_Relevant_Genes.csv", row.names = FALSE)
#
### MEGAX: rRNA genes #####################################################################################################################################
rel_rRNA_TIF <- data.frame(Gene = c("38003_rsmG", "38735_rsmJ", "39935_rsmD", "38555_rsmB", "39397_rsmA", "39418_rsmH", "40004_rpsL", "40005_rpsG", 
                                    "38586_rpsJ", "38581_rpsS", "38579_rpsC", "38576_rpsQ", "38572_rpsN", "38571_rpsH", "38568_rpsE", "38561_rpsD", 
                                    "38441_rpsI", "38498_rpsO", "37935_rpsB", "38812_rpsA", "40035_rpsP", "37746_rpsU", "37893_rpsR", "37891_rpsF",
                                    "38585_rplC", "38584_rplD", "38583_rplW", "38582_rplB", "38580_rplV", "38578_rplP", "38577_rpmC", "38575_rplN",
                                    "38574_rplX", "38573_rplE", "38570_rplF", "38567_rpmD", "38566_rplO", "38559_rplQ", "38413_prmA", "38440_rplM", 
                                    "38470_rpmA", "38613_rpmF", "38895_rpmI", "38896_rplT", "40096_rplY", "37894_rplI", "40394_rplL", "40393_rplJ", 
                                    "40392_rplA", "40391_rplK", "39858_rpmE", "39815_rpmG", "39814_rpmB", "38495_infB"),
                           Product = c(rep("16S rRNA guanine methyltransferase", 3), "16S rRNA cytosine methyltransferase", 
                                       "16S rRNA adenine dimethyltransferase", "16S rRNA cytosine methyltransferase", rep("30S ribosomal protein", 18), 
                                       rep("50S ribosomal protein", 14), "50S ribosomal protein L11 methyltransferase", rep("50S ribosomal protein", 14),
                                       "Translation initiation factor"))
rel_rRNA_TIF <- rel_rRNA_TIF[order(rel_rRNA_TIF$Gene), ]

rlvnt_genes$rRNA_TIF <- case_when(
  rlvnt_genes$Gene %in% rel_rRNA_TIF$Gene ~ TRUE,
  TRUE ~ FALSE
)

rlvnt_rRNA_TIF <- subset(rlvnt_genes, rRNA_TIF == TRUE, select = Gene:Same_rel_patt)

unique(rlvnt_rRNA_TIF$Gene == rel_rRNA_TIF$Gene) #If TRUE, then continue
rlvnt_rRNA_TIF$Product <- rel_rRNA_TIF$Product

write.csv(x = rlvnt_rRNA_TIF, file = "8Results/rRNA-TIF_Genes.csv", row.names = FALSE)

BM$rRNA <- case_when(
  BM$Gene %in% rlvnt_rRNA_TIF$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_rRNA <- subset(BM, rRNA == TRUE, select = Gene:ModelCode)

write.csv(x = BM_rRNA, file = "8Results/Best_Model_rRNA-TIF.csv", row.names = FALSE)

### MEGAX: MLSA genes ######################################################################################################################################
rel_MLSA <- data.frame(Gene = c("38009_atpA", "38005_atpB", "38012_atpC", "38006_atpE", "38007_atpF", "38010_atpG", "38008_atpH", "38004_atpI", 
                                "39382_dnaK", "39185_leuS", "38495_infB", "40044_recA", "39904_recF", "38360_recN", "38375_recO", "37521_recR", 
                                "40043_recX", "38560_rpoA", "40395_rpoB", "40396_rpoC", "37744_rpoD", "39939_rpoH", "37484_rpoZ"),
                       Product = c("F0F1 ATP synthase", "F0F1 ATP synthase", "F0F1 ATP synthase", "F0F1 ATP synthase", "F0F1 ATP synthase", 
                                   "F0F1 ATP synthase", "F0F1 ATP synthase", "F0F1 ATP synthase", 
                                   
                                   "Molecular chaperone", 
                                   
                                   "leucine--tRNA ligase",
                                   
                                   "Translation initiation factor", 
                                   
                                   "Recombinase", "DNA replication/repair protein", "DNA repair protein", "DNA repair protein", 
                                   "Recombination protein", "Regulatory protein", 
                                   
                                   "DNA-directed RNA polymerase", "DNA-directed RNA polymerase", "DNA-directed RNA polymerase", 
                                   "RNA polymerase sigma factor", "RNA polymerase sigma factor", "DNA-directed RNA polymerase"))
rel_MLSA <- rel_MLSA[order(rel_MLSA$Gene), ]

rlvnt_genes$MLSA <- case_when(
  rlvnt_genes$Gene %in% rel_MLSA$Gene ~ TRUE,
  TRUE ~ FALSE
)

rlvnt_MLSA <- subset(rlvnt_genes, MLSA == TRUE, select = Gene:Same_rel_patt)

unique(rlvnt_MLSA$Gene == rel_MLSA$Gene) # If TRUE, then continue
rlvnt_MLSA$Product <- rel_MLSA$Product

write.csv(x = rlvnt_MLSA, file = "8Results/MLSA_Genes.csv", row.names = FALSE)

BM$MLSA <- case_when(
  BM$Gene %in% rlvnt_MLSA$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_MLSA <- subset(BM, MLSA == TRUE, select = Gene:ModelCode)

write.csv(x = BM_MLSA, file = "8Results/Best_Model_MLSA.csv", row.names = FALSE)

### MEGAX: Genera ##########################################################################################################################################
gnr <- subset(rlvnt_genes, select = Gene:Same_rel_patt)
gnr$first_cal <- substring(gnr$cal_rel, first = 1, last = 5)
gnr$first_gav <- substring(gnr$gav_rel, first = 1, last = 5)
gnr$same_first <- case_when(
  gnr$first_cal %in% c("CF_EC", "EC_CF") & gnr$first_gav %in% c("CF_EC", "EC_CF") ~ "Outgroup",
  gnr$first_cal %in% c("PA_PS", "PS_PA") & gnr$first_gav %in% c("PA_PS", "PS_PA") ~ "Pantoea",
  gnr$first_cal %in% c("EA_ET", "ET_EA") & gnr$first_gav %in% c("EA_ET", "ET_EA") ~ "Erwinia",
  gnr$first_cal %in% c("TP_TS", "TS_TP") & gnr$first_gav %in% c("TP_TS", "TS_TP") ~ "Tatumella",
  TRUE ~ "Nope"
)

gnr_same <- subset(gnr, same_first != "Nope")

EA_ET <- subset(gnr_same, first_cal %in% c("EA_ET", "ET_EA"))
EA_ET$random <- sample(rep_len(1:(nrow(EA_ET)/5), nrow(EA_ET))) # NOTE: will give different results each time this line is run
EA_ET <- subset(EA_ET, random == round(nrow(EA_ET)/5))

PA_PS <- subset(gnr_same, first_cal %in% c("PA_PS", "PS_PA"))
PA_PS$random <- sample(rep_len(1:(nrow(PA_PS)/5), nrow(PA_PS)))
PA_PS <- subset(PA_PS, random == round(nrow(PA_PS)/5))

TP_TS <- subset(gnr_same, first_cal %in% c("TP_TS", "TS_TP"))

EC_CF <- subset(gnr_same, first_cal %in% c("EC_CF", "CF_EC"))
EC_CF$random <- sample(rep_len(1:(nrow(EC_CF)/5), nrow(EC_CF)))
EC_CF <- subset(EC_CF, random == round(nrow(EC_CF)/5))

Genes <- data.frame(Gene = c("37985_ilvA", "38387_mltF_2", "39108_moaE", "39192_mrdB", "39776_thiL", "37400_yffB", "37940_uppS", "38441_rpsI", 
                             "38681_ppiC", "39090_rlmF", "37936_tsf", "38581_rpsS", "37335_iscS"),
                    Product = c("Threonine ammonia-lyase, biosynthetic", "Membrane-bound lytic murein transglycosylase",
                                "Molybdopterin synthase catalytic subunit", "Peptidoglycan glycosyltransferase", "Thiamine-phosphate kinase", 
                                "ArsC family reductase", "(2E,6E)-farnesyl-diphosphate-specific ditrans, polycis-undecaprenyl-diphosphate synthase",
                                "30S ribosomal protein", "peptidylprolyl isomerase", "23S rRNA (adenine(1618)-N(6))-methyltransferase",
                                "Elongation factor Ts", "30S ribosomal protein S19", "IscS subfamily cysteine desulfurase"),
                    Genus = c(rep("Pantoea", 5), rep("Erwinia", 5), rep("Tatumella", 2), "Outgroup"))

Genes <- Genes[order(Genes$Gene), ]

rlvnt_genes$genera <- case_when(
  rlvnt_genes$Gene %in% Genes$Gene ~ TRUE,
  TRUE ~ FALSE
)

rlvnt_Genera <- subset(rlvnt_genes, genera == TRUE, select = Gene:Same_rel_patt)

unique(rlvnt_Genera$Gene == Genes$Gene) # If TRUE, then continue
rlvnt_Genera$Product <- Genes$Product

write.csv(x = rlvnt_Genera, file = "8Results/Genera_Genes.csv", row.names = FALSE)

BM$Genera <- case_when(
  BM$Gene %in% rlvnt_Genera$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_Genera <- subset(BM, Genera == TRUE, select = Gene:ModelCode)

write.csv(x = BM_Genera, file = "8Results/Best_Model_Genera.csv", row.names = FALSE)

### MEGAX: First species is significantly different from next ##############################################################################################
M_cal_SR <- read.csv(file = "8Results/M_calida_Sig_Rel.csv", stringsAsFactors = FALSE)
M_cal_SR <- subset(M_cal_SR, First_Species == TRUE)
M_cal_SR <- M_cal_SR[order(M_cal_SR$Gene_name), ]

M_cal_RL <- read.csv(file = "8Results/M_calida_Relatives_Length.csv", stringsAsFactors = FALSE)
M_cal_RL <- subset(M_cal_RL, select = c(Gene, One, Two, Three, Four, Mean, Rela_Pattern))

M_cal_RL$Same_Species <- case_when(
  M_cal_RL$Gene %in% M_cal_SR$Gene_name ~ TRUE,
  TRUE ~ FALSE
)
M_cal_RL <- subset(M_cal_RL, Same_Species == TRUE)
M_cal_RL$first <- substring(M_cal_RL$Three, first = 1, last = 7)
M_cal_RL$second <- substring(M_cal_RL$Four, first = 1, last = 7)
M_cal_RL$Same <- case_when(
  M_cal_RL$first == M_cal_RL$second ~ TRUE,
  TRUE ~ FALSE
)
M_cal_RL <- subset(M_cal_RL, Same == TRUE)

rlvnt_genes$Sig_Species <- case_when(
  rlvnt_genes$Gene %in% M_cal_RL$Gene ~ TRUE,
  TRUE ~ FALSE
)

rlvnt_Sig_Species <- subset(rlvnt_genes, Sig_Species == TRUE, select = Gene:Same_rel_patt)

rlvnt_Sig_Species <- mutate(rlvnt_Sig_Species,
                            Product = c("Histidine--tRNA ligase", "Exopolyphosphatase", "Glucokinase", "Miniconductance mechanosensitive channel",
                                        "DNA topoisomerase IV", "Phosphoglucosamine mutase", "Superoxide dismutase [Mn]", "MFS transporter",
                                        "Leucine--tRNA ligase", "Peptidoglycan glycosyltransferase", "2-isopropylmalate synthase", 
                                        "UDP-3-O-acyl-N-acetylglucosamine deacetylase", "Pyridoxal phosphate-dependent aminotransferase", 
                                        "NADH-quinone oxidoreductase subunit", "NADH-quinone oxidoreductase subunit", 
                                        "NADH-quinone oxidoreductase subunit", "Ribonucleoside-diphosphate reductase subunit", "YdbH family protein",
                                        "Aspartate-semialdehyde dehydrogenase", "Two-component system sensor histidine kinase", 
                                        "Signal recognition particle protein", 
                                        "Trifunctional transcriptional regulator/proline dehydrogenase/L-glutamate gamma-semialdehyde dehydrogenase",
                                        "Carboxy terminal-processing peptidase"),
                            Gene_check = c("37350_hisS", "37381_ppx", "37443_glk", "37870_mscM", "38070_parC", "38491_glmM", "38689_sodA", "38820_ycaD_2",
                                           "39185_leuS", "39192_mrdB", "39413_leuA", "39431_lpxC", "39534_alaA", "39540_nuoF", "39541_nuoG", "39547_nuoM",
                                           "39557_nrdA", "39682_hypothetical_protein", "39951_asd", "39964_envZ", "40036_ffh", "40248_putA", "40292_prc"))
unique(rlvnt_Sig_Species$Gene == rlvnt_Sig_Species$Gene_check) # If TRUE, then continue
rlvnt_Sig_Species <- subset(rlvnt_Sig_Species, select = Gene:Product)

write.csv(x = rlvnt_Sig_Species, file = "8Results/Sig_Species_Genes_cal.csv", row.names = FALSE)

BM$Sig_Species <- case_when(
  BM$Gene %in% rlvnt_Sig_Species$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_Sig_Species <- subset(BM, Sig_Species == TRUE, select = Gene:ModelCode)

write.csv(x = BM_Sig_Species, file = "8Results/Best_Model_Sig_Species_cal.csv", row.names = FALSE)


# M. gaviniae #
M_gav_SR <- read.csv(file = "8Results/M_gaviniae_Sig_Rel.csv", stringsAsFactors = FALSE)
M_gav_SR <- subset(M_gav_SR, First_Species == TRUE)
M_gav_SR <- M_gav_SR[order(M_gav_SR$Gene_name), ]

M_gav_RL <- read.csv(file = "8Results/M_gaviniae_Relatives_Length.csv", stringsAsFactors = FALSE)
M_gav_RL <- subset(M_gav_RL, select = c(Gene, One, Two, Three, Four, Mean, Rela_Pattern))

M_gav_RL$Same_Species <- case_when(
  M_gav_RL$Gene %in% M_gav_SR$Gene_name ~ TRUE,
  TRUE ~ FALSE
)
M_gav_RL <- subset(M_gav_RL, Same_Species == TRUE)

M_gav_RL$first <- substring(M_gav_RL$Three, first = 1, last = 7)
M_gav_RL$second <- substring(M_gav_RL$Four, first = 1, last = 7)
M_gav_RL$Same <- case_when(
  M_gav_RL$first == M_gav_RL$second ~ TRUE,
  TRUE ~ FALSE
)
M_gav_RL <- subset(M_gav_RL, Same == TRUE)

rlvnt_genes$Sig_Species_gav <- case_when(
  rlvnt_genes$Gene %in% M_gav_RL$Gene ~ TRUE,
  TRUE ~ FALSE
)

rlvnt_Sig_Species_g <- subset(rlvnt_genes, Sig_Species_gav == TRUE, select = Gene:Same_rel_patt)

rlvnt_Sig_Species_g <- mutate(rlvnt_Sig_Species_g,
                              Product = c("Exopolyphosphatase", "Miniconductance mechanosensitive channel", "DNA topoisomerase IV subunit", 
                                          "Superoxide dismutase [Mn]","23S rRNA (guanine(2445)-N(2))-methyltransferase", 
                                          "Replication-associated recombination protein A", "Glutamate 5-kinase", "Class II glutamine amidotransferase", 
                                          "Glycine betaine/L-proline ABC transporter ATP-binding protein", "Phenylalanine--tRNA ligase subunit",
                                          "Leucine--tRNA ligase", "Peptidoglycan glycosyltransferase", "2-isopropylmalate synthase", 
                                          "Undecaprenyldiphospho-muramoylpentapeptide beta-N-acetylglucosaminyltransferase",
                                          "Bifunctional aconitate hydratase 2/2-methylisocitrate dehydratase", "ABC transporter permease",
                                          "Ribonucleoside-diphosphate reductase subunit", "Ferredoxin--NADP(+) reductase", "Porphobilinogen synthase", 
                                          "HlyC/CorC family transporter", "Histidinol dehydrogenase", "ATP phosphoribosyltransferase",
                                          "Trifunctional transcriptional regulator/proline dehydrogenase/L-glutamate gamma-semialdehyde dehydrogenase"),
                              Gene_check = c("37381_ppx", "37870_mscM",
                                             "38070_parC", "38689_sodA",
                                             "38788_rlmL", "38824_rarA",
                                             "38841_proB", "38851_yafJ",
                                             "38877_proV", "38898_pheT",
                                             "39185_leuS", "39192_mrdB",
                                             "39413_leuA", "39426_murG",
                                             "39453_acnB", "39527_hisM",
                                             "39557_nrdA", "39847_fpr",
                                             "39927_hemB", "40038_hypothetical_protein",
                                             "40162_hisD", "40163_hisG",
                                             "40248_putA"))

unique(rlvnt_Sig_Species_g$Gene == rlvnt_Sig_Species_g$Gene_check) # If TRUE, then continue
rlvnt_Sig_Species_g <- subset(rlvnt_Sig_Species_g, select = Gene:Product)

write.csv(x = rlvnt_Sig_Species_g, file = "8Results/Sig_Species_Genes_gav.csv", row.names = FALSE)

BM$Sig_Species <- case_when(
  BM$Gene %in% rlvnt_Sig_Species_g$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_Sig_Species <- subset(BM, Sig_Species == TRUE, select = Gene:ModelCode)

write.csv(x = BM_Sig_Species, file = "8Results/Best_Model_Sig_Species_gav.csv", row.names = FALSE)


test <- read.csv(file = "8Results/Sig_Species_Genes_gav.csv", stringsAsFactors = FALSE)
test2 <- paste(test$Gene, ".fasta", sep = ""); test2


nucl_g <- read.csv(file = "8Results/M_gaviniae_Nucleotide.csv", stringsAsFactors = FALSE)
nucl_g$Gene_name <- substring(nucl_g$gene, 7)
test <- M_gav_RL$Gene

### MEGAX: First genus is significantly different from next ###############################################################################################
rlvnt_genes <- read.csv(file = "8Results/MEGAX_Relevant_Genes.csv", stringsAsFactors = FALSE)
BM <- read.csv(file = "8Results/Best_Model.csv", stringsAsFactors = FALSE)

M_cal_SR <- read.csv(file = "8Results/M_calida_Sig_Rel.csv", 
                     stringsAsFactors = FALSE)
M_cal_SR <- subset(M_cal_SR, First_Genus == TRUE)   
M_cal_SR <- M_cal_SR[order(M_cal_SR$Gene_name), ]

M_cal_SR$weird <- case_when(
  M_cal_SR$Gene_name %in% c("38235_csdA", "38434_hypothetical_protein", "38746_fhuA_1", "38985_ttuB", "39549_hypothetical_protein") ~ TRUE,
  TRUE ~ FALSE
)

M_cal_SR <- subset(M_cal_SR, weird == FALSE, select = Gene_name:First_Genus)

M_cal_RL <- read.csv(file = "8Results/M_calida_Relatives_Length.csv", stringsAsFactors = FALSE)
M_cal_RL <- subset(M_cal_RL, select = c(Gene, One, Two, Three, Four, Mean, Rela_Pattern))

M_cal_RL$Same_Genus <- case_when(
  M_cal_RL$Gene %in% M_cal_SR$Gene_name ~ TRUE,
  TRUE ~ FALSE
)
M_cal_RL <- subset(M_cal_RL, Same_Genus == TRUE)

rlvnt_genes$Sig_Genus <- case_when(
  rlvnt_genes$Gene %in% M_cal_RL$Gene ~ TRUE,
  TRUE ~ FALSE
)

rlvnt_Sig_Genus <- subset(rlvnt_genes, Sig_Genus == TRUE, select = Gene:Same_rel_patt)

rlvnt_Sig_Genus <- mutate(rlvnt_Sig_Genus,
                          Product = c("Membrane-bound lytic murein transglycosylase", "Protein-disulfide reductase", "dGTPase", 
                                      "ADP-ribose diphosphatase", "Unknown*", "Murein L,D-transpeptidase","YdgA family protein", 
                                      "DUF1615 domain-containing protein", "PTS transporter subunit", "LPS assembly lipoprotein", 
                                      "Proline-specific permease", "Branched-chain amino acid transporter carrier protein", "EamA family transporter", 
                                      "NCS2 family permease"),
                          Gene_check = c("37594_emtA", "37857_dsbD", "37929_dgt", "38062_nudF", "38683_pdeK", "38850_hypothetical_protein", "38936_ydgA", 
                                         "39107_hypothetical_protein", "39169_nagE", "39186_lptE", "39562_proY", "39563_brnQ", "39916_eamA", 
                                         "40196_adeP"))

unique(rlvnt_Sig_Genus$Gene == rlvnt_Sig_Genus$Gene_check) # If TRUE, then continue
rlvnt_Sig_Genus <- subset(rlvnt_Sig_Genus, select = Gene:Product)

write.csv(x = rlvnt_Sig_Genus, file = "8Results/Sig_Genus_Genes_cal.csv", row.names = FALSE)

BM$Sig_Genus <- case_when(
  BM$Gene %in% rlvnt_Sig_Genus$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_Sig_Genus <- subset(BM, Sig_Genus == TRUE, select = Gene:ModelCode)

write.csv(x = BM_Sig_Genus, file = "8Results/Best_Model_Sig_Genus_cal.csv", row.names = FALSE)


# M. gaviniae #
M_gav_SR <- read.csv(file = "8Results/M_gaviniae_Sig_Rel.csv", 
                     stringsAsFactors = FALSE)
M_gav_SR <- subset(M_gav_SR, First_Genus == TRUE)   
M_gav_SR <- M_gav_SR[order(M_gav_SR$Gene_name), ]

M_gav_SR$weird <- case_when(
  M_gav_SR$Gene_name %in% c("38235_csdA", "38434_hypothetical_protein", "38746_fhuA_1", "38985_ttuB", "39549_hypothetical_protein") ~ TRUE,
  TRUE ~ FALSE
)

M_gav_SR <- subset(M_gav_SR, weird == FALSE, select = Gene_name:First_Genus)

M_gav_RL <- read.csv(file = "8Results/M_gaviniae_Relatives_Length.csv", stringsAsFactors = FALSE)
M_gav_RL <- subset(M_gav_RL, select = c(Gene, One, Two, Three, Four, Mean, Rela_Pattern))

M_gav_RL$Same_Genus <- case_when(
  M_gav_RL$Gene %in% M_gav_SR$Gene_name ~ TRUE,
  TRUE ~ FALSE
)
M_gav_RL <- subset(M_gav_RL, Same_Genus == TRUE)

rlvnt_genes$Sig_Genus <- case_when(
  rlvnt_genes$Gene %in% M_gav_RL$Gene ~ TRUE,
  TRUE ~ FALSE
)

rlvnt_Sig_Genus <- subset(rlvnt_genes, Sig_Genus == TRUE, select = Gene:Same_rel_patt)

rlvnt_Sig_Genus <- mutate(rlvnt_Sig_Genus,
                          Product = c("YchO/YchP family invasin", "dGTPase", "ADP-ribose diphosphatase", "Lipoprotein", "Unknown*", 
                                      "DUF1615 domain-containing protein", "PTS transporter subunit", "LPS assembly protein", 
                                      "Type II secretion system protein", "Proline-specific permease", 
                                      "Branched-chain amino acid transporter carrier protein", "Phosphate regulon sensor histidine kinase",
                                      "EamA family transporter", "NCS2 family permease", "Carboxy terminal-processing peptidase"),
                          Gene_check = c("37315_hypothetical_protein", "37929_dgt", "38062_nudF", "38500_nlpI", "38683_pdeK", "39107_hypothetical_protein",
                                         "39169_nagE", "39400_lptD", "39439_gspE", "39562_proY", "39563_brnQ", "39565_phoR", "39916_eamA", "40196_adeP",
                                         "40292_prc"))

unique(rlvnt_Sig_Genus$Gene == rlvnt_Sig_Genus$Gene_check) # If TRUE, then continue
rlvnt_Sig_Genus <- subset(rlvnt_Sig_Genus, select = Gene:Product)

write.csv(x = rlvnt_Sig_Genus, file = "8Results/Sig_Genus_Genes_gav.csv", row.names = FALSE)

BM$Sig_Genus <- case_when(
  BM$Gene %in% rlvnt_Sig_Genus$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_Sig_Genus <- subset(BM, Sig_Genus == TRUE, select = Gene:ModelCode)

write.csv(x = BM_Sig_Genus, file = "8Results/Best_Model_Sig_Genus_gav.csv", row.names = FALSE)

#
### MEGAX: First two relatives are NOT of the same genus ###################################################################################################
# Somehow, somewehere, I screwed up which genes had first relatives being of the same genus and which didn't

rlvnt_genes <- read.csv(file = "8Results/MEGAX_Relevant_Genes.csv", stringsAsFactors = FALSE)
BM <- read.csv(file = "8Results/Best_Model.csv", stringsAsFactors = FALSE)

rlvnt_Same_Genus <- rlvnt_genes
rlvnt_Same_Genus$cal_two <- substring(rlvnt_Same_Genus$cal_rel, first = 1, last = 5)
rlvnt_Same_Genus$gav_two <- substring(rlvnt_Same_Genus$gav_rel, first = 1, last = 5)

rlvnt_Same_Genus$cal_check <- case_when(
  rlvnt_Same_Genus$cal_two %in% c("PS_PA", "PA_PS", "ET_EA", "EA_ET", "TS_TP", "TP_TS", "EC_CF", "CF_EC") ~ FALSE,
  TRUE ~ TRUE
)

rlvnt_Same_Genus <- subset(rlvnt_Same_Genus, cal_check == TRUE)

rlvnt_Same_Genus$gav_check <- case_when(
  rlvnt_Same_Genus$gav_two %in% c("PS_PA", "PA_PS", "ET_EA", "EA_ET", "TS_TP", "TP_TS", "EC_CF", "CF_EC") ~ FALSE,
  TRUE ~ TRUE
)

rlvnt_Same_Genus <- subset(rlvnt_Same_Genus, gav_check == TRUE, select = Gene:gav_two)

rlvnt_Same_Genus$Same <- case_when(
  rlvnt_Same_Genus$cal_two %in% c("CF_ET", "ET_CF") & rlvnt_Same_Genus$gav_two %in% c("CF_ET", "ET_CF") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("CF_PA", "PA_CF") & rlvnt_Same_Genus$gav_two %in% c("CF_PA", "PA_CF") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("CF_PS", "PS_CF") & rlvnt_Same_Genus$gav_two %in% c("CF_PS", "PS_CF") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("CF_TS", "TS_CF") & rlvnt_Same_Genus$gav_two %in% c("CF_TS", "TS_CF") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("EA_PA", "PA_EA") & rlvnt_Same_Genus$gav_two %in% c("EA_PA", "PA_EA") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("EA_PS", "PS_EA") & rlvnt_Same_Genus$gav_two %in% c("EA_PS", "PS_EA") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("EC_EA", "EA_EC") & rlvnt_Same_Genus$gav_two %in% c("EC_EA", "EA_EC") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("EC_ET", "ET_EC") & rlvnt_Same_Genus$gav_two %in% c("EC_ET", "ET_EC") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("EC_PS", "PS_EC") & rlvnt_Same_Genus$gav_two %in% c("EC_PS", "PS_EC") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("ET_PA", "PA_ET") & rlvnt_Same_Genus$gav_two %in% c("ET_PA", "PA_ET") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("ET_PS", "PS_ET") & rlvnt_Same_Genus$gav_two %in% c("ET_PS", "PS_ET") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("PA_EC", "EC_PA") & rlvnt_Same_Genus$gav_two %in% c("PA_EC", "EC_PA") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("PS_TP", "TP_PS") & rlvnt_Same_Genus$gav_two %in% c("PS_TP", "TP_PS") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("PS_TS", "TS_PS") & rlvnt_Same_Genus$gav_two %in% c("PS_TS", "TS_PS") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("TP_EA", "EA_TP") & rlvnt_Same_Genus$gav_two %in% c("TP_EA", "EA_TP") ~ TRUE,
  rlvnt_Same_Genus$cal_two %in% c("TP_ET", "ET_TP") & rlvnt_Same_Genus$gav_two %in% c("TP_ET", "ET_TP") ~ TRUE,
  
  TRUE ~ FALSE
)

same <- subset(rlvnt_Same_Genus, Same == TRUE)

### Citrobacter freundii and Pantoea septica ####
CF_PS <- subset(same, cal_two %in% c("CF_PS", "PS_CF"), select = Gene:Same_rel_patt)
CF_PS <- mutate(CF_PS, 
                Product = c("Trigger factor", "30S ribosomal protein", "Acetyl-CoA carboxylase biotin carboxyl carrier protein", 
                            "Polyribonucleotide nucleotidyltransferase", "30S ribosomal protein", "30S ribosomal protein", "30S ribosomal protein", 
                            "Beta-ketoacyl-ACP synthase I", "DoxX family protein", "50S ribosomal protein"),
                Gene_check = c("37551_tig", "37891_rpsF", "38417_accB", "38499_pnp", "38568_rpsE", "38576_rpsQ", "38812_rpsA", "39508_fabB", "40023_yphA", 
                               "40394_rplL"))

unique(CF_PS$Gene == CF_PS$Gene_check) # If TRUE, then continue
CF_PS <- subset(CF_PS, select = Gene:Product)

write.csv(x = CF_PS, file = "8Results/CF_PS.csv", row.names = FALSE)

BM$CF_PS <- case_when(
  BM$Gene %in% CF_PS$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_CF_PS <- subset(BM, CF_PS == TRUE, select = Gene:ModelCode)

write.csv(x = BM_CF_PS, file = "8Results/Best_Model_CF_PS.csv", row.names = FALSE)

### Enterobacter cloacae and Erwinia amylovora ####
EC_EA <- subset(same, cal_two %in% c("EC_EA", "EA_EC"), select = Gene:Same_rel_patt)

EC_EA <- mutate(EC_EA, 
                Product = c("50S ribosomal protein", "ATP-dependent dethiobiotin synthetase"),
                Gene_check = c("38574_rplX", "39664_bioD1"))

unique(EC_EA$Gene == EC_EA$Gene_check) # If TRUE, then continue
EC_EA <- subset(EC_EA, select = Gene:Product)

write.csv(x = EC_EA, file = "8Results/EC_EA.csv", row.names = FALSE)

BM$EC_EA <- case_when(
  BM$Gene %in% EC_EA$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_EC_EA <- subset(BM, EC_EA == TRUE, select = Gene:ModelCode)

write.csv(x = BM_EC_EA, file = "8Results/Best_Model_EC_EA.csv", row.names = FALSE)

### Enterobacter cloacae and Erwinia tasmaniensis ####
EC_ET <- subset(same, cal_two %in% c("EC_ET", "ET_EC"), select = Gene:Same_rel_patt)

EC_ET <- mutate(EC_ET, 
                Product = c("Ribosome-associated protein", "GNAT family N-acetyltransferase"),
                Gene_check = c("37505_ybcJ", "37628_aaaT"))

unique(EC_ET$Gene == EC_ET$Gene_check) # If TRUE, then continue
EC_ET <- subset(EC_ET, select = Gene:Product)

write.csv(x = EC_ET, file = "8Results/EC_ET.csv", row.names = FALSE)

BM$EC_ET <- case_when(
  BM$Gene %in% EC_ET$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_EC_ET <- subset(BM, EC_ET == TRUE, select = Gene:ModelCode)

write.csv(x = BM_EC_EA, file = "8Results/Best_Model_EC_ET.csv", row.names = FALSE)

### Enterobacter cloacae and Pantoea agglomerans ####
EC_PA <- subset(same, cal_two %in% c("PA_EC", "EC_PA"), select = Gene:Same_rel_patt)

EC_PA <- mutate(EC_PA, 
                Product = c("Unknown*", "6-carboxytetrahydropterin synthase", "50S ribosomal protein"),
                Gene_check = c("38045_nudK", "38254_queD", "38573_rplE"))

unique(EC_PA$Gene == EC_PA$Gene_check) # If TRUE, then continue
EC_PA <- subset(EC_PA, select = Gene:Product)

write.csv(x = EC_PA, file = "8Results/EC_PA.csv", row.names = FALSE)

BM$EC_PA <- case_when(
  BM$Gene %in% EC_PA$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_EC_PA <- subset(BM, EC_PA == TRUE, select = Gene:ModelCode)

write.csv(x = BM_EC_PA, file = "8Results/Best_Model_EC_PA.csv", row.names = FALSE)

### Enterobacter cloacae and Pantoea septica ####
EC_PS <- subset(same, cal_two %in% c("EC_PS", "PS_EC"), select = Gene:Same_rel_patt)

EC_PS <- mutate(EC_PS, 
                Product = c("Recombination protein", "Spermidine/putrescine ABC transporter permease", 
                            "Spermidine/putrescine ABC transporter substrate-binding protein", "VOC family protein", 
                            "Class II fructose-bisphosphate aldolase", "Septum formation inhibitor", "30S ribosomal protein", "50S ribosomal protein", 
                            "30S ribosomal protein", "23S rRNA pseudouridine(955/2504/2580) synthase", "ATP-dependent RNA helicase", 
                            "50S ribosomal protein", "Phenylalanine--tRNA ligase subunit", "Ribonuclease T", "Pirin family protein", 
                            "Two-component system response regulator", "Peptidoglycan-binding protein", "50S ribosomal protein", 
                            "Carboxy-S-adenosyl-L-methionine synthase", "Transcription termination/antitermination protein"),
                Gene_check = c("37521_recR", "37719_ydcV", "37720_potD", "37952_cetB", "38145_fbaA", "38426_yhdE", "38579_rpsC", "38580_rplV", 
                               "38586_rpsJ", "38617_rluC", "38678_rhlB", "38895_rpmI", "38897_pheS", "38917_rnt", "39300_yhaK", "39373_arcA", "39771_ygaU", 
                               "39814_rpmB", "40266_cmoA", "40390_nusG"))






unique(EC_PS$Gene == EC_PS$Gene_check) # If TRUE, then continue
EC_PS <- subset(EC_PS, select = Gene:Product)

write.csv(x = EC_PS, file = "8Results/EC_PS.csv", row.names = FALSE)

BM$EC_PS <- case_when(
  BM$Gene %in% EC_PS$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_EC_PS <- subset(BM, EC_PS == TRUE, select = Gene:ModelCode)

write.csv(x = BM_EC_PS, file = "8Results/Best_Model_EC_PS.csv", row.names = FALSE)

### Erwinia amylovora and Pantoea agglomerans ####
EA_PA <- subset(same, cal_two %in% c("EA_PA", "PA_EA"), select = Gene:Same_rel_patt)
EA_PA <- mutate(EA_PA, 
                Product = c("DUF3561 family protein", "Succinate dehydrogenase membrane anchor subunit", "2Fe-2S ferredoxin-like protein"),
                Gene_check = c("38262_ygbE", "39149_sdhD", "39555_hypothetical_protein"))

unique(EA_PA$Gene == EA_PA$Gene_check) # If TRUE, then continue
EA_PA <- subset(EA_PA, select = Gene:Product)

write.csv(x = EA_PA, file = "8Results/EA_PA.csv", row.names = FALSE)

BM$EA_PA <- case_when(
  BM$Gene %in% EA_PA$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_EA_PA <- subset(BM, EA_PA == TRUE, select = Gene:ModelCode)

write.csv(x = BM_EA_PA, file = "8Results/Best_Model_EA_PA.csv", row.names = FALSE)

### Erwinia amylovora and Pantoea septica ####
EA_PS <- subset(same, cal_two %in% c("EA_PS", "PS_EA"), select = Gene:Same_rel_patt)
EA_PS <- mutate(EA_PS, 
                Product = c("PTS glucose transporter subunit", "Nitrogen regulation protein", "Glucose-1-phosphatase", 
                            "Phosphate ABC transporter ATP-binding protein", "Conjugal transfer protein", "4-hydroxybenzoate octaprenyltransferase",
                            "Quinone oxidoreductase", "Adenylosuccinate synthase", "Bifunctional glycosyl transferase/transpeptidase", 
                            "D-glycero-beta-D-manno-heptose 1,7-bisphosphate 7-phosphatase", "Dihydroxy-acid dehydratase", "F0F1 ATP synthase subunit",
                            "Arginine exporter", "Phosphoglycerate dehydrogenase", "2-octaprenyl-6-methoxyphenyl hydroxylase", 
                            "Amino-acid N-acetyltransferase", "DNA polymerase III subunit", "NAD(+) kinase", "Nucleotide exchange factor", 
                            "Acetyl-CoA carboxylase biotin carboxylase", "Metalloprotease", "Transcriptional regulator", 
                            "Transcription-repair coupling factor", "3'-5' ssDNA/RNA exonuclease","Protoheme IX biogenesis protein", 
                            "16S rRNA (guanine(1516)-N(2))-methyltransferase", "Transcriptional repressor", "YdiU family protein", 
                            "D-lactate dehydrogenase", "Manganese-binding transcriptional regulator", "DNA-binding protein", 
                            "Molybdopterin synthase sulfur carrier subunit", "GTP 3',8-cyclase", "Divisome-associated lipoprotein", 
                            "2,3-diphosphoglycerate-dependent phosphoglycerate mutase", "Molecular chaperone", 
                            "Division/cell wall cluster transcriptional repressor", "16S rRNA (cytosine(1402)-N(4))-methyltransferase", 
                            "Dephospho-CoA kinase", "TSUP family transporter", "UbiX family flavin prenyltransferase", "DUF412 domain-containing protein", 
                            "Bifunctional 2-polyprenyl-6-hydroxyphenol methylase/3-demethylubiquinol 3-O-methyltransferase", "Type I DNA topoisomerase", 
                            "Oligopeptide ABC transporter permease", "YicC family protein", "Rhodanese-like domain-containing protein", 
                            "YheU family protein", "Glutathione-regulated potassium-efflux system ancillary protein", "Sulfurtransferase complex subunit", 
                            "Transcriptional regulator", "30S ribosomal protein", "GrxA family glutaredoxin", "Two-component system response regulator", 
                            "Bifunctional phosphoribosyl-AMP cyclohydrolase/phosphoribosyl-ATP diphosphatase", "Hypothetical protein", 
                            "DUF1315 family protein", "Glucans biosynthesis protein", "DNA-directed RNA polymerase subunit"),
                Gene_check = c("37425_crr", "37471_glnG", "37475_yihX", "37489_pstB", "37515_hypothetical_protein", "37788_ubiA", "37796_qorA", 
                               "37885_purA", "37924_mrcB", "37964_gmhB","37986_ilvD", "38007_atpF", "38147_argO", "38151_serA", "38156_ubiH", "38231_argA", 
                               "38286_dnaQ", "38361_nadK", "38362_grpE", "38416_accC", "38429_tldD", "38435_argR", "38594_mfd", "38630_tatD", "38666_hemY",
                               "38735_rsmJ", "38880_mprA", "38901_hypothetical_protein", "38954_dld", "39082_mntR", "39091_ybiB", "39109_moaD", "39112_moaA",
                               "39299_osmY_2", "39370_cobC", "39382_dnaK", "39417_mraZ", "39418_rsmH", "39437_coaE", "39504_yfcA", "39524_ubiX", 
                               "39532_hypothetical_protein", "39558_ubiG", "39717_topA_2", "39747_oppB", "39805_hypothetical_protein", "39832_yibN", 
                               "39990_hypothetical_protein", "39994_ywrO", "40002_tusC", "40014_rcsB", "40035_rpsP", "40047_grxA", "40128_baeR", 
                               "40156_hisI", "40294_hypothetical_protein", "40327_hypothetical_protein", "40377_mdoG", "40395_rpoB"))

unique(EA_PS$Gene == EA_PS$Gene_check) # If TRUE, then continue
EA_PS <- subset(EA_PS, select = Gene:Product)

write.csv(x = EA_PS, file = "8Results/EA_PS.csv", row.names = FALSE)

BM$EA_PS <- case_when(
  BM$Gene %in% EA_PS$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_EA_PS <- subset(BM, EA_PS == TRUE, select = Gene:ModelCode)

write.csv(x = BM_EA_PS, file = "8Results/Best_Model_EA_PS.csv", row.names = FALSE)

### Erwinia amylovora and Tatumella ptyseos ####
EA_TP <- subset(same, cal_two %in% c("TP_EA", "EA_TP"), select = Gene:Same_rel_patt)
EA_TP <- mutate(EA_TP, 
                Product = c("50S ribosomal protein"),
                Gene_check = c("38470_rpmA"))

unique(EA_TP$Gene == EA_TP$Gene_check) # If TRUE, then continue
EA_TP <- subset(EA_TP, select = Gene:Product)

write.csv(x = EA_TP, file = "8Results/EA_TP.csv", row.names = FALSE)

BM$EA_TP <- case_when(
  BM$Gene %in% EA_TP$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_EA_TP <- subset(BM, EA_TP == TRUE, select = Gene:ModelCode)

write.csv(x = BM_EA_TP, file = "8Results/Best_Model_EA_TP.csv", row.names = FALSE)

### Erwinia tasmaniensis and Pantoea amylovora ####
ET_PA <- subset(same, cal_two %in% c("ET_PA", "PA_ET"), select = Gene:Same_rel_patt)
ET_PA <- mutate(ET_PA, 
                Product = c("Cytochrome o ubiquinol oxidase subunit", "YggL family protein", "NADH-quinone oxidoreductase subunit"),
                Gene_check = c("37558_cyoD", "38125_hypothetical_protein", "39545_nuoK"))

unique(ET_PA$Gene == ET_PA$Gene_check) # If TRUE, then continue
ET_PA <- subset(ET_PA, select = Gene:Product)

write.csv(x = ET_PA, file = "8Results/ET_PA.csv", row.names = FALSE)

BM$ET_PA <- case_when(
  BM$Gene %in% ET_PA$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_ET_PA <- subset(BM, ET_PA == TRUE, select = Gene:ModelCode)

write.csv(x = BM_ET_PA, file = "8Results/Best_Model_ET_PA.csv", row.names = FALSE)

#### Erwinia tasmaniensis and Pantoea septica ####
ET_PS <- subset(same, cal_two %in% c("ET_PS", "PS_ET"), select = Gene:Same_rel_patt)
ET_PS <- mutate(ET_PS, 
                Product = c("Inositol-1-monophosphatase", "YfgM family protein", 
                            "YpfN family protein", "Sulfate/thiosulfate ABC transporter permease", 
                            "Glutamate--tRNA ligase", 
                            
                            "YihA family ribosome biogenesis GTP-binding protein", "AsmA family protein", 
                            "Phosphate ABC transporter permease", "5-(carboxyamino)imidazole ribonucleotide mutase", 
                            "Efflux RND transporter periplasmic adaptor subunit", 
                            
                            "YajQ family cyclic di-GMP-binding protein", "Small ribosomal subunit biogenesis GTPase", 
                            "30S ribosomal protein", "Ribosome-associated protein", 
                            "1-deoxy-D-xylulose-5-phosphate reductoisomerase", 
                            
                            "Acetyl-CoA carboxylase carboxyl transferase subunit", 
                            "YifB family Mg chelatase-like AAA ATPase", 
                            "3',5'-cyclic-AMP phosphodiesterase", 
                            "Deoxyribose-phosphate aldolase", 
                            "Purine-nucleoside phosphorylase", 
                            
                            
                            ),
                Gene_check = c("37332_suhB", "37351_hypothetical_protein", "37398_hypothetical_protein", "37411_cysT", "37438_gltX", 
                               
                               "37467_engB", "37479_hypothetical_protein", "37488_pstA", "37509_purE", "37528_acrA", 
                               
                               "37561_hypothetical_protein", "37872_rsgA", "37893_rpsR", "37906_hypothetical_protein", "37939_dxr", 
                               
                               "37951_accA", 
                               "37990_comM", 
                               "38064_cpdA_1", 
                               "38090_deoC", 
                               "38093_deoD", 
                               ))








unique(ET_PS$Gene == ET_PS$Gene_check) # If TRUE, then continue
ET_PS <- subset(ET_PS, select = Gene:Product)

write.csv(x = ET_PS, file = "8Results/ET_PS.csv", row.names = FALSE)

BM$ET_PS <- case_when(
  BM$Gene %in% ET_PS$Gene ~ TRUE,
  TRUE ~ FALSE
)

BM_ET_PS <- subset(BM, ET_PS == TRUE, select = Gene:ModelCode)

write.csv(x = BM_ET_PS, file = "8Results/Best_Model_ET_PS.csv", row.names = FALSE)

### Erwinia tasmaniensis and Tatumella ptyseos ####
ET_TP <- subset(same, cal_two %in% c("TP_ET", "ET_TP"), select = Gene:Same_rel_patt)

### Pantoea septica and Tatumella ptyseos ####
PS_TP <- subset(same, cal_two %in% c("PS_TP", "TP_PS"), select = Gene:Same_rel_patt)

### Pantoea septica and Tatumella saanichensis ####
PS_TS <- subset(same, cal_two %in% c("PS_TS", "TS_PS"), select = Gene:Same_rel_patt)








#################


nuc_cal <- read.csv(file = "8Results/M_calida_Nucleotide.csv", stringsAsFactors = FALSE)
test <- ET_PS$Gene
length(test)
nuc_cal[which(nuc_cal$gene == test[1]), ]


check <- read.csv(file = "8Results/EA_PS.csv", stringsAsFactors = FALSE)
paste(check$Gene, ".fasta", sep = "")

#
### Kittens ###############################################################################################################################################
showmekittens()