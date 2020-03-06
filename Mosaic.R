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

M_calida <- as.data.frame(matrix(ncol = 11, nrow = 0))                    # Dataframe for M. calida's closest relatives
colnames(M_calida) <- c("Gene", "Itself_check", "First_rel", "Second_rel", "Third_rel", "Fourth_rel", "Fifth_rel", "Sixth_rel", "Seventh_rel",
                        "Eighth_rel", "Ninth_rel")                        # Changes the column names

M_gaviniae <- as.data.frame(matrix(ncol = 11, nrow = 0))                  # Dataframe for M. gaviniae's closest relatives
colnames(M_gaviniae) <- c("Gene", "Itself_check", "First_rel", 
                          "Second_rel", "Third_rel", "Fourth_rel", 
                          "Fifth_rel", "Sixth_rel", "Seventh_rel",
                          "Eighth_rel", "Ninth_rel")                      # Changes the column names

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
  
  M_cal <- as.numeric(rbind(t(dist[6, 1:6]), dist[7, 6], dist[8, 6],      # Grabs the distances for each species relative to M. calida
                            dist[9, 6], dist[10, 6]))
  M_gav <- as.numeric(rbind(t(dist[7, 1:7]), dist[8, 7], dist[9, 7],      # Grabs the distances for each species relative to M. gaviniae
                            dist[10, 7]))
  
  MCclorel <- close_relative(gene, M_cal)                                 # Calls the close_relative function
  MGclorel <- close_relative(gene, M_gav)
  
  M_calida <- rbind(M_calida, MCclorel)                                   # Combines everything together
  M_gaviniae <- rbind(M_gaviniae, MGclorel)
}
beep(8)
rm(dist, MCclorel, mega, mega2, MGclorel, gene, M_cal, M_gav, path)

write.csv(x = M_calida, file = "8Results/M_calida_Relatives.csv", row.names = FALSE)
write.csv(x = M_gaviniae, file = "8Results/M_gaviniae_Relatives.csv", row.names = FALSE)

### Categorizing closest relative #########################################################################################################################
# Weird ones for M_calida: 
#   37893_rpsR ---> 30S ribosomal subunit protein
#   38470_rpmA ---> 50S ribosomal subunit protein
#   38576_rpsQ ---> 30S ribosomal subunit protein
#   38613_rpmF ---> 50S ribosomal subunit protein
# These genes had M. gaviniae as "Itself". Therefore, M. calida and M. gaviniae have identical genes for these four.

# Weird ones for M_gaviniae:
#   38559_rplQ ---> 50S ribosomal subunit protein
#   39814_rpmB ---> 50S ribosomal subunit protein
# Same thing as above.
## These 6 genes are identical in both M. calida and M. gaviniae.

M_calida <- read.csv(file = "8Results/M_calida_Relatives.csv",            # Reads in the M. calida results
                     stringsAsFactors = FALSE)

for(row in 1:nrow(M_calida)) {                                            # To check if the first closest relative is the other Mixta species
  M_calida$Results_Mixta[row] <- case_when(
    M_calida$One[row] == "Mixta_calida" ~ M_calida$Two[row],
    M_calida$One[row] == "Mixta_gaviniae" ~ M_calida$One[row]
  )
}

for(row in 1:nrow(M_calida)) {                                            # Check the next (non-Mixta) relative
  M_calida$Results_Other[row] <- case_when(
    M_calida$Two[row] %in% c("Mixta_calida", "Mixta_gaviniae") ~ 
      M_calida$Three[row],
    TRUE ~ M_calida$Two[row]
  )
}

write.csv(x = M_calida, file = "8Results/M_calida_Relatives.csv", 
          row.names = FALSE)                                              # Write the results to a csv file



M_gaviniae <- read.csv(file = "8Results/M_gaviniae_Relatives.csv",        # Reads in the M. gaviniae results
                       stringsAsFactors = FALSE)

for(row in 1:nrow(M_gaviniae)) {                                          # To check if the first closest relative is the other Mixta species
  M_gaviniae$Results_Mixta[row] <- case_when(
    M_gaviniae$One[row] == "Mixta_gaviniae" ~ M_gaviniae$Two[row],
    M_gaviniae$One[row] == "Mixta_calida" ~ M_gaviniae$One[row]
  )
}

for(row in 1:nrow(M_gaviniae)) {                                          # Check the next (non-Mixta) relative
  M_gaviniae$Results_Other[row] <- case_when(
    M_gaviniae$Two[row] %in% c("Mixta_calida", "Mixta_gaviniae") ~ 
      M_gaviniae$Three[row],
    TRUE ~ M_gaviniae$Two[row]
  )
}

write.csv(x = M_gaviniae, file = "8Results/M_gaviniae_Relatives.csv", 
          row.names = FALSE)                                              # Write the results to a csv file



length(which(M_calida$Results_Other == "Pantoea_septica"))                # 828
length(which(M_calida$Results_Other == "Pantoea_agglomerans"))            # 81
length(which(M_calida$Results_Other == "Erwinia_amylovora"))              # 49
length(which(M_calida$Results_Other == "Erwinia_tasmaniensis"))           # 53
length(which(M_calida$Results_Other == "Tatumella_ptyseos"))              # 4
length(which(M_calida$Results_Other == "Tatumella_saanichensis"))         # 2
length(which(M_calida$Results_Other == "Citrobacter_freundii"))           # 8
length(which(M_calida$Results_Other == "Enterobacter_cloacae"))           # 10

length(which(M_gaviniae$Results_Other == "Pantoea_septica"))              # 838
length(which(M_gaviniae$Results_Other == "Pantoea_agglomerans"))          # 68
length(which(M_gaviniae$Results_Other == "Erwinia_amylovora"))            # 50
length(which(M_gaviniae$Results_Other == "Erwinia_tasmaniensis"))         # 55
length(which(M_gaviniae$Results_Other == "Tatumella_ptyseos"))            # 3
length(which(M_gaviniae$Results_Other == "Tatumella_saanichensis"))       # 3
length(which(M_gaviniae$Results_Other == "Citrobacter_freundii"))         # 7
length(which(M_gaviniae$Results_Other == "Enterobacter_cloacae"))         # 11

### Distances and standard errors #########################################################################################################################
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

for(row in 1:nrow(ID_calida)) {                                           # Set to TRUE if gene name is one of the 1035 genes in the project
  ID_calida$Check <- case_when(
    ID_calida$Gene %in% M_calida_dist$Gene ~ TRUE,
    TRUE ~ FALSE
  )
}
rm(row)

for(row in 1:nrow(ID_gaviniae)) {                                         # Set to TRUE if gene name is one of the 1035 genes in the project
  ID_gaviniae$Check <- case_when(
    ID_gaviniae$Gene %in% M_gaviniae_dist$Gene ~ TRUE,
    TRUE ~ FALSE
  )
}
rm(row)

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
M_calida_dist <- read.csv(file = "8Results/M_calida_Distances.csv", 
                          stringsAsFactors = FALSE)                       # Read in M. calida's distances
M_gaviniae_dist <- read.csv(file = "8Results/M_gaviniae_Distances.csv", 
                            stringsAsFactors = FALSE)                     # Read in M. gaviniae's distances

tidy_dist_C <- subset(M_calida_dist, select = c("ID", "Gene", "Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", 
                                                "Erwinia_amylovora", "Erwinia_tasmaniensis", "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", 
                                                "Pantoea_septica", 
                                                "Tatumella_ptyseos"))     # Subset for only the distances
tidy_dist_C <- tidy_dist_C %>%                                            # Gather distances into one column according to species and by ID
  pivot_longer(cols = Tatumella_saanichensis:Tatumella_ptyseos,
               names_to = "Species", values_to = "Distance")

tidy_ster_C <- subset(M_calida_dist, select = c("ID", "Gene", "TS_Error", "CF_Error", "EC_Error", "EA_Error", "ET_Error", "MC_Error", "MG_Error", 
                                                "PA_Error", "PS_Error", 
                                                "TP_Error"))              # Subset for only the standard errors
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
tidy_dist_G <- subset(M_gaviniae_dist, select = c("ID", "Gene", "Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", 
                                                  "Erwinia_amylovora", "Erwinia_tasmaniensis", "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", 
                                                  "Pantoea_septica", 
                                                  "Tatumella_ptyseos"))   # Subset for only the distances
tidy_dist_G <- tidy_dist_G %>%                                            # Gather distances into one column according to species and by ID
  pivot_longer(cols = Tatumella_saanichensis:Tatumella_ptyseos,
               names_to = "Species", values_to = "Distance")

tidy_ster_G <- subset(M_gaviniae_dist, select = c("ID", "Gene", "TS_Error", "CF_Error", "EC_Error", "EA_Error", "ET_Error", "MC_Error", "MG_Error", 
                                                  "PA_Error", "PS_Error", 
                                                  "TP_Error"))            # Subset for only the standard errors
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
### Mixta calida ###
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

### Mixta gaviniae ###
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

### M. calida plot ###
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

### M. gaviniae plot ###
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
### Min, Mean, and Max ###
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
M_calida_sdist <- read.csv(file = "8Results/M_calida_Sort_Dist.csv",
                           stringsAsFactors = FALSE)    
M_calida_sdist$Min_ConfInt <- M_calida_sdist$Distance - 
  (M_calida_sdist$Std_Errors * 1.96)                                      # Lower confidence interval
M_calida_sdist$Max_ConfInt <- M_calida_sdist$Distance + 
  (M_calida_sdist$Std_Errors * 1.96)                                      # Upper confidence interval

M_calida_sdist <- subset(M_calida_sdist, 
                         Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                        "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                        "Tatumella_ptyseos"))             # Remove Mixta species

M_gaviniae_sdist <- read.csv(file = "8Results/M_gaviniae_Sort_Dist.csv",
                             stringsAsFactors = FALSE)
M_gaviniae_sdist$Min_ConfInt <- M_gaviniae_sdist$Distance - 
  (M_gaviniae_sdist$Std_Errors * 1.96)                                    # Lower confidence interval
M_gaviniae_sdist$Max_ConfInt <- M_gaviniae_sdist$Distance + 
  (M_gaviniae_sdist$Std_Errors * 1.96)                                    # Upper confidence interval

M_gaviniae_sdist <- subset(M_gaviniae_sdist, 
                           Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                          "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                          "Tatumella_ptyseos"))           # Remove Mixta species

uniq_gene <- as.data.frame(unique(M_calida_sdist$Gene))                   # Dataframe with the 1035 gene names
colnames(uniq_gene) <- "Gene_name"

### Significant genes by species ###
# M. calida #
sig_genes <- as.data.frame(matrix(nrow = 0, ncol = 0))
for(row in 1:nrow(uniq_gene)) {                                           # First relative's distance is significantly different from next species
  gene <- subset(M_calida_sdist, Gene == uniq_gene$Gene_name[row])
  
  x <- Rfast::nth(x = gene$Distance, k = 1, descending = FALSE,           # First relative
                  index.return = TRUE)
  y <- Rfast::nth(x = gene$Distance, k = 2, descending = FALSE,           # Second relative
                  index.return = TRUE)
  
  sig <- as.data.frame(matrix(nrow = 1, ncol = 0))
  sig <- mutate(sig,                                                      # Creates a dataframe with gene names and if significantly different
                Gene_name = uniq_gene$Gene_name[row],
                Check = case_when(
                  gene$Max_ConfInt[x] < gene$Min_ConfInt[y] ~ TRUE,       # Distances are significantly different
                  TRUE ~ FALSE                                            # Distances are not significantly different
                ))
  
  sig_genes <- rbind(sig_genes, sig)
}
rm(gene, sig, row, x, y)

Mcalida_spp_sig_genes <- subset(sig_genes, Check == TRUE,
                                select = Gene_name)                       # Subset for significantly different distances, keep only gene names

write.csv(x = Mcalida_spp_sig_genes, file = "8Results/M_calida_Species_Sig.csv", 
          row.names = FALSE)

# M.gaviniae #
sig_genes <- as.data.frame(matrix(nrow = 0, ncol = 0))
for(row in 1:nrow(uniq_gene)) {                                           # First relative's distance is significantly different from next species
  gene <- subset(M_gaviniae_sdist, Gene == uniq_gene$Gene_name[row])
  
  x <- Rfast::nth(x = gene$Distance, k = 1, descending = FALSE,           # First relative
                  index.return = TRUE)
  y <- Rfast::nth(x = gene$Distance, k = 2, descending = FALSE,           # Second relative
                  index.return = TRUE)
  
  sig <- as.data.frame(matrix(nrow = 1, ncol = 0))
  sig <- mutate(sig,                                                      # Creates a dataframe with gene names and if significantly different
                Gene_name = uniq_gene$Gene_name[row],
                Check = case_when(
                  gene$Max_ConfInt[x] < gene$Min_ConfInt[y] ~ TRUE,       # Distances are significantly different
                  TRUE ~ FALSE                                            # Distances are not significantly different
                ))
  
  sig_genes <- rbind(sig_genes, sig)
}
rm(gene, sig, row, x, y)

Mgaviniae_spp_sig_genes <- subset(sig_genes, Check == TRUE, 
                                  select = Gene_name)                     # Subset for significantly different distances, keep only gene names

write.csv(x = Mgaviniae_spp_sig_genes, file = "8Results/M_gaviniae_Species_Sig.csv", 
          row.names = FALSE)

### Significant genes by genus ###
# M. calida #
sig_genes <- as.data.frame(matrix(nrow = 0, ncol = 0))
for(row in 1:nrow(uniq_gene)) {                                           # First relative's distance is significantly different from next genus
  gene <- subset(M_calida_sdist, Gene == uniq_gene$Gene_name[row])
  
  x <- Rfast::nth(x = gene$Distance, k = 1, descending = FALSE,           # First relative
                  index.return = TRUE)
  y <- Rfast::nth(x = gene$Distance, k = 2, descending = FALSE,           # Second relative
                  index.return = TRUE)
  z <- Rfast::nth(x = gene$Distance, k = 3, descending = FALSE,           # Third relative
                  index.return = TRUE)
  
  genus_check <- case_when(
    x %in% c(1, 8) & y %in% c(1, 8) ~ TRUE,                               # Tatumella
    x %in% c(6, 7) & y %in% c(6, 7) ~ TRUE,                               # Pantoea
    x %in% c(4, 5) & y %in% c(4, 5) ~ TRUE,                               # Erwinia
    TRUE ~ FALSE
  )
  
  sig <- as.data.frame(matrix(nrow = 1, ncol = 0))
  sig <- mutate(sig,                                                      # Creates a dataframe with gene names, if first two are same genus, and if 
                Gene_name = uniq_gene$Gene_name[row],                     # significantly different
                Same_Genus = genus_check,
                Check = case_when(
                  genus_check == TRUE ~ case_when(                        # If first two are same genus, compare first and third relatives
                    gene$Max_ConfInt[x] < gene$Min_ConfInt[z] ~ TRUE,     # Distances are significantly different
                    TRUE ~ FALSE                                          # Distances are not significantly different
                  ),
                  genus_check == FALSE ~ case_when(                       # If first two are not in the same genus, compare first and second relatives
                    gene$Max_ConfInt[x] < gene$Min_ConfInt[y] ~ TRUE,     # Distances are significantly different
                    TRUE ~ FALSE                                          # Distances are not significantly different
                  )
                ))
  
  sig_genes <- rbind(sig_genes, sig)
}
rm(gene, sig, row, genus_check, x, y, z)

M_calida_gen_sig <- subset(sig_genes, Check == TRUE, select = Gene_name)  # Subset for significantly different distances, keep only gene names
write.csv(x = M_calida_gen_sig, file = "8Results/M_calida_Genus_Sig.csv", 
          row.names = FALSE)

M_calida_same_gen <- subset(sig_genes, select = Gene_name:Same_Genus)     # Keep gene names and whether the first two relatives are from the same genus
write.csv(x = M_calida_same_gen, file = "8Results/M_calida_Same_Genus.csv", 
          row.names = FALSE)

# M. gaviniae #
sig_genes <- as.data.frame(matrix(nrow = 0, ncol = 0))
for(row in 1:nrow(uniq_gene)) {                                           # First relative's distance is significantly different from next genus
  gene <- subset(M_gaviniae_sdist, Gene == uniq_gene$Gene_name[row])
  
  x <- Rfast::nth(x = gene$Distance, k = 1, descending = FALSE,           # First relative
                  index.return = TRUE)
  y <- Rfast::nth(x = gene$Distance, k = 2, descending = FALSE,           # Second relative
                  index.return = TRUE)
  z <- Rfast::nth(x = gene$Distance, k = 3, descending = FALSE,           # Third relative
                  index.return = TRUE)
  
  genus_check <- case_when(
    x %in% c(1, 8) & y %in% c(1, 8) ~ TRUE,                               # Tatumella
    x %in% c(6, 7) & y %in% c(6, 7) ~ TRUE,                               # Pantoea
    x %in% c(4, 5) & y %in% c(4, 5) ~ TRUE,                               # Erwinia
    TRUE ~ FALSE
  )
  
  sig <- as.data.frame(matrix(nrow = 1, ncol = 0))
  sig <- mutate(sig,                                                      # Creates a dataframe with gene names, if first two are same genus, and if 
                Gene_name = uniq_gene$Gene_name[row],                     # significantly different
                Same_Genus = genus_check,
                Check = case_when(
                  genus_check == TRUE ~ case_when(                        # If first two are same genus, compare first and third relatives
                    gene$Max_ConfInt[x] < gene$Min_ConfInt[z] ~ TRUE,     # Distances are significantly different
                    TRUE ~ FALSE                                          # Distances are not significantly different
                  ),
                  genus_check == FALSE ~ case_when(                       # If first two are not in the same genus, compare first and second relatives
                    gene$Max_ConfInt[x] < gene$Min_ConfInt[y] ~ TRUE,     # Distances are significantly different
                    TRUE ~ FALSE                                          # Distances are not significantly different
                  )
                ))
  
  sig_genes <- rbind(sig_genes, sig)
}
rm(gene, sig, row, genus_check, x, y, z)

M_gaviniae_gen_sig <- subset(sig_genes, Check == TRUE, 
                             select = Gene_name)                          # Subset for significantly different distances, keep only gene names
write.csv(x = M_gaviniae_gen_sig, file = "8Results/M_gaviniae_Genus_Sig.csv", row.names = FALSE)

M_gaviniae_same_gen <- subset(sig_genes, select = Gene_name:Same_Genus)   # Keep gene names and whether the first two relatives are from the same genus
write.csv(x = M_gaviniae_same_gen, file = "8Results/M_gaviniae_Same_Genus.csv", row.names = FALSE)

### Plots ###
M_calida_S_sig <- read.csv(file = "8Results/M_calida_Species_Sig.csv", 
                           stringsAsFactors = FALSE)
M_calida_G_sig <- read.csv(file = "8Results/M_calida_Genus_Sig.csv", 
                           stringsAsFactors = FALSE)
M_calida_Same_Genus <- read.csv(file = "8Results/M_calida_Same_Genus.csv", 
                                stringsAsFactors = FALSE)
M_calida_Same_Genus <- subset(M_calida_Same_Genus, Same_Genus == TRUE)


M_gaviniae_S_sig <- read.csv(file = "8Results/M_gaviniae_Species_Sig.csv", 
                             stringsAsFactors = FALSE)
M_gaviniae_G_sig <- read.csv(file = "8Results/M_gaviniae_Genus_Sig.csv", 
                             stringsAsFactors = FALSE)
M_gaviniae_Same_Genus <- read.csv(file = "8Results/M_gaviniae_Same_Genus.csv", 
                                  stringsAsFactors = FALSE)
M_gaviniae_Same_Genus <- subset(M_gaviniae_Same_Genus, Same_Genus == TRUE)


M_cal_dist <- read.csv(file = "8Results/M_calida_Sort_Dist.csv", stringsAsFactors = FALSE)

M_cal_dist <- subset(M_cal_dist, Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                                "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                                "Tatumella_ptyseos"))     # Removes the Mixta species since they are most likely ~ 0

M_cal_dist <- mutate(M_cal_dist,
                     Sig_Species = case_when(
                       M_cal_dist$Gene %in% M_calida_S_sig$Gene_name ~ TRUE,
                       TRUE ~ FALSE
                     ),
                     Sig_Genus = case_when(
                       M_cal_dist$Gene %in% M_calida_G_sig$Gene_name ~ TRUE,
                       TRUE ~ FALSE
                     ),
                     spp_gen_check = case_when(
                       Sig_Species == TRUE & Sig_Genus == TRUE ~ 1,
                       Sig_Species == TRUE & Sig_Genus == FALSE ~ 2,
                       Sig_Species == FALSE & Sig_Genus == TRUE ~ 3,
                       TRUE ~ 4
                     ))

length(which(M_cal_dist$spp_gen_check == 0)) == 0 # If TRUE, then continue

M_cal_dist <- subset(M_cal_dist, spp_gen_check %in% c(1, 3), select = ID:Sig_Genus)

M_cal_dist <- mutate(M_cal_dist,
                     DistanceN = Distance * -1,
                     Both = case_when(
                       M_cal_dist$Sig_Species == M_cal_dist$Sig_Genus ~ "Both",
                       TRUE ~ "Genus Only"
                     ),
                     upper = DistanceN + (Std_Errors * 1.96),
                     lower = DistanceN - (Std_Errors * 1.96))

# M_cal_dist <- subset(M_cal_dist, Distance <= 1)

Mcal_sig <- ggplot(data = M_cal_dist, aes(x = ID, y = DistanceN)) +
  # geom_point(aes(colour = Species), show.legend = TRUE) +
  # geom_errorbar(aes(ymin = lower, ymax = upper, colour = Species), width = 0.2, position = position_dodge(0.05)) +
  # scale_colour_manual("Species",
  #                     values = alpha(c("red", "orange", "darkgreen", "green3", "blue3", "dodgerblue2", "darkorchid", "violetred1")),
  #                     labels = c("Citrobacter freundii", "Enterobacter cloacae", "Erwinia amylovora", "Erwinia tasmaniensis",
  #                                "Pantoea agglomerans", "Pantoea septica", "Tatumella ptyseos", "Tatumella saanichensis")) +
  # theme(legend.position = "bottom", text = element_text(size = 9),
  #       legend.text = element_text(face = "italic")) +
  # labs(x = expression(paste(italic("M. calida"), " Gene ID")),
  #      y = expression(paste("Negative Distance from ", italic("M. calida")))) +
  # scale_x_continuous(breaks = round(seq(min(M_cal_dist$ID), max(M_cal_dist$ID), by = 1000), -2),
  #                    limits = c(0, 4032), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(-5.7, 0.1), expand = c(0, 0)) +
  geom_bar(data = test, aes(x = ID, y = rnorm, fill = Both), stat = "identity", inherit.aes = FALSE); Mcal_sig
  # scale_fill_manual("Sig. Genus and Species", values = rep(1, 2), guide = guide_legend(override.aes = list(fill = c("black", "pink"), colour = c("black", "pink"))))
#

test <- as.data.frame(matrix(ncol = 0, nrow = 153))
test <- mutate(test,
               ID = unique(M_cal_dist$ID),
               Gene = unique(M_cal_dist$Gene))

unique(M_calida_G_sig$Gene_name == test$Gene) # If TRUE, then continue

test <- mutate(test,
               Sig_Species = case_when(
                 test$Gene %in% M_calida_S_sig$Gene_name ~ TRUE,
                 TRUE ~ FALSE
               ),
               Sig_Genus = case_when(
                 test$Gene %in% M_calida_G_sig$Gene_name ~ TRUE,
                 TRUE ~ FALSE
               ))

test <- mutate(test,
               Both = case_when(
                 test$Sig_Species == test$Sig_Genus ~ "Both",
                 TRUE ~ "Genus Only"
               ))

test$rnorm <- rnorm(n = 153, mean = 5, sd = 1)


# ggsave(Mcal_sig, file = "9_1Plots_calida/MC_dist_sigSG.png", 
#        width = 16.51, height = 12.38, units = "cm")







M_gav_dist <- read.csv(file = "8Results/M_gaviniae_Sort_Dist.csv", stringsAsFactors = FALSE)

M_gav_dist <- subset(M_gav_dist, Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                                "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                                "Tatumella_ptyseos"))     # Removes the Mixta species since they are most likely ~ 0

M_gav_dist <- mutate(M_gav_dist,
                     Sig_Species = case_when(
                       M_gav_dist$Gene %in% M_gaviniae_S_sig$Gene_name ~ TRUE,
                       TRUE ~ FALSE
                     ),
                     Sig_Genus = case_when(
                       M_gav_dist$Gene %in% M_gaviniae_G_sig$Gene_name ~ TRUE,
                       TRUE ~ FALSE
                     ),
                     spp_gen_check = case_when(
                       Sig_Species == TRUE & Sig_Genus == TRUE ~ 1,
                       Sig_Species == TRUE & Sig_Genus == FALSE ~ 2,
                       Sig_Species == FALSE & Sig_Genus == TRUE ~ 3,
                       TRUE ~ 4
                     ))

length(which(M_gav_dist$spp_gen_check == 0)) == 0 # If TRUE, then continue

M_gav_dist <- subset(M_gav_dist, spp_gen_check %in% c(1, 3), select = ID:Sig_Genus)

M_gav_dist <- mutate(M_gav_dist,
                     DistanceN = Distance * -1,
                     Both = case_when(
                       M_gav_dist$Sig_Species == M_gav_dist$Sig_Genus ~ "Both",
                       TRUE ~ "Genus Only"
                     ),
                     upper = DistanceN + (Std_Errors * 1.96),
                     lower = DistanceN - (Std_Errors * 1.96))
# dtst
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
### Kittens ###############################################################################################################################################
showmekittens()


# c("Tatumella_saanichensis", "Citrobacter_freundii",	"Enterobacter_cloacae", "Erwinia_amylovora", "Erwinia_tasmaniensis", 
#                               "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", "Pantoea_septica", "Tatumella_ptyseos")
