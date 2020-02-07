### Information ############################################################################################################################################
# Author: Kim Hinz
# Date of study: 2019-09-05 -- 2020-04-01
# Purpose: Phylogenetic analyses of Mixta genes to determine origination.
# Name of study: Mosaic Mixta


# The following code does different things for the phylogenetic study of two Mixta species (bacteria). Depending on the genes, model, and statistical 
# method chosen for the phylogenetic tree, the Mixta species group with different genera (see Palmer et al. 2018 and Rezzonico et al. 2016 for an example).
# Primarily, Mixta appeas to be a close relative to Pantoea with some leaning towards Erwinia. The purpose of my research is to determine why this might be
# the case by performing distance matrix analyses for the homologous genes between two Mixta species, two Pantoea species, two Erwinia species, two 
# Tatumella species, one Citrobacter species, and one Enterobacter species. The Citrobacter and Enterobacter species form the outgroup.


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

### Packages ###############################################################################################################################################
library("seqinr")

library("plyr")
library("msa")
library("beepr")


library("ape")
library("adegenet")
library("Rfast")
library("dplyr")
library("plyr")
library("tidyr")
library("ggplot2")
theme_set(theme_bw())

setwd("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/")

#
### Selection for genes ###################################################################################################################################
# This portion of the code acts as filter; it passes forward files that have full sequences (genes aren't split) and that don't have truncated sequences
# (genes must be at least 90% of the length of the longest gene in each file). For example, if the longest sequence is 1000 bp, then the rest of the 
# sequences in that file must be at least 900 bp. If at least one is shorter than 900 bp, the whole file is excluded.

fastaFiles <- as.data.frame(list.files(path = "3Homologues_10/",
                                       pattern = ".fasta"))               # Makes a dataframe where the first column is a list of fasta gene files
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("3Homologues_10/", fastaFiles$File_name, 
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
  gene_file <- read.table(file = path, header = FALSE, sep = "\n", 
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
# This section aligns the genes that passed the filter using ClustalW through the R package msa.

fastaFilesOrg <- as.data.frame(list.files(path = "4Organize/", 
                                          pattern = ".fasta"))            # Makes a dataframe listing the fasta files in the folder
colnames(fastaFilesOrg) <- "File_name"                                    # Changes the column name
fastaFilesOrg$Path_name <- paste("4Organize/", fastaFilesOrg$File_name, 
                                 sep = "")                                # Creates a file pathway for each gene

align_gene <- function(gene_file) {                                       # Aligns the sequences in the file and converts to a dataframe
  align <- msa::msaClustalW(inputSeqs = gene_file, maxiters = 100, 
                            type = "dna", order = "input")                # Aligns the genes
  
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
              open = "w", nbchar = 10000, as.string = TRUE)               # Creates a fasta file for each gene, will continue on to distance matrices
}
beep(8)
rm(alignFA, gene_file, path, row)

### Best model for genes ##################################################################################################################################
fastaFilesModel <- as.data.frame(list.files(path = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/6Model/",
                                            pattern = ".csv"))            # Make a dataframe listing the csv files in the folder
colnames(fastaFilesModel) <- "File_name"                                  # Changes the column name
fastaFilesModel$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/6Model/",
                                   fastaFilesModel$File_name, sep = "")   # Creates a file pathway for each gene

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
  
  model$ModelCode <- case_when(                                           # Get the model codes
    model$Model1 %in% c("JC", "JC+I") ~ "JC",
    model$Model1 %in% c("JC+G", "JC+G+I") ~ "JC_G",
    model$Model1 %in% c("K2", "K2+I") ~ "K2",
    model$Model1 %in% c("K2+G", "K2+G+I") ~ "K2_G",
    model$Model1 %in% c("T92", "T92+I") ~ "T92",
    model$Model1 %in% c("T92+G", "T92+G+I") ~ "T92_G",
    model$Model1 %in% c("TN93", "TN93+I") ~ "TN93",
    model$Model1 %in% c("TN93+G", "TN93+G+I") ~ "TN93_G",
  )
  best_model <- rbind(best_model, model)                                  # Combine all best models to one dataset
}
rm(gene_model_test, model, path, row)                                     # Remove unneeded variables from the for loop

best_model$Gene <- gsub(pattern = "-4212.csv", replacement = "",          # Creates a column with just the gene name
                        x = best_model$File_name)
best_model$Path_Name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/5Aligned/",
                              best_model$Gene, ".fasta", sep = "")        # New pathway to aligned fasta file


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
megFiles <- as.data.frame(list.files(path = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/7Distance/",
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
  
  rela <- data.frame(Gene = gen,
                     One = relative(min1),
                     Two = relative(min2),
                     Three = relative(min3))
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

M_calida <- as.data.frame(matrix(ncol = 4, nrow = 0))                     # Dataframe for M. calida's closest relatives
colnames(M_calida) <- c("Gene", "Itself_check", "First_relative",
                        "Second_relative")                                # Changes the column names

M_gaviniae <- as.data.frame(matrix(ncol = 4, nrow = 0))                   # Dataframe for M. gaviniae's closest relatives
colnames(M_gaviniae) <- c("Gene", "Itself_check", "First_relative", 
                          "Second_relative")                              # Changes the column names

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

M_calida <- read.csv(file = "8Results/M_calida.csv",                      # Reads in the M. calida results
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

write.csv(x = M_calida, file = "8Results/M_calida.csv")                   # Write the results to a csv file



M_gaviniae <- read.csv(file = "8Results/M_gaviniae.csv",                  # Reads in the M. gaviniae results
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

write.csv(x = M_gaviniae, file = "8Results/M_gaviniae.csv")               # Write the results to a csv file



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
megFiles <- as.data.frame(list.files(path = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/7Distance/",
                                     pattern = ".meg"))                   # Makes a list of all .meg file in this diretory
colnames(megFiles) <- "File_name"                                         # Changes the column name
megFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/7Distance/",
                            megFiles$File_name, sep = "")                 # Adds the path name for each gene
megFiles$Gene <- best_model$Gene                                          # Adds the gene name (no extension)


M_calida_dist <- as.data.frame(matrix(ncol = 21, nrow = 0))               # Dataframe for M. calida's closest relatives
colnames(M_calida_dist) <- c("Gene", "Tatumella_saanichensis", "TS_Error", "Citrobacter_freundii", "CF_Error", "Enterobacter_cloacae", "EC_Error",
                             "Erwinia_amylovora", "EA_Error", "Erwinia_tasmaniensis", "ET_Error", "Mixta_calida", "MC_Error", "Mixta_gaviniae", 
                             "MG_Error", "Pantoea_agglomerans", "PA_Error", "Pantoea_septica", "PS_Error", 
                             "Tatumella_ptyseos", "TP_Error")             # Changes the column names

M_gaviniae_dist <- as.data.frame(matrix(ncol = 21, nrow = 0))             # Dataframe for M. gaviniae's closest relatives
colnames(M_gaviniae_dist) <- c("Gene", "Tatumella_saanichensis", "TS_Error", "Citrobacter_freundii", "CF_Error", "Enterobacter_cloacae", "EC_Error",
                               "Erwinia_amylovora", "EA_Error", "Erwinia_tasmaniensis", "ET_Error", "Mixta_calida", "MC_Error", "Mixta_gaviniae", 
                               "MG_Error", "Pantoea_agglomerans", "PA_Error", "Pantoea_septica", "PS_Error", 
                               "Tatumella_ptyseos", "TP_Error")           # Changes the column names

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

write.csv(x = M_calida_dist, file = "8Results/M_calida_dist.csv", row.names = FALSE)
write.csv(x = M_gaviniae_dist, file = "8Results/M_gaviniae_dist.csv", row.names = FALSE)

### Retrieve gene IDs #####################################################################################################################################
M_calida_dist <- read.csv(file = "8Results/M_calida_dist.csv")
M_gaviniae_dist <- read.csv(file = "8Results/M_gaviniae_dist.csv")

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

write.csv(x = M_calida_dist, file = "8Results/M_calida_dist.csv",         # Save distances, standard errors, and gene IDs
          row.names = FALSE)
write.csv(x = M_gaviniae_dist, file = "8Results/M_gaviniae_dist.csv", 
          row.names = FALSE)

### Distance plot #########################################################################################################################################
M_calida_dist <- read.csv(file = "8Results/M_calida_dist.csv", 
                          stringsAsFactors = FALSE)                       # Read in M. calida's distances
M_gaviniae_dist <- read.csv(file = "8Results/M_gaviniae_dist.csv", 
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

write.csv(x = M_calida_sort, file = "8Results/M_calida_sort.csv",         # Write this to a csv file
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

write.csv(x = M_gaviniae_sort, file = "8Results/M_gaviniae_sort.csv", 
          row.names = FALSE)                                              # Write this to a csv file


### M. calida plots #######################################################################################################################################
M_calida_sort <- read.csv(file = "8Results/M_calida_sort.csv", 
                          stringsAsFactors = FALSE)                       # Read in the tidied and sorted distances for M. calida
M_cal1 <- subset(M_calida_sort, Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                               "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                               "Tatumella_ptyseos"))      # Removes the Mixta species since they are most likely ~ 0
M_cal1$DistanceN <- M_cal1$Distance * -1                                  # Creates a column with negative distances (so 0 will be at top of plot)

png("9_1Plots_calida/MC_Full_dist.png", width = 2000, height = 1200)
ggplot(data = M_cal1, aes(x = ID, y = DistanceN)) +                       # Full plot
  geom_point(aes(colour = M_cal1$Species), size = 2, alpha = 0.75) +
  geom_line(aes(colour = M_cal1$Species), linetype = "dotted") +
  theme(legend.position = "bottom", text = element_text(size = 36)) +
  labs(x = "M. calida Gene ID", y = "Distance from M. calida", colour = "Species")
dev.off()

# Separated into 8 groups
dist_plot <- function(dtst) {
  ggplot(data = dtst, aes(x = ID, y = DistanceN)) +
    geom_point(aes(colour = dtst$Species), size = 3, alpha = 0.75) +
    geom_line(aes(colour = dtst$Species), linetype = "dotted") +
    theme(legend.position = "bottom", text = element_text(size = 20)) +
    labs(x = "M. calida Gene ID", y = "Distance from M. calida", 
         colour = "Species") +
    ylim(-5.01, 0) +
    geom_errorbar(aes(ymin = DistanceN - Std_Errors, 
                      ymax = DistanceN + Std_Errors, colour = Species), 
                  width = 0.2, position = position_dodge(0.05))
}

png("9_1Plots_calida/MC_dist_1_8.png", width = 1000, height = 750)
dist_plot(subset(M_cal1, ID <= 510))
dev.off()

png("9_1Plots_calida/MC_dist_2_8.png", width = 1000, height = 750)
dist_plot(subset(M_cal1, ID %in% 511:1021))
dev.off()

png("9_1Plots_calida/MC_dist_3_8.png", width = 1000, height = 750)
dist_plot(subset(M_cal1, ID %in% 1021:1531))
dev.off()

png("9_1Plots_calida/MC_dist_4_8.png", width = 1000, height = 750)
dist_plot(subset(M_cal1, ID %in% 1532:2042))
dev.off()

png("9_1Plots_calida/MC_dist_5_8.png", width = 1000, height = 750)
dist_plot(subset(M_cal1, ID %in% 2043:2552))
dev.off()

png("9_1Plots_calida/MC_dist_6_8.png", width = 1000, height = 750)
dist_plot(subset(M_cal1, ID %in% 2553:3063))
dev.off()

png("9_1Plots_calida/MC_dist_7_8.png", width = 1000, height = 750)
dist_plot(subset(M_cal1, ID %in% 3064:3573))
dev.off()

png("9_1Plots_calida/MC_dist_8_8.png", width = 1000, height = 750)
dist_plot(subset(M_cal1, ID %in% 3574:4084))
dev.off()

# Separated into 8 groups and distances no greater than 1
M_cal3 <- subset(M_cal1, Distance < 1)

dist_plot_one <- function(dtst) {
  ggplot(data = dtst, aes(x = ID, y = DistanceN)) +
    geom_point(aes(colour = dtst$Species), size = 3, alpha = 0.75) +
    # geom_line(aes(colour = dtst$Species), linetype = "dotted") +
    theme(legend.position = "bottom", text = element_text(size = 20)) +
    labs(x = "M. calida Gene ID", y = "Distance from M. calida", 
         colour = "Species") +
    geom_errorbar(aes(ymin = DistanceN - Std_Errors, 
                      ymax = DistanceN + Std_Errors, colour = Species), 
                  width = 0.2, position = position_dodge(0.05))
}

png("9_1Plots_calida/MC_one_dist_1_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_cal3, ID <= 510))
dev.off()

png("9_1Plots_calida/MC_one_dist_2_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_cal3, ID %in% 511:1021))
dev.off()

png("9_1Plots_calida/MC_one_dist_3_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_cal3, ID %in% 1022:1531))
dev.off()

png("9_1Plots_calida/MC_one_dist_4_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_cal3, ID %in% 1532:2042))
dev.off()

png("9_1Plots_calida/MC_one_dist_5_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_cal3, ID %in% 2043:2552))
dev.off()

png("9_1Plots_calida/MC_one_dist_6_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_cal3, ID %in% 2553:3063))
dev.off()

png("9_1Plots_calida/MC_one_dist_7_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_cal3, ID %in% 3064:3573))
dev.off()

png("9_1Plots_calida/MC_one_dist_8_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_cal3, ID %in% 3574:4084))
dev.off()



png("9_1Plots_calida/MC_one_dist_1_4.png", width = 2000, height = 1000)
dist_plot_one(subset(M_gav3, ID %in% 1:1021))
dev.off()

png("9_1Plots_calida/MC_one_dist_2_4.png", width = 2000, height = 1000)
dist_plot_one(subset(M_gav3, ID %in% 1022:2042))
dev.off()

png("9_1Plots_calida/MC_one_dist_3_4.png", width = 2000, height = 1000)
dist_plot_one(subset(M_gav3, ID %in% 2043:3063))
dev.off()

png("9_1Plots_calida/MC_one_dist_4_4.png", width = 2000, height = 1000)
dist_plot_one(subset(M_gav3, ID %in% 3064:4084))
dev.off()

### M. gaviniae plots #####################################################################################################################################
M_gaviniae_sort <- read.csv(file = "8Results/M_gaviniae_sort.csv", 
                            stringsAsFactors = FALSE)                     # Read in the tidied and sorted distances for M. gaviniae

M_gav1 <- subset(M_gaviniae_sort, Species %in% c("Tatumella_saanichensis", "Citrobacter_freundii", "Enterobacter_cloacae", "Erwinia_amylovora", 
                                                 "Erwinia_tasmaniensis", "Pantoea_agglomerans", "Pantoea_septica", 
                                                 "Tatumella_ptyseos"))      # Removes the Mixta species since they are most likely ~ 0
M_gav1$DistanceN <- M_gav1$Distance * -1                                  # Creates a column with negative distances (so 0 will be at top of plot)

png("9_2Plots_gaviniae/MG_Full_dist.png", width = 2000, height = 1200)
ggplot(data = M_gav1, aes(x = ID, y = DistanceN)) +                       # Full plot
  geom_point(aes(colour = M_gav1$Species), size = 2, alpha = 0.75) +
  geom_line(aes(colour = M_gav1$Species), linetype = "dotted") +
  theme(legend.position = "bottom", text = element_text(size = 36)) +
  labs(x = "M. gaviniae Gene ID", y = "Distance from M. gaviniae", colour = "Species")
dev.off()

# Separated into 8 groups
dist_plot <- function(dtst) {
  ggplot(data = dtst, aes(x = ID, y = DistanceN)) +
    geom_point(aes(colour = dtst$Species), size = 3, alpha = 0.75) +
    geom_line(aes(colour = dtst$Species), linetype = "dotted") +
    theme(legend.position = "bottom", text = element_text(size = 20)) +
    labs(x = "M. gaviniae Gene ID", y = "Distance from M. gaviniae", 
         colour = "Species") +
    ylim(-5.25, 0) +
    geom_errorbar(aes(ymin = DistanceN - Std_Errors, 
                      ymax = DistanceN + Std_Errors, colour = Species), 
                  width = 0.2, position = position_dodge(0.05))
}

png("9_2Plots_gaviniae/MG_dist_1_8.png", width = 1000, height = 750)
dist_plot(subset(M_gav1, ID <= 531))
dev.off()

png("9_2Plots_gaviniae/MG_dist_2_8.png", width = 1000, height = 750)
dist_plot(subset(M_gav1, ID %in% 532:1063))
dev.off()

png("9_2Plots_gaviniae/MG_dist_3_8.png", width = 1000, height = 750)
dist_plot(subset(M_gav1, ID %in% 1064:1595))
dev.off()

png("9_2Plots_gaviniae/MG_dist_4_8.png", width = 1000, height = 750)
dist_plot(subset(M_gav1, ID %in% 1596:2127))
dev.off()

png("9_2Plots_gaviniae/MG_dist_5_8.png", width = 1000, height = 750)
dist_plot(subset(M_gav1, ID %in% 2128:2659))
dev.off()

png("9_2Plots_gaviniae/MG_dist_6_8.png", width = 1000, height = 750)
dist_plot(subset(M_gav1, ID %in% 2660:3191))
dev.off()

png("9_2Plots_gaviniae/MG_dist_7_8.png", width = 1000, height = 750)
dist_plot(subset(M_gav1, ID %in% 3192:3723))
dev.off()

png("9_2Plots_gaviniae/MG_dist_8_8.png", width = 1000, height = 750)
dist_plot(subset(M_gav1, ID %in% 3724:4255))
dev.off()

# Separated into 8 groups and distances no greater than 1
M_gav3 <- subset(M_gav1, Distance < 1)

dist_plot_one <- function(dtst) {
  ggplot(data = dtst, aes(x = ID, y = DistanceN)) +
    geom_point(aes(colour = dtst$Species), size = 3, alpha = 0.75) +
    # geom_line(aes(colour = dtst$Species), linetype = "dotted") +
    theme(legend.position = "bottom", text = element_text(size = 20)) +
    labs(x = "M. gaviniae Gene ID", y = "Distance from M. gaviniae", 
         colour = "Species") +
    geom_errorbar(aes(ymin = DistanceN - Std_Errors, 
                      ymax = DistanceN + Std_Errors, colour = Species), 
                  width = 0.2, position = position_dodge(0.05))
}

png("9_2Plots_gaviniae/MG_one_dist_1_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_gav3, ID <= 531))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_2_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_gav3, ID %in% 532:1063))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_3_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_gav3, ID %in% 1064:1595))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_4_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_gav3, ID %in% 1596:2127))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_5_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_gav3, ID %in% 2128:2659))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_6_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_gav3, ID %in% 2660:3191))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_7_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_gav3, ID %in% 3192:3723))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_8_8.png", width = 1000, height = 750)
dist_plot_one(subset(M_gav3, ID %in% 3724:4255))
dev.off()



png("9_2Plots_gaviniae/MG_one_dist_1_4.png", width = 2000, height = 1000)
dist_plot_one(subset(M_gav3, ID %in% 1:1063))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_2_4.png", width = 2000, height = 1000)
dist_plot_one(subset(M_gav3, ID %in% 1064:2127))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_3_4.png", width = 2000, height = 1000)
dist_plot_one(subset(M_gav3, ID %in% 2128:3191))
dev.off()

png("9_2Plots_gaviniae/MG_one_dist_4_4.png", width = 2000, height = 1000)
dist_plot_one(subset(M_gav3, ID %in% 3192:4255))
dev.off()

### Circular plot #########################################################################################################################################
# Mixta calida
M_calida <- read.csv(file = "8Results/M_calida.csv",                      # Reads in the M. calida results
                     stringsAsFactors = FALSE)
M_calida <- subset(M_calida, select = Gene:Results_Other)

M_calida_dist <- read.csv(file = "8Results/M_calida_dist.csv",
                          stringsAsFactors = FALSE)

M_calida$ID <- M_calida_dist$ID

M_calida$Results_Number <- case_when(
  M_calida$Results_Other == "Pantoea_agglomerans" ~ 1,
  M_calida$Results_Other == "Pantoea_septica" ~ 2,
  M_calida$Results_Other == "Erwinia_amylovora" ~ 3,
  M_calida$Results_Other == "Erwinia_tasmaniensis" ~ 4,
  M_calida$Results_Other == "Tatumella_ptyseos" ~ 5,
  M_calida$Results_Other == "Tatumella_saanichensis" ~ 6,
  M_calida$Results_Other == "Citrobacter_freundii" ~ 7,
  M_calida$Results_Other == "Enterobacter_cloacae" ~ 8
)

# Mixta gaviniae #
M_gaviniae <- read.csv(file = "8Results/M_gaviniae.csv",                  # Reads in the M. gaviniae results
                       stringsAsFactors = FALSE)
M_gaviniae <- subset(M_gaviniae, select = Gene:Results_Other)

M_gaviniae_dist <- read.csv(file = "8Results/M_gaviniae_dist.csv",
                            stringsAsFactors = FALSE)

M_gaviniae$ID <- M_gaviniae_dist$ID

M_gaviniae$Results_Number <- case_when(
  M_gaviniae$Results_Other == "Pantoea_agglomerans" ~ 1,
  M_gaviniae$Results_Other == "Pantoea_septica" ~ 2,
  M_gaviniae$Results_Other == "Erwinia_amylovora" ~ 3,
  M_gaviniae$Results_Other == "Erwinia_tasmaniensis" ~ 4,
  M_gaviniae$Results_Other == "Tatumella_ptyseos" ~ 5,
  M_gaviniae$Results_Other == "Tatumella_saanichensis" ~ 6,
  M_gaviniae$Results_Other == "Citrobacter_freundii" ~ 7,
  M_gaviniae$Results_Other == "Enterobacter_cloacae" ~ 8
)

# M. calida plot #
png("9_1Plots_calida/MC_categ_results.png", width = 1000, height = 725)
ggplot(data = M_calida, aes(xmin = ID - 5, xmax = ID, ymin = 0, ymax = Results_Number, fill = Results_Other)) +
  geom_rect() +
  coord_polar() +
  scale_y_continuous(limits = c(0, 8)) +
  theme(legend.position = "right", text = element_text(size = 20), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  labs(fill = "Species")
dev.off()

# M. gaviniae plot
png("9_2Plots_gaviniae/MG_categ_results.png", width = 1000, height = 725)
ggplot(data = M_gaviniae, aes(xmin = ID - 5, xmax = ID, ymin = 0, ymax = Results_Number, fill = Results_Other)) +
  geom_rect() +
  coord_polar() +
  scale_y_continuous(limits = c(0, 8)) +
  theme(legend.position = "right", text = element_text(size = 20), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  labs(fill = "Species")
dev.off()

#
### Kittens ###############################################################################################################################################
showmekittens()


# c("Tatumella_saanichensis", "Citrobacter_freundii",	"Enterobacter_cloacae", "Erwinia_amylovora", "Erwinia_tasmaniensis", 
#                               "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", "Pantoea_septica", "Tatumella_ptyseos")
