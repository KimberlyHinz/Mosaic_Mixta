## I used the terminal in BioLinux to change the extensions from .fna to .fasta
#       for f in *.fna
#       do
#       [ -f "$f" ] && mv "$f" "${f%fna}fasta"
#       done

library("seqinr")
library("msa")
library("ape")
library("adegenet")
library("Rfast")
library("beepr")
library("dplyr")

# If on my own computer:
setwd("D:/Honours/")
fastaFiles <- as.data.frame(list.files(path = "C:/Users/Kim/Documents/School/2019_3Fall/Biology_498/Mosaic_Mixta/Genes", 
                                       pattern = ".fasta"))               # Makes a dataframe listing the fasta files in the folder
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("C:/Users/Kim/Documents/School/2019_3Fall/Biology_498/Mosaic_Mixta/Genes/", 
                              fastaFiles$File_name, sep = "")             # Creates a file pathway for each gene

# If in the office:
setwd("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/")
fastaFiles <- as.data.frame(list.files(path = "3Homologues_10/",
                                       pattern = ".fasta"))               # Makes a dataframe where the first column is a list of fasta gene files
colnames(fastaFiles) <- "File_name"                                       # Changes the column name
fastaFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/3Homologues_10/", 
                              fastaFiles$File_name, sep = "")             # Creates a file pathway for each gene

### Selection for genes ###################################################################################################################################
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
# Using ClustalW in R
fastaFilesOrg <- as.data.frame(list.files(path = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/4Organize/", 
                                          pattern = ".fasta"))            # Makes a dataframe listing the fasta files in the folder
colnames(fastaFilesOrg) <- "File_name"                                    # Changes the column name
fastaFilesOrg$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/4Organize/", 
                                 fastaFilesOrg$File_name, sep = "")      # Creates a file pathway for each gene

for(row in 1:nrow(fastaFilesOrg)) {                                       # Aligns genes that passed the filter
  path <- fastaFilesOrg$Path_name[row]
  gene_file <- readDNAStringSet(path)                                     # Reads in the genes (this way ensures the names are kept)
  
  align <- msa::msaClustalW(inputSeqs = gene_file, maxiters = 100, 
                            type = "dna", order = "input")                # Aligns the genes
  alignConv <- msaConvert(align, type = "seqinr::alignment")              # Converts the aligned genes into a readable form
  
  alignFA <- base::as.data.frame(matrix(ncol = 2, nrow = 10))             # Create the necessary format for write.fasta()
  colnames(alignFA) <- c("sequences", "species")
  alignFA$sequences <- alignConv$seq                                      # Copy aligned genes into new dataframe
  alignFA$species <- alignConv$nam
  
  write.fasta(sequences = as.list(alignFA$sequences),
              names = alignFA$species,
              file.out = paste("5Aligned/", fastaFilesOrg$File_name[row],
                               sep = ''),
              open = "w", nbchar = 10000, as.string = TRUE)               # Creates a fasta file for each gene, will continue on to distance matrices
}
beep(8)
rm(align, alignConv, alignFA, gene_file, path, row)

### Best model for genes ##################################################################################################################################
fastaFilesModel <- as.data.frame(list.files(path = "D:/Honours/6Model/",
                                            pattern = ".csv"))            # Make a dataframe listing the csv files in the folder
colnames(fastaFilesModel) <- "File_name"                                  # Changes the column name
fastaFilesModel$Path_name <- paste("D:/Honours/6Model/",
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

### Distance matrices #####################################################################################################################################
megFiles <- as.data.frame(list.files(path = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/7Distance/",
                                     pattern = ".meg"))                   # Makes a list of all .meg file in this diretory
colnames(megFiles) <- "File_name"                                         # Changes the column name
megFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/7Distance/",
                            megFiles$File_name, sep = "")                 # Adds the path name for each gene
megFiles$Gene <- best_model$Gene                                          # Adds the gene name (no extension)

relative <- function(number) {                                            # Returns the relative name given row number
  if(number == 1) {
    rltv <- "Tatumella_saanichensis"
  } else if(number == 2) {
    rltv <- "Citrobacter_freundii"
  } else if(number == 3) {
    rltv <- "Enterobacter_cloacae"
  } else if(number == 4) {
    rltv <- "Erwinia_amylovora"
  } else if(number == 5) {
    rltv <- "Erwinia_tasmaniensis"
  } else if(number == 6) {
    rltv <- "Mixta_calida"
  } else if(number == 7) {
    rltv <- "Mixta_gaviniae"
  } else if(number == 8) {
    rltv <- "Pantoea_agglomerans"
  } else if(number == 9) {
    rltv <- "Pantoea_septica"
  } else if(number == 10) {
    rltv <- "Tatumella_ptyseos"
  }
}

relatives <- function()
for(row in 1:nrow(megFiles)) {                                            # Finds the two closest relatives to Mixta calida
  path <- megFiles$Path_name[row]                                         # Path name
  gene <- megFiles$Gene[row]                                              # Gene name
  
  if(gene %in% c("37818_hypothetical_protein", "38262_ygbE", 
                 "38956_hypothetical_protein", "39709_yciH", 
                 "39916_eamA")) {                                         # These five genes had to be run manually (therefore different format)
    mega <- read.table(file = path, stringsAsFactors = FALSE, skip = 37, fill = TRUE)
  } else {                                                                # For the rest
    mega <- read.table(file = path, stringsAsFactors = FALSE, skip = 45, fill = TRUE)
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
  colnames(dist) <- c("Tatumella_saanichensis", "Citrobacter_freundii",	
                      "Enterobacter_cloacae", "Erwinia_amylovora", 
                      "Erwinia_tasmaniensis", "Mixta_calida", 
                      "Mixta_gaviniae", "Pantoea_agglomerans", 
                      "Pantoea_septica", "Tatumella_ptyseos")             # Changes the column names

  MCmin1 <- which(dist$Mixta_calida == nth(x = dist$Mixta_calida, k = 1, descending = FALSE, arr.ind = TRUE))
}
rm(dist, mega, mega2, gene, path)



min2 <- which(dist_matrix$Mixta_calida ==                               # Finds second closest relative (first is itself)
                nth(x = dist_matrix$Mixta_calida, k = 2, 
                    descending = FALSE), arr.ind = TRUE)                # k = 2 is second smallest number
min3 <- which(dist_matrix$Mixta_calida ==                               # Finds third closest relative
                nth(x = dist_matrix$Mixta_calida, k = 3,
                    descending = FALSE), arr.ind = TRUE)                # k = 3 is third smallest number

rela <- data.frame(gene = csvFiles$Gene[row],                           # Dataframe with first and second relative
                   two = relative(min2),
                   third = relative(min3))

M_calida <- rbind(M_calida, rela)



# c("Tatumella_saanichensis", "Citrobacter_freundii",	"Enterobacter_cloacae", "Erwinia_amylovora", "Erwinia_tasmaniensis", 
#                               "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", "Pantoea_septica", "Tatumella_ptyseos")

### Categorizing closest relative #########################################################################################################################
csvFiles <- as.data.frame(list.files(path = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/6Distance/",
                                     pattern = ".csv"))                   # Makes a dataframe listing the csv files in the folder
colnames(csvFiles) <- "File_name"                                         # Changes the column name
csvFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/6Distance/",
                            csvFiles$File_name, sep = "")                 # Creates a file pathway for each distance matrix
csvFiles$Gene <- gsub(pattern = ".csv", replacement = "", 
                      x = csvFiles$File_name)                             # Creates a column with just the gene name

relative <- function(number) {                                            # Returns the relative name given row number
  if(number == 1) {
    rltv <- "Tatumella_saanichensis"
  } else if(number == 2) {
    rltv <- "Citrobacter_freundii"
  } else if(number == 3) {
    rltv <- "Enterobacter_cloacae"
  } else if(number == 4) {
    rltv <- "Erwinia_amylovora"
  } else if(number == 5) {
    rltv <- "Erwinia_tasmaniensis"
  } else if(number == 6) {
    rltv <- "Mixta_calida"
  } else if(number == 7) {
    rltv <- "Mixta_gaviniae"
  } else if(number == 8) {
    rltv <- "Pantoea_agglomerans"
  } else if(number == 9) {
    rltv <- "Pantoea_septica"
  } else if(number == 10) {
    rltv <- "Tatumella_ptyseos"
  }
}

M_calida <- as.data.frame(matrix(ncol = 3, nrow = 0))                     # Dataframe for M. calida's closest relatives
colnames(M_calida) <- c("Gene", "First_relative", "Second_relative")      # Changes the column names

for(row in 1:nrow(csvFiles)) {                                            # Finds the two closest relatives to Mixta calida
  path <- csvFiles$Path_name[row]
  dist_matrix <- read.csv(file = path, header = TRUE)                     # Reads in the csv file
  
  min2 <- which(dist_matrix$Mixta_calida ==                               # Finds second closest relative (first is itself)
                  nth(x = dist_matrix$Mixta_calida, k = 2, 
                      descending = FALSE), arr.ind = TRUE)                # k = 2 is second smallest number
  min3 <- which(dist_matrix$Mixta_calida ==                               # Finds third closest relative
                  nth(x = dist_matrix$Mixta_calida, k = 3,
                      descending = FALSE), arr.ind = TRUE)                # k = 3 is third smallest number
  
  rela <- data.frame(gene = csvFiles$Gene[row],                           # Dataframe with first and second relative
                     two = relative(min2),
                     third = relative(min3))
  
  M_calida <- rbind(M_calida, rela)
}
rm(dist_matrix, rela, min2, min3, path, row)

for(row in 1:nrow(M_calida)) {                                            # Writes the closest, non-Mixta relative into Results column
  if(M_calida$two[row] == "Mixta_gaviniae") {
    M_calida$Results[row] <- as.character(M_calida$third[row])
  } else {
    M_calida$Results[row] <- as.character(M_calida$two[row])
  }
}
rm(row)

M_gaviniae <- as.data.frame(matrix(ncol = 3, nrow = 0))                   # Dataframe for M. gaviniae's closest relatives
colnames(M_gaviniae) <- c("Gene", "First_relative", "Second_relative")    # Changes the column names

for(row in 1:nrow(csvFiles)) {                                            # Finds the two closest relatives to Mixta gaviniae
  path <- csvFiles$Path_name[row]
  dist_matrix <- read.csv(file = path, header = TRUE)                     # Reads in the csv file
  
  min2 <- which(dist_matrix$Mixta_gaviniae ==                             # Finds second closest relative (first is itself)
                  nth(x = dist_matrix$Mixta_gaviniae, k = 2, 
                      descending = FALSE), arr.ind = TRUE)                # k = 2 is second smallest number
  min3 <- which(dist_matrix$Mixta_gaviniae ==                             # Finds third closest relative
                  nth(x = dist_matrix$Mixta_gaviniae, k = 3,
                      descending = FALSE), arr.ind = TRUE)                # k = 3 is third smallest number
  
  rela <- data.frame(gene = csvFiles$Gene[row],                           # Dataframe with first and second relative
                     two = relative(min2),
                     third = relative(min3))
  
  M_gaviniae <- rbind(M_gaviniae, rela)
}
rm(rela, min2, min3, path, row)

for(row in 1:nrow(M_gaviniae)) {                                          # Writes the closest, non-Mixta relative into Results column
  if(M_gaviniae$two[row] == "Mixta_calida") {
    M_gaviniae$Results[row] <- as.character(M_gaviniae$third[row])
  } else {
    M_gaviniae$Results[row] <- as.character(M_gaviniae$two[row])
  }
}
rm(row)

write.csv(x = M_calida, file = "7Results/M_calida.csv")                   # Write results into csv
write.csv(x = M_gaviniae, file = "7Results/M_gaviniae.csv")               # Write results into csv

### Kittens ###############################################################################################################################################
showmekittens()
