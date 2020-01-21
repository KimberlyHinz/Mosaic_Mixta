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

close_relative <- function(gen, spcs) {
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
  
  M_cal <- as.numeric(rbind(t(dist[6, 1:6]), dist[7, 6], dist[8, 6], 
                            dist[9, 6], dist[10, 6]))
  M_gav <- as.numeric(rbind(t(dist[7, 1:7]), dist[8, 7], dist[9, 7],
                            dist[10, 7]))
  
  MCclorel <- close_relative(gene, M_cal)
  MGclorel <- close_relative(gene, M_gav)
  
  M_calida <- rbind(M_calida, MCclorel)
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

write.csv(x = M_calida, file = "8Results/M_calida.csv")
write.csv(x = M_gaviniae, file = "8Results/M_gaviniae.csv")

M_calida <- read.csv(file = "8Results/M_calida.csv")
M_gaviniae <- read.

for(row in 1:nrow(M_calida)) {
  if(M_calida$One[row] %in% c("Mixta_calida", "Mixta_gaviniae")) {
    if(M_calida$Two[row] %in% c("Mixta_calida", "Mixta_gaviniae")) {
      M_calida$Results[row] <- as.character(M_calida$Three[row])
    } else {
      M_calida$Results[row] <- as.character(M_calida$Two[row])
    }
  } else {
    M_calida$Results[row] <- as.character(M_calida$One[1])
  }
}
if(M_calida$One[1] %in% c("Mixta_calida", "Mixta_gaviniae")) {
  if(M_calida$Two[1] %in% c("Mixta_calida", "Mixta_gaviniae")) {
    test <- as.character(M_calida$Three[1])
  } else {
    test <- as.character(M_calida$Two[1])
  }
} else {
  test <- as.character(M_calida$One[1])
}





write.csv(x = M_calida, file = "7Results/M_calida.csv")                   # Write results into csv
write.csv(x = M_gaviniae, file = "7Results/M_gaviniae.csv")               # Write results into csv

### Kittens ###############################################################################################################################################
showmekittens()


# c("Tatumella_saanichensis", "Citrobacter_freundii",	"Enterobacter_cloacae", "Erwinia_amylovora", "Erwinia_tasmaniensis", 
#                               "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", "Pantoea_septica", "Tatumella_ptyseos")
