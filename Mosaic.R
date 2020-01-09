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

# If on my own computer:
setwd("C:/Users/Kim/Documents/School/2019_3Fall/Biology_498/Mosaic_Mixta/")
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

best_model <- as.data.frame(matrix(ncol = 3, nrow = 0))                   # Dataframe for each gene's best model
for(row in 1:nrow(fastaFilesModel)) {
  path <- fastaFilesModel$Path_name[row]
  gene_model_test <- read.csv(file = path)
  
  model <- as.data.frame(gene_model_test$Model[1])
  colnames(model) <- "Model"
  model$File_name <- fastaFilesModel$File_name[row]
  model$ModelA <- "hello"
  if(model$Model == "T92+G+I") {
    model$ModelA == "T92_G"
  } else if (model$Model == "TN93+G+I") {
    model$ModelA == "TN93_G"
  } else if (model$Model == "T92+G") {
    model$ModelA == "T92_G"
  } else if (model$Model == "GTR+G+I") {
    model$ModelA == "GTR_G"
  } else if (model$Model == "HKY+G+I") {
    model$ModelA == "HKY_G"
  } else if (model$Model == "K2+G") {
    model$ModelA == "K2_G"
  } else if (model$Model == "TN93+G") {
    model$ModelA == "TN93_G"
  } else if (model$Model == "K2+I") {
    model$ModelA == "K2"
  } else if (model$Model == "TN93+I") {
    model$ModelA == "TN93"
  } else if (model$Model == "GTR+G") {
    model$ModelA == "GTR_G"
  } else if (model$Model == "K2+G+I") {
    model$ModelA == "K2_G"
  } else if (model$Model == "HKY+G") {
    model$ModelA == "HKY_G"
  } else if (model$Model == "JC+G") {
    model$ModelA == "JC_G"
  } else if (model$Model == "T92+I") {
    model$ModelA == "T92"
  }
  
  best_model <- rbind(best_model, model)
}
rm(gene_model_test, model, path, row)

best_model$Gene <- gsub(pattern = "-4212.csv", replacement = "", 
                             x = best_model$File_name)                    # Creates a column with just the gene name
best_model$Path_Name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/5Aligned/",
                              best_model$Gene, ".fasta", sep = "")

for(row in 1:nrow(best_model)) {
  
}



Uniq_mods <- as.data.frame(unique(best_model$Model))
colnames(Uniq_mods) <- "Model_Para"
Uniq_mods$Model_Name <- c("T92_G_I", "TN93_G_I", "T92_G", "GTR_G_I", 
                              "HKY_G_I", "K2_G", "TN93_G", "K2_I", 
                              "TN93_I", "GTR_G", "K2_G_I", "HKY_G", 
                              "JC_G", "T92_I")
for(row in 1:nrow(Uniq_mods)) {
  Para <- Uniq_mods$Model_Para[row]
  Name <- Uniq_mods$Model_Name[row]
  
  datframe <- subset(best_model, Model == Para)
  datframe <- datframe$Path_Name
  
  write.table(datframe, file = paste(Name, ".txt", sep = ""), sep = "\n",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
rm(datframe, Name, Para, row)

### Distance matrices #####################################################################################################################################
fastaFilesAlign <- as.data.frame(list.files(path = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/5Aligned/",
                                            pattern = ".fasta"))          # Makes a dataframe listing the fasta files in the folder
colnames(fastaFilesAlign) <- "File_name"                                  # Changes the column name
fastaFilesAlign$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/5Aligned/",
                                   fastaFilesAlign$File_name, sep = "")   # Creates a file pathway for each gene
fastaFilesAlign$Gene <- gsub(pattern = ".fasta", replacement = "", 
                             x = fastaFilesAlign$File_name)               # Creates a column with just the gene name


for(row in 1:nrow(fastaFilesAlign)) {                                     # Creates a distance matrix for each aligned gene
  path <- fastaFilesAlign$Path_name[row]
  gene_file <- fasta2DNAbin(file = path)                                  # Reads in fasta file as a DNAbin
  
  gene_dist <- dist.dna(x = gene_file, model = "T92",                     # Tamura 1992 model
                        pairwise.deletion = TRUE,                         # Does pairwise deletions
                        as.matrix = TRUE)                                 # Returns matrix format
  
  gene_dist_df <- as.data.frame(as.matrix(gene_dist))                     # Converts to dataframe
  colnames(gene_dist_df) <- c("Tatumella_saanichensis", "Citrobacter_freundii",	"Enterobacter_cloacae", "Erwinia_amylovora", "Erwinia_tasmaniensis", 
                              "Mixta_calida", "Mixta_gaviniae", "Pantoea_agglomerans", "Pantoea_septica", "Tatumella_ptyseos")
  # Change column names (useful for next step)
  write.csv(x = gene_dist_df, 
            file = paste("6Distance/", fastaFilesAlign$Gene[row], 
                         ".csv", sep = ""))                               # Write distance matrices into csv
}
rm(gene_dist, gene_dist_df, gene_file, path, row)

### Categorizing closest relative #########################################################################################################################
csvFiles <- as.data.frame(list.files(path = "C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/6Distance/",
                                     pattern = ".csv"))                   # Makes a dataframe listing the csv files in the folder
colnames(csvFiles) <- "File_name"                                         # Changes the column name
csvFiles$Path_name <- paste("C:/Users/officePC/Documents/Kim_Honours/Mixta_Mosaic/6Distance/",
                            csvFiles$File_name, sep = "")                 # Creates a file pathway for each distance matrix
csvFiles$Gene <- gsub(pattern = ".csv", replacement = "", 
                      x = csvFiles$File_name)                      # Creates a column with just the gene name

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
