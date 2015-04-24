# Here, we read the table into R. Always remember that when working with text that strings should not be factors! 
tab <- read.table("master.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Now, here we prepare for the cytb gene excluding NAs
reduced_table <- tab[-which(is.na(tab$nucleotides_CYTB)), ]

# set nomenclature in our submission
gene_name <- "CYTB"
prod_name <- "cytochrome b"

# Here we write a GenBank format fasta file containing critical source modifying annotations
fasta_description <- paste0(">", paste0(reduced_table$otherCatalogNumbers, "_", gene_name), " ", "[organism=", reduced_table$genus, " ", reduced_table$specificEpithet, "]", " ", "[Bio_material=", reduced_table$otherCatalogNumbers, "]", " ", "[Specimen-voucher=", reduced_table$institutionCode, ":", reduced_table$catalogNumber, "]", " ", "[location=mitochondrion] [mgcode=2]")
fasta_complete <- paste(fasta_description, reduced_table$nucleotides_CYTB, sep="\n")# add data to fasta
write(fasta_complete, file="sequences.fsa", append=FALSE)# write out the fasta file

# Now we create the feature table containing the locations of attributes of the sequence.
feature_tab <- paste0(paste0(">Feature", " ", reduced_table$otherCatalogNumbers, "_", gene_name),"\n", # 
    "1", "\t", ">", nchar(reduced_table$nucleotides_CYTB), "\t", "gene", "\n", #
    "\t", "\t", "\t", "gene", "\t", gene_name, "\n", #
    "1", "\t", ">", nchar(reduced_table$nucleotides_CYTB), "\t", "CDS", "\t", "\t", "\n", #
    "\t", "\t", "\t", "product", "\t", prod_name, "\n", #
    "\t", "\t", "\t", "codon_start", "\t", "1")#
write(feature_tab, file="features.tbl", append=FALSE)# write out


# Now we repeat for our 16S sample. Remember we are overwriting our previous objects, but it's okay, as these were written to disk.
reduced_table <- tab[-which(is.na(tab$nucleotides_16S)), ]
gene_name <- "rRNA"
prod_name <- "16S ribosomal RNA"

# Note that we deleted the genetic code element, as this is not a coding sequence
# also note that 'append' is now set to 'TRUE' to add the data to previously written files
fasta_description <- paste0(">", paste0(reduced_table$otherCatalogNumbers, "_", gene_name), " ", "[organism=", reduced_table$genus, " ", reduced_table$specificEpithet, "]", " ", "[Bio_material=", reduced_table$otherCatalogNumbers, "]", " ", "[Specimen-voucher=", reduced_table$institutionCode, ":", reduced_table$catalogNumber, "]", " ", "[location=mitochondrion]")
fasta_complete <- paste(fasta_description, reduced_table$nucleotides_16S, sep="\n")# add data to fasta
write(fasta_complete, file="sequences.fsa", append=TRUE)# write out the fasta file
#
feature_tab <- paste0(paste0(">Feature", " ", reduced_table$otherCatalogNumbers, "_", gene_name),"\n", # 
    "<1", "\t", ">", nchar(reduced_table$nucleotides_16S), "\t", gene_name, "\n", #
    "\t", "\t", "\t", "product", "\t", prod_name, "\n") #
write(feature_tab, file="features.tbl", append=TRUE)# write