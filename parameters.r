## This file is sourced by spacer_db.r, which shoud contain parameters.

####################
## file path and names

## threads used in bowtie searching for off-target sites
threads <- 15
## temporary file directory
tmpdir <- './tmp'
## species information
species.name <- 'Aedes aegypti'
species.id <- '7159'
genome.build <- 'AaegL5.0'
genome.build.accession <- 'NCBI_Assembly:GCF_002204515.2'
genome.build.url <- 'https://www.ncbi.nlm.nih.gov/assembly/GCF_002204515.2'
species.url <- 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7159'
## file path of reference fasta file
fa.file <- 'reference/Aa.genome.fa'
## prefix of bowtie index
index.prefix <- 'reference/Aa.genome'
## file path of reference gff file
gff.file <- 'GCF_002204515.2_AaegL5.0_genomic.gff'
## file path of database to store spacer information
db.file <- "7159.spacer.sqlite"
## file path of database to store TxDb by R GenomicFeatures package
TxDb.file <- '7159.TxDb.sqlite'

####################
## CRISPR
spacer.length <- 20
pam.length <- 3
pam <- 'NGG'
####################
