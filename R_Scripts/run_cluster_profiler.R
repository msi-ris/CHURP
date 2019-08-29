#!/usr/bin/env Rscript


#######################################################################
# Install and load required R packages
#######################################################################

req_pkgs <- list("tidyr","openxlsx","clusterProfiler","DOSE","RColorBrewer","glue","msigdbr","pathview","readr","stringr")

for (i in req_pkgs){
  if (i %in% rownames(installed.packages())){
    library(i,character.only = TRUE)
  }else{
    BiocManager::install(i,character.only = TRUE)
    library(i,character.only = TRUE)
  }
}


# show loaded packages
#(.packages())
#unloadNamespace("org.Cf.eg.db")


#######################################################################
# Parse options
#######################################################################

# Clear out any objects from workspace
rm(list = ls(all.names = TRUE))

library("optparse")
option_list = list(
  make_option(c("-f", "--file"), 
              type="character", 
              default=NULL, 
              help="DE expression result file name",
              metavar="character"),
  make_option(c("-s","--species"),
              type="character",
              default=NULL,
              help="species name from mSigDB; underscores should replace spaces. eg. Canis_lupus_familiaris ", 
              metavar="character"),
  make_option(c("-d","--dir"),
              type="character",
              default=NULL,
              help="project directory",metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt


setwd(opt$dir)

# show species that are available in the msigdb
msigdbr_show_species()
# insert species here
species = opt$species
species
species = gsub("_"," ",species)
species
#species = "Mus musculus"

#------------species in molecular signatures databases-----------#
dbs <- read.delim("species_databases.txt",sep="\t",row.names = "msig_species")
# shows all species
print(rownames(dbs))
# save msig species library as "orgdb"
orgdb <- as.character(dbs[rownames(dbs)==species,]$db)
orgdb
# save kegg species identifier as "kegg_sp" (used in )
kegg_sp <- as.character(dbs[rownames(dbs)==species,]$KEGG_species)
kegg_sp


# Install and load the Entrez species genome annotation library for the species; or load it if it is already installed
if (orgdb %in% rownames(installed.packages())){
  library(orgdb,character.only = TRUE)
} else {
  BiocManager::install(orgdb,character.only = TRUE)
  library(orgdb,character.only = TRUE)
}


# Load the edgeR DE results
d1 <- read.delim(opt$file,header = T, sep = '\t', as.is = T)
colnames(d1)[1] <- "gene"
# Add Entrez IDs
d1_ids <- bitr(d1$gene, fromType = 'ENSEMBL', toType=c('ENTREZID'), OrgDb=orgdb)
d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = 'ENSEMBL')
names(d1_merged)[2] <- 'log2FC'


# Set R working directory
project_dir <- opt$dir
setwd(project_dir)


# ---------------------------------------------------------------------
# Read in R objects
# ---------------------------------------------------------------------


# load custom functions
source(glue("{project_dir}/tk_custom_functions_mmedit.R"))

# create a directory for species msigdb
out_dir <- "034_msigdb"
prefix <- "034_"
if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
source("tk_custom_functions_mmedit.R")

# retrieve all gene sets associated with species
msigdb <- msigdbr(species = species)
# convert colons to underscores in gene set subcategory column
colon_to_underscore <- gsub(":", "_", msigdb$gs_subcat)
database <- tolower(paste0(msigdb$gs_cat, "_", colon_to_underscore))
database <- sub("_$","",database) # MM modification -- this is needed; otherwise the "hallmark gene set" shows 0 genes when you run the cluster profiler wrapper function
table(database)
msigdb$database <- database

# Gene sets retrived
#H: hallmark gene sets
#C1: positional gene sets
#C2: curated gene sets
#C3: motif gene sets
#C4: computational gene sets
#C5: GO gene sets
#C6: oncogenic signatures
#C7: immunologic signatures


# Add the links
msigdb$description <- paste0("http://www.broadinstitute.org/gsea/msigdb/cards/", msigdb[[1]])
# Create a data frame for clusterProfiler (for genes as Entrez Gene IDs).
msigdb_long <- msigdb %>%
  dplyr::select(gs_name, entrez_gene, human_gene_symbol, database, description) %>% 
  as.data.frame()
# change column names
colnames(msigdb_long) <- c("ont", "gene", "symbol", "database", "description")
write_tsv(msigdb_long, paste0(out_dir, "/", prefix, "msigdb_long.tsv"))
msigdb_long$database <- factor(msigdb_long$database)


# Set up output directory and filename prefix
out_dir <- glue("{project_dir}/cluster_profiler_results")

# Create output directory
if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}


# ---------------------------------------------------------------------
# Run cluster profiler wrapper function
# ---------------------------------------------------------------------


setwd(out_dir)
cluster_profiler <- tk_cluster_profiler(de = d1_merged, 
      log2fc_cutoff = 0.2, 
      padj_cutoff = 0.05, 
      logFC_col = which(colnames(d1_merged)=="log2FC"),
      padj_col = which(colnames(d1_merged)=="FDR"), 
      ensembl_id_col = which(colnames(d1_merged)=="gene"), 
      entrez_id_col = which(colnames(d1_merged)=="ENTREZID"), 
      msigdb_long = msigdb_long, 
      database_types = c("h", "c2_cp_kegg", "c2_cp_reactome","c7"), 
      #plot_title = "", 
      run_ora = TRUE, 
      run_go = F, 
      run_kegg = TRUE, 
      run_gsea = TRUE)
setwd(project_dir)











