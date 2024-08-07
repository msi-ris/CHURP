############################
# CHURP 1.0.0 bulk RNA-seq summary script
# See run_summary_stats.sh for additional run context.
# See https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf for the edgeR user manual.
# Contact help@msi.umn.edu for questions
############################


# Prepend the library path to libPaths
# Load libraries
library('limma')
library('edgeR')
library('ggplot2')
library('reshape2')
library('gplots')
library('gtools')
library('grid')
library('readxl')
library('tools')

#grab the working and output directories, as well as the sample sheet,
# and merged raw counts matrix, and the groupsheet
args <- commandArgs(trailingOnly = T)
out_dir <- args[1]
work_dir <- args[2]
samp_sheet <- args[3]
fc_mat <- args[4]
min_len <- as.numeric(args[5])
min_cts <- as.numeric(args[6])
group_sheet_loc <- args[7]

setwd(work_dir)

############################
# Read in data and prep relevant data and set output files
############################

# Get the sample sheet to grab group membership downstream 
sample_sheet <- read.table(samp_sheet, sep = "|", header = F, comment.char = "#")

# parse the excel spreadsheet. The first sheet has group (and batch)
# information. The second sheet, if present, has DEG testing contrast information.
# the file itself should always exist, because it is tested for earlier. However,
# there hasn't yet been a test for the DEG groups (sheet 2), so here we'll pay
# careful attention to that.
group_sheet <- readxl::read_excel(group_sheet_loc,
                                  sheet =1)
comparison_sheet <- tryCatch({
  readxl::read_excel(group_sheet_loc,
                     sheet =2)
}, 
error = function(e){
  NULL
}
)
if (is.null(comparison_sheet) || nrow(comparison_sheet) == 0){
  write(
    paste0("No contrasts found in sheet 2 of ",
           basename(group_sheet_loc),
           ", will not perform DEG testing."),
    stderr()
  )
  do_deg = FALSE
}else{
  do_deg = TRUE
}


# Because there may be cases where a subset of individuals in the samplesheet are run. We'll pull in the featureCounts matrix early and grab the relevant IDs
raw_mat <- read.table(fc_mat, header = T, sep = '\t', comment.char = '#')
samp_ids <- names(raw_mat)[-(1:6)]
sample_sheet <- sample_sheet[make.names(sample_sheet$V1) %in% samp_ids,]
group_sheet <- group_sheet[make.names(group_sheet$SampleName) %in% samp_ids,]
group_sheet <- group_sheet[match(samp_ids, make.names(group_sheet$SampleName)),]

# Now we can proceed with group determinations. We have to call make.names()
# for the groups because level names must be valid R object names.
# The call to make.names() converts the 'NULL' group designation to 'NULL.', which 
# we have to account for in our true_groups assignments.
groups <- make.names(as.vector(group_sheet$Group))
uniq_groups <- unique(groups)
true_groups <- groups[groups != 'NULL.']
n_true_groups <- length(unique(true_groups))

# We optionally allow the user to add a third column to their group sheet, which
# is assumed to be a batch variable. Here, we check for the third column, and 
# whether it is categorical-like. We do a simple check to ensure there are fewer
# batch categories than their are total samples
has_batch_variable <- FALSE
if (ncol(group_sheet) > 2){
  batch_name <- colnames(group_sheet)[3]
  batches <- make.names(as.vector(group_sheet[,3]))
  uniq_batches <- unique(batches)
  if (length(uniq_batches) < nrow(group_sheet)){
    has_batch_variable <- TRUE
  }
}


# Then check if the number of replicates meets our criteria.
if (length(true_groups[which(table(true_groups) < 3)]) > 0){
  write("At least one of the groups has fewer than three replicates, which makes reliable statistical interpretation difficult and therefore DEG testing will not be done.", stderr())
  writeLines(text = paste("This group has fewer than three replicates", unique(groups)[which(table(groups) < 3)], sep=" : "), con = stderr())
  do_deg <- FALSE
}

# Set filename variables
mds_plot <- paste(out_dir, "Plots/mds_plot.pdf", sep = "/")
counts_plot <- paste(out_dir, "Plots/cpm_plot.pdf", sep = "/")
counts_list <- paste(out_dir, "Counts/cpm_list.txt", sep = "/")
hmap <- paste(out_dir, "Plots/high_variance_heatmap.pdf", sep = "/")

# Filter out genes that are below the length threshold
raw_mat <- raw_mat[which(raw_mat$Length >= min_len),]

# Check for any library sizes of zero and exit with 1 if found
lib_sizes <- colSums(raw_mat[,seq(-1,-6)])
if (any(lib_sizes == 0)){
  write(paste0("summarize_bulk_rnaseq.R: ERROR\n",
               "The following samples had zero counts after removing genes shorter than min_len or with zero variance, exiting early:\n", 
               paste(samp_ids[lib_sizes == 0], collapse = ","),
               "\nPerhaps your GTF doesn't match your species?"),
        stderr())
  quit(status = 1, save = "no")
}

# Convert the raw matrix into a DGE object. Column 1 is Geneid, columns 7+ are the sample counts. The groups list  is generated above. 
edge_mat <- DGEList(counts = raw_mat[,seq(-1,-6)], genes = raw_mat[,1], group = groups)

############################
# Generate descriptive accounts of the data (MDS, Normalized Counts, Counts Distributions, and a Heatmap) 
# Note that we DO NOT apply the minimum count filters prior to these descriptive summaries.
############################

# This is THE ONLY colorblind acceptable palette with 4 colors from colorbrewer: http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=4
if(length(uniq_groups) > 8) {
    write("There are more than eight groups, so all groups will be colored dark blue. We plot a maximum of eight colors for colorblindness considerations.", stderr())
    pal <- rep("navy", length.out=length(uniq_groups))
} else {
  # Colorblind friendly palette
  pal <- c("#5330A0","#F21245","#1E88E5","#FFC107",
           "#7AF3DE","#81C5E6","#004D40","#A43CB3")
}
col_vec <- pal[match(groups,uniq_groups)]

# legend code adapted from https://support.bioconductor.org/p/101530/
# Set the MDS plot pdf and write the plot
pdf(mds_plot)
# If there are fewer than 3 samples, we still want to write a PDF for the MDS,
# but we can't actually generate an MDS plot. The edgeR function dies with
# less than 3 samples
if(length(samp_ids) < 3) {
  plot(c(0, 1), c(0, 1), ann=F, bty="n", type="n", xaxt="n", yaxt="n")
  text(x=0.5, y=0.5, "Less than 3 samples;\nMDS not possible", cex=1, col="black")
} else if (has_batch_variable){
  pchs <- 21:(20+length(uniq_batches))
  pch_vec <- pchs[match(batches, uniq_batches)]
  
  opar <- par(no.readonly = TRUE)
  par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
  plotMDS(edge_mat, cex = 0.75, col = col_vec, bg= col_vec, pch = pch_vec)
  # we do this to get the x,y locs, because we can't plot labels and points in a single plotMDS call
  xy <- as.data.frame(plotMDS(edge_mat, cex = 0.75, plot = FALSE))
  text(xy$x, xy$y, label = row.names(xy), pos = 1, cex = 0.5)
  legend(par("usr")[2], mean(par("usr")[3:4])+.25, legend = c('Group', uniq_groups), text.col = c('black', unique(col_vec)), bty = "n")
  legend(par("usr")[2], mean(par("usr")[3:4])-.25, legend = c(batch_name, uniq_batches), pch = c(26, unique(pch_vec)), bty = "n")
  
  par(opar)
} else {
  opar <- par(no.readonly = TRUE)
  par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
  plotMDS(edge_mat, cex = 0.75, col = col_vec)
  legend(par("usr")[2], mean(par("usr")[3:4]), legend = c('Group', uniq_groups), text.col = c('black', unique(col_vec)), bty = "n")
  par(opar)
}
dev.off()


# Set a variable holding the log2(1+CPM) counts.
cpm_counts <- cpm(edge_mat, log = T, prior.count = 1)

# Create a dataframe of cpm_counts and gene IDs in 'wide' format. Melt into 'long' format for ggplot
cdf <- data.frame(edge_mat$genes, cpm_counts)
write.table(cdf, file = counts_list, sep = '\t', quote = FALSE, row.names = FALSE)
tidy_cdf <- melt(cdf, id.vars = "genes", variable.name = "sample_id", value.name = "per_feature_count")

# manually add group information to tidy_cdf
tidy_cdf$group <- factor(rep(uniq_groups[match(groups,uniq_groups)], each = nrow(edge_mat$genes)))

# Set the counts plot pdf and write the violin plot of normalized counts per sample
pdf(counts_plot)
p <- ggplot(tidy_cdf, aes(x = sample_id, y = per_feature_count, fill = group)) + 
  geom_violin(trim = F) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7)) + 
  scale_fill_manual("Group", values = pal)
p + labs(x = "Sample ID", y = "Feature count -- log(1+cpm)", fill = "Group")
dev.off()


# Calculate count variance across samples and select the top 500 variance features.
#cpm_counts <- cpm(edge_mat, log = T, prior.count = 1)
# For some reason, sometimes there are fewer than 500 genes that pass filtering
n_genes <- min(500, nrow(cpm_counts))
if(n_genes < 500) {
  write("There are fewer than 500 genes that pass variance filtering for the clustering heatmap. This is not an error, but you should be aware of it.", stderr())
}
if(length(samp_ids) == 1) {
  write("There is only one sample, so we will not try to generate a clustering heatmap. This is not an error.", stderr())
  pdf(hmap)
  plot(c(0, 1), c(0, 1), ann=F, bty="n", type="n", xaxt="n", yaxt="n")
  text(x=0.5, y=0.5, "1 sample;\nClustering heatmap not possible", cex=1, col="black")
  dev.off()
} else {
  gene_var <- apply(cpm_counts, 1, var)
  # Need to run a check here to see that they are not all 0 variance
  if(all(gene_var == 0)) {
    write("All genes have 0 variance, so we will not try to generate a clustering heatmap. This is not an error.", stderr())
    pdf(hmap)
    plot(c(0, 1), c(0, 1), ann=F, bty="n", type="n", xaxt="n", yaxt="n")
    text(x=0.5, y=0.5, "All genes have 0 variance in expression;\nClustering heatmap not possible", cex=1, col="black")
    dev.off()
  } else {
    select_var <- names(sort(gene_var, decreasing=TRUE))[1:n_genes]
    high_var <- cpm_counts[select_var,]
    # Set the heatmap pdf and plot the normalized counts heatmap
    pdf(hmap)
    # different cases for when there are group variables or not
    if (n_true_groups > 0){
      # pheatmap uses a dataframe for variable annotation where the rownames match the matrix samplenames
      annotation <- as.data.frame(group_sheet)
      row.names(annotation) <- make.names(group_sheet$SampleName)
      annotation[,1] <- NULL
      # to specify annotation colors for pheatmap, use a named list with named color vector
      colors <- col_vec
      colors <- unique(colors)
      # A quick fix - if the length of the color vector is 1, then we
      # either have 1 group or >4 groups. We will overwrite the color
      # vector in this case
      if(length(colors) == 1) {
        colors <- rep(colors, length(unique(annotation[,1])))
      }
      names(colors) <- unique(annotation[,1])
      
      color_list <- list()
      color_list[[colnames(annotation)[1]]] <- colors
      pheatmap::pheatmap(high_var,
                         treeheight_row = 0,
                         show_rownames = FALSE,
                         annotation_col = annotation,
                         annotation_colors = color_list)
    }else{
      pheatmap::pheatmap(high_var,
                         treeheight_row = 0,
                         show_rownames = FALSE)
    }
    dev.off()
  }
}
############################
# Differential expression testing and summaries
############################
if(!do_deg) {
  quit(status = 0, save = "no")
}
n_groups <- length(uniq_groups)

# If there is only one grouping in the data, we don't need to run the subsequent tests.
if (n_groups == 1){
  write("Only 1 grouping present, skipping differential expression tests.", stderr())
  quit(status = 0, save = "no")
}


# Subset the data object to get rid of samples with a 'NULL' group
edge_mat <- edge_mat[,edge_mat$samples$group %in% true_groups]

# Filter out genes wtih low expression. We employ the following filtering
# scheme, which is similar to what edgeR's `filterByExpr()` function does, but
# with explicit statements:
#   1: Calculate the median library size across all samples (C)
#   2: Calcualte the CPM (K, not on log scale) corresponding to `min_cts` in C
#   3: Calculate the size of the smallest group (G)
#   4: Keep genes where at least G samples have CPM of K.
med_lib <- median(edge_mat$samples$lib.size) / 1000000
min_cpm <- as.numeric(min_cts) / med_lib
min_grp <- min(table(true_groups))
filter_low_expression <- function(gene_row, min_expr, min_samples) {
    num_expr <- sum(as.numeric(gene_row) >= min_expr)
    if(sum(num_expr) >= min_samples) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
keep <- apply(
    cpm(edge_mat, normalized=TRUE, log=FALSE),
    1,
    filter_low_expression,
    min_cpm,
    min_grp)
edge_mat <- edge_mat[keep, ,keep.lib.sizes = FALSE]
# Calculate the normalization factors
edge_mat <- calcNormFactors(edge_mat)

# Print some diagnostic info
print(paste("Median library size in millions of fragments: ", med_lib, sep=""))
print(paste("CPM threshold for 'unexpressed' (not log scale): ", min_cpm, sep=""))
print(paste("Size of smallest group: ", min_grp, sep=""))
print(paste("Number of retained genes: ", nrow(edge_mat), sep=""))

# Generate the design matrix for GLM fitting and estimate common and tag-wise dispersion in one go.
design <- model.matrix(~0+group, data = edge_mat$samples)
edge_mat <- estimateDisp(edge_mat, design = design)

# Fit the per-feature negative binomial GLM.
fit <- glmQLFit(edge_mat, design)


# Check if the Reference and Test Groups listed in the comparison CSV file are 
# present within the edgeR sample groups. If not, skip testing for that comparison.
for (i in 1:dim(comparison_sheet)[1]){
  # check if the groups in the comparison match what is present in the sample sheet
  comparison <- comparison_sheet$Comparison_Name[i]
  ref_group <- comparison_sheet$Reference_Group[i]
  test_group <- comparison_sheet$Test_Group[i]
  if (ref_group %in% true_groups & test_group %in% true_groups){
    comp <- paste0("group",test_group,"-group",ref_group)
    comp_var <- makeContrasts(comp, levels = design)
    qlf <- glmQLFTest(fit, contrast  = comp_var)
    tags <- topTags(qlf, n = nrow(qlf$genes))
    comp <- gsub("group","",comp)
    de_file <- paste(out_dir, "/DEGs/DE_", comp, "_list.txt", sep = "")
    write.table(tags$table, file = de_file, sep = '\t', quote = FALSE, row.names = FALSE)
  }else{
    #print missing a group. or group misspelled
    write(paste0("Missing a group in Comparison: ",comparison,". A Reference and/or Test group does not match the groups listed in the Sample Sheet. Check the spelling of the group names to make sure that they match. "), stderr())
  }
}

