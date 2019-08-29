#!/usr/bin/env Rscript


#######################################################################
# Todd's custom R functions, all start with "tk_"
#######################################################################






# ---------------------------------------------------------------------
# Run enricher hypergeometric test
# ---------------------------------------------------------------------

tk_enricher_core <- function(geneset, msigdb_long, database_types = levels(msigdb_long$database)) {
	# all databases
	# database_types <- levels(msigdb_long$database)
	# OR MANUAL EDIT TO CHOOSE FOR YOUR FAVORITE DATABASES
	# database_types <- c("c2_cgp", "c2_cp_reactome", "c2_cp_biocarta", "c2_cp_kegg", "h", "c6")
	# database_types <- levels(msigdb_long$database)[grepl("_all", levels(msigdb_long$database))]
	hypergeo_msigdb_df <- data.frame()
	hypergeo_msigdb_list <- list()
	for (i in seq_along(database_types)) {
		# get only one database type at a time
		# Print database name
		print(database_types[i])
		msigdb_long_curr <- filter(msigdb_long, database == database_types[i])
		# Print the number of gene sets in database
		print(length(levels(factor(msigdb_long_curr$ont))))
		hypergeo_msigdb_curr <- enricher(geneset, TERM2GENE = msigdb_long_curr, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2)
		if (!is.null(hypergeo_msigdb_curr)) {
			hypergeo_msigdb_curr <- setReadable(hypergeo_msigdb_curr, orgdb, keyType = "ENTREZID")
		}
		hypergeo_msigdb_list[[database_types[i]]] <- hypergeo_msigdb_curr
		# Clean up results
		if (!is.null(hypergeo_msigdb_curr)) {
			print("appending all list")
			hypergeo_msigdb_curr_df <- hypergeo_msigdb_curr@result
			hypergeo_msigdb_curr_df <- hypergeo_msigdb_curr_df[hypergeo_msigdb_curr_df$Count > 1, ]
			hypergeo_msigdb_curr_df$Description <- NULL
			hypergeo_msigdb_curr_df$gene_set_database <- database_types[i]
			hypergeo_msigdb_df <- rbind(hypergeo_msigdb_df, hypergeo_msigdb_curr_df)
		} else { 
			print("NOT appending to list") 
		}
		# NOTE:
		# if you get some error like this:
		# --> No gene can be mapped....
		# --> Expected input gene ID: 235293,17125,18710,225326,225010,104709
		# --> return NULL...

		# All your genes in your geneset of interest (upregulated) must be present in the gene sets of interest. If you select
		# a list of gene sets, find the total uniqe set of genes available, the genes you're interested in testing ALL must be present
		# in that gene set list.
		# up %in% unique(msigdb_long_curr$gene)
	}
	core_list <- list(hypergeo_msigdb_df = hypergeo_msigdb_df, hypergeo_msigdb_list = hypergeo_msigdb_list)
	return(core_list)
}





# ---------------------------------------------------------------------
# Gene Ontology - Group info
# ---------------------------------------------------------------------

tk_go_group <- function(geneset, go_level = c(2, 3, 4)) {
	if (length(geneset) == 0) {
		print("geneset is empty")
		go_list <- list()
		return(go_list)
	} else {
		ontology <- c("CC", "MF", "BP")
		go_list <- list()
		for (i in seq_along(ontology)) {
			for (j in seq_along(go_level)) {
			name <- paste0("go_", ontology[i], "_", go_level[j])
			go <- groupGO(gene = geneset, OrgDb = orgdb, ont = ontology[i], level = go_level[j], readable = TRUE)
			# Sort results
			go_sorted <- go@result[order(go@result$Count, decreasing = TRUE), ]
			go_list[[name]] <- go_sorted
			}
		}
		return(go_list)
	}
}


# ---------------------------------------------------------------------
# Gene Ontology - Group enrichment
# ---------------------------------------------------------------------

tk_go_enrich <- function(geneset, universe) {
	if (length(geneset) == 0) {
		print("geneset is empty")
		go_list <- list()
		return(go_list)
	} else {
		ontology <- c("CC", "MF", "BP")
		go_list <- list()
		for (i in seq_along(ontology)) {
			name <- paste0("go_", ontology[i])
			go <- enrichGO(gene = geneset, OrgDb = orgdb, ont = ontology[i], readable = TRUE, universe = universe, pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
			# Sort the results
			go_sorted <- go@result[go@result$qvalue <= 0.05, ]
			go_sorted <- go_sorted[order(go_sorted$qvalue, decreasing = FALSE), ]
			go_list[[name]] <- go_sorted
		}
		return(go_list)
	}
}

###### bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
###### enrichMap(bp)






# ---------------------------------------------------------------------
# KEGG
# ---------------------------------------------------------------------


tk_kegg_enrich <- function(geneset, prerank) {
	if (length(geneset) == 0) {
		print("geneset is empty")
		kegg_enrich_list <- list()
		return(kegg_enrich_list)
	} else {
		# Make sure gene set is only entrez_ids, (digits)
		geneset <- geneset[str_detect(geneset, "\\d")]
		# Run KEGG enrichment test
		kegg_enrich <- enrichKEGG(gene = geneset,
							organism = kegg_sp,
							keyType = "kegg",
       						pvalueCutoff = 0, 
       						pAdjustMethod = "BH",
       						minGSSize = 10, 
       						maxGSSize = 500, 
       						qvalueCutoff = 0)
		kegg_enrich <- setReadable(kegg_enrich, OrgDb = orgdb, keyType = "ENTREZID")
		kegg_enrich_all <- kegg_enrich@result
		kegg_enrich_all$database_link <- paste0("https://www.genome.jp/dbget-bin/www_bget?", kegg_enrich_all$ID)
		kegg_enrich_signif <- kegg_enrich_all[kegg_enrich_all$p.adjust < 0.05 & kegg_enrich_all$qvalue < 0.2, ]
		# Make nice names for filenames
		description_clean <- str_replace_all(kegg_enrich_signif$Description, " ", "_")
		# Create entrez_id named vector of log2FC values for all genes
		#entrez_id <- de[[entrez_id_col]]
		#log2FoldChange <- de[[logFC_col]]
		print(head(prerank))
		gene_list_tb <- prerank %>%
						dplyr::select(ENTREZID, log2FC) %>%
						filter(str_detect(ENTREZID, "\\d")) %>%
						arrange(desc(log2FC))
		gene_list <- gene_list_tb[["log2FC"]]
		names(gene_list) <- gene_list_tb[["ENTREZID"]]
		# Print out all the figs
		orig_wd <- getwd()
		results_out_dir <- glue("{orig_wd}/kegg_enrich")
		if (!dir.exists(results_out_dir)) {dir.create(results_out_dir, recursive = TRUE)}
		setwd(results_out_dir)
		# TEMP TEST
		options(bitmapType = 'cairo')
		for (i in seq_along(kegg_enrich_signif$ID)) {
			# Print KEGG plots
			p1 <- pathview(gene.data  = gene_list,
						pathway.id = kegg_enrich_signif$ID[i],
						species    = kegg_sp,
						low = list(gene = "blue", cpd = "green"),
						mid = list(gene = "white", cpd = "white"),
						high = list(gene = "firebrick1", cpd = "yellow"),
						out.suffix = description_clean[i])
			# Print pathview PDFs
			p2 <- pathview(gene.data  = gene_list,
						pathway.id = kegg_enrich_signif$ID[i],
						species    = kegg_sp,
						low = list(gene = "blue", cpd = "green"),
						mid = list(gene = "white", cpd = "white"),
						high = list(gene = "firebrick1", cpd = "yellow"),
						out.suffix = description_clean[i],
						kegg.native = FALSE)
		}
		# Delete intermediate files
		unlink("*.xml")
		unlink("*[[:digit:]].png")
		setwd(orig_wd)
	}
	return(list(kegg_enrich_all = kegg_enrich_all, kegg_enrich_signif = kegg_enrich_signif))
}






# ---------------------------------------------------------------------
# Create pre-ranked gene list for GSEA
# ---------------------------------------------------------------------

tk_create_prerank <- function(tibble, method = NULL, logFC_col = NULL, padj_col = NULL, sort = TRUE) {
	# method must be one of: "pval_sign" or "lfc"
	# Check input parameters
	if (is.null(logFC_col)) {
		print("You need to provide the column index for log fold change values")
		stop()
	}
	if (is.null(padj_col)) {
		print("You need to provide the column index for adjusted p-values")
		stop()
	}
	# Create prerank table
	if (method == "pval_sign") {
		tibble <- add_column(tibble, "log2FoldChange_sign" = sign(tibble[[logFC_col]]))
		tibble <- add_column(tibble, "minuslog10pval" = -log10(tibble[[padj_col]]))
		# If any are Inf, set to very high value (e.g. 300)
		tibble$minuslog10pval[which(tibble$minuslog10pval == "Inf")] <- 300
		# Create the prerank_metric: -log10pval * fold change sign
		tibble <- add_column(tibble, "prerank_metric" = (tibble$minuslog10pval * tibble$log2FoldChange_sign))
		# Remove genes that have p values with "NA"
		tibble <- tibble[!(is.na(tibble$minuslog10pval)), ]
		if (sort == TRUE) {
			tibble <- tibble[order(tibble$prerank_metric, decreasing = TRUE), ]
		}
		return(tibble)
	} else if (method == "lfc") {
		tibble <- tibble %>% 
					add_column(prerank_metric = tibble[[logFC_col]]) %>%
					drop_na(padj_col) %>%
					filter(ENTREZID != "unknown")
		if (sort == TRUE) {
			tibble <- tibble %>% 
					arrange(desc(tibble[[logFC_col]]))
		}	
		return(tibble)
	} else {
		print("You need to provide a method argument that is either pval_sign or lfc")
		stop()
	}
}





# ---------------------------------------------------------------------
# GSEA function
# ---------------------------------------------------------------------

tk_enricher_gsea <- function(sorted_prerank_df, msigdb_long, plot_title = "Plot Title", qvalue_threshold = 0.05, database_types = levels(fmsigdb_long$database)) {

	prerank_gene_list <- sorted_prerank_df$prerank_metric
	names(prerank_gene_list) <- as.character(sorted_prerank_df$ENTREZID)

	# DATABASES
	# all:
	# database_types <- levels(msigdb_long$database)
	# OR MANUAL EDIT TO CHOOSE FOR YOUR FAVORITE DATABASES
	# database_types <- levels(msigdb_long$database)[grepl("_all", levels(msigdb_long$database))]
	# database_types <- c("c2_cgp", "c2_cp_reactome", "c2_cp_biocarta", "c2_cp_kegg", "h", "c6")
	
	
	gsea_msigdb_df <- data.frame()
	gsea_msigdb_list <- list()
	for (i in seq_along(database_types)) {
		# get only one database type at a time
		# Print database name
		print(database_types[i])
		msigdb_long_curr <- filter(msigdb_long, database == database_types[i])
		# Print the number of gene sets in database
		print(length(levels(factor(msigdb_long_curr$ont))))
		gsea_msigdb_curr <- clusterProfiler::GSEA(prerank_gene_list, TERM2GENE = msigdb_long_curr, nPerm = 10000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH")
		if (length(gsea_msigdb_curr$ID) > 0) {
			gsea_msigdb_curr <- setReadable(gsea_msigdb_curr, OrgDb = orgdb, keyType = "ENTREZID")
			gsea_msigdb_curr@result$gene_set_database <- database_types[i]
			gsea_msigdb_curr@result$gene_set_description <- msigdb_long$description[match(gsea_msigdb_curr$ID, msigdb_long$ont)]
			# Sometimes with only 1 signif geneset, the qvalues are logical with NA. This causes problems when plotting
			if (!is.logical(gsea_msigdb_curr@result$qvalues)) {
				gsea_msigdb_curr_signif_df <- gsea_msigdb_curr[gsea_msigdb_curr@result$qvalues <= qvalue_threshold, ]
				gsea_msigdb_curr_signif_df <- gsea_msigdb_curr_signif_df[order(gsea_msigdb_curr_signif_df$qvalues, decreasing = FALSE), ]
			} else {
				# Else, filter on adjusted p values (which are similar)
				gsea_msigdb_curr_signif_df <- gsea_msigdb_curr[gsea_msigdb_curr@result$p.adjust <= qvalue_threshold, ]
				gsea_msigdb_curr_signif_df <- gsea_msigdb_curr_signif_df[order(gsea_msigdb_curr_signif_df$p.adjust, decreasing = FALSE), ]
			}
			# Combine results from curr df into full df
			gsea_msigdb_df <- rbind(gsea_msigdb_df, gsea_msigdb_curr_signif_df)
			# Combine results from curr database to all others
			gsea_msigdb_list[[database_types[i]]] <- gsea_msigdb_curr
			# Make individual plots for all signif gene sets
			for (j in seq_along(gsea_msigdb_curr_signif_df$ID)) {
				if (gsea_msigdb_curr_signif_df$NES[j] > 0) {
					nes_direction <- "pos"
				} else {
					nes_direction <- "neg"
				}
				if (!dir.exists(paste0("msigdb_gsea/", database_types[i]))) {dir.create(paste0("msigdb_gsea/", database_types[i]), recursive = TRUE)}
				pdf(paste0("msigdb_gsea/", database_types[i], "/", str_pad(j, 4, pad = "0"), "_", nes_direction, "_", gsea_msigdb_curr_signif_df$ID[j], ".pdf"))
				print(
					gseaplot(gsea_msigdb_curr, gsea_msigdb_curr_signif_df$ID[j], title = paste0(plot_title, "\n", gsea_msigdb_curr_signif_df$ID[j], "\n", "NES = ", round(gsea_msigdb_curr_signif_df$NES[j], 4), "\n", "qvalue = ", round(gsea_msigdb_curr_signif_df$qvalues[j], 4)))
				)
				dev.off()
			}
		}
	}
	# Clean up and sort data frame, export to Excel file
	gsea_msigdb_df <- gsea_msigdb_df[order(gsea_msigdb_df$qvalues, decreasing = FALSE), ]
	gsea_msigdb_df$Description <- NULL
	class(gsea_msigdb_df$gene_set_description) <- "hyperlink"
	
	prerank_excel <- as.data.frame(sorted_prerank_df)
	print(head(prerank_excel))
	print(table(prerank_excel$prerank_metric))
	class(prerank_excel$ncbi_gene_db) <- "hyperlink"
	class(prerank_excel$ensembl_gene_db) <- "hyperlink"
	prerank_excel_short <- prerank_excel[, c("gene", "prerank_metric")]
		
		
	excel_list <- list(gsea_msigdb_df, prerank_excel_short, prerank_excel)
	names(excel_list) <- c("gsea_signif", "preranked_gene_order", "preranked_all_columns")
	wb <- openxlsx::write.xlsx(excel_list, file = paste0("msigdb_gsea.xlsx"))
	setColWidths(wb, sheet = 1, cols = 1, widths = 60)
	saveWorkbook(wb, paste0("msigdb_gsea.xlsx"), overwrite = TRUE)
	
	# Create dir for text file output
	if (!dir.exists("suppl_files")) {dir.create("suppl_files")}
	write_tsv(gsea_msigdb_df, "suppl_files/gsea_msigdb.tsv")
	write_tsv(prerank_excel_short, "suppl_files/prerank_less_cols.tsv")
	write_tsv(prerank_excel, "suppl_files/prerank.tsv")
	return(gsea_msigdb_list)
}










# ---------------------------------------------------------------------
# Enrichment testing 
# ---------------------------------------------------------------------


# Wrapper function for the core enricher function
# Takes a de list, runs the hypergeometric test for all, up, and down genes
tk_cluster_profiler <- function(de, log2fc_cutoff = 0.2, padj_cutoff = 0.05, logFC_col = NULL, padj_col = NULL, ensembl_id_col = NULL, entrez_id_col = NULL, msigdb_long, database_types = c("h", "c6"), plot_title = "Plot Title", run_ora = TRUE, run_go = TRUE, run_kegg = TRUE, run_gsea = TRUE) {
	
	# Create dir for text file output
	if (!dir.exists("suppl_files")) {dir.create("suppl_files")}


	# Make the de table nicer -- add links to outside databases
	tk_print_ncbi <- function(entrez_id) {
		if (str_detect(entrez_id, "[[:digit:]]")) {
			x <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", entrez_id)
		} else {
			x <- NA
		}
		return(x)
	}
	de$ncbi_gene_db <- sapply(de[[entrez_id_col]], function(x) {tk_print_ncbi(x)})
	class(de$ncbi_gene_db) <- "hyperlink"
	de$ensembl_gene_db <- sapply(de[[ensembl_id_col]], function(x) {paste0("http://www.ensembl.org/id/", x)})
	class(de$ensembl_gene_db) <- "hyperlink"


	# Get only significant/substantial genes 
	de_signif <- de[abs(de[[logFC_col]]) > log2fc_cutoff & (!is.na(de[[padj_col]]) & de[[padj_col]] < padj_cutoff), ]
	
	
	de_signif_list <- list()
	# Get all entrez ids that are known
	de_signif_list[["all"]] <- de_signif[[entrez_id_col]][str_detect(de_signif[[entrez_id_col]], "[[:digit:]]")]
	de_signif_list[["up"]] <- de_signif[[entrez_id_col]][str_detect(de_signif[[entrez_id_col]], "[[:digit:]]") & de_signif[[logFC_col]] > 0]
	de_signif_list[["dn"]] <- de_signif[[entrez_id_col]][str_detect(de_signif[[entrez_id_col]], "[[:digit:]]") & de_signif[[logFC_col]] < 0]
	

	
	# OVER-REPRESENTATION ANALYSIS
	# For all, up, down lists
	if (run_ora == TRUE) {
		enricher_list <- list()
		enricher_list_df <- list()
		enricher_list_list <- list()
		for (i in seq_along(names(de_signif_list))) {
			if (length(de_signif_list[[i]]) > 2) {
				print(names(de_signif_list)[i])
				enricher_list[[names(de_signif_list)[i]]] <- tk_enricher_core(de_signif_list[[i]], msigdb_long, database_types)
				enricher_list_df[[names(de_signif_list)[i]]] <- enricher_list[[names(de_signif_list)[i]]][[1]]
				enricher_list_list[[names(de_signif_list)[i]]] <- enricher_list[[names(de_signif_list)[i]]][[2]]
			} else {
				print(paste("The length of ", de_signif_list[[i]], " is less than 2."))
				# enricher_list_df[[names(de_signif_list)[i]]] <- data.frame("Not enough genes for analysis")
				# enricher_list_list[[names(de_signif_list)[i]]] <- list("Not enough genes for analysis")
			}
		}
		# Clean up and sort the enrichment results
		for (j in seq_along(names(enricher_list_df))) {
			if (dim(enricher_list_df[[j]])[1] == 0) {
				enricher_list_df[[j]] <- data.frame("Zero pathways enriched.")
			} else {
				enricher_list_df[[j]] <- enricher_list_df[[j]][order(enricher_list_df[[j]]$qvalue, decreasing = FALSE), ]
				enricher_list_df[[j]]$gene_set_description <- msigdb_long$description[match(enricher_list_df[[j]]$ID, msigdb_long$ont)]
				class(enricher_list_df[[j]]$gene_set_description) <- "hyperlink"
			}
		}
		# Write out text files
		for (i in seq_along(enricher_list_df)) {
			write_tsv(enricher_list_df[[i]], path = paste0("suppl_files/msigdb_enrichment_", names(enricher_list_df)[i], ".tsv"))
		}
		# Write out enrichment excel file
		excel_list <- enricher_list_df
		names(excel_list) <- paste0("msigdb_enrichment_", names(enricher_list_df))
		wb <- openxlsx::write.xlsx(excel_list, file = paste0("msigdb_enrichment.xlsx"))
		setColWidths(wb, sheet = 1, cols = 1, widths = 60)
		setColWidths(wb, sheet = 2, cols = 1, widths = 60)
		setColWidths(wb, sheet = 3, cols = 1, widths = 60)
		saveWorkbook(wb, paste0("msigdb_enrichment.xlsx"), overwrite = TRUE)		
		
		# ORA Plots????
	}


	
	
	
	# GO ANALYSIS
	if (run_go == TRUE) {
		print("running GO Group Analysis")
		go_group_list <- list()
		for (i in seq_along(de_signif_list)) {
			go_group_list[[names(de_signif_list)[i]]] <- tk_go_group(geneset = de_signif_list[[i]], go_level = c(4))
		}
		# Write out multi-tab excel
		for (i in seq_along(go_group_list)) {		
			wb <- openxlsx::write.xlsx(go_group_list[[i]], file = paste0("go_groups_", names(go_group_list)[i], ".xlsx"))
			for (j in seq_along(go_group_list[[i]])) {
				setColWidths(wb, sheet = j, cols = 1, widths = 13)
				setColWidths(wb, sheet = j, cols = 2, widths = 40)
				saveWorkbook(wb, paste0("go_groups_", names(go_group_list)[i], ".xlsx"), overwrite = TRUE)
			}
		}
		# Write out text files
		for (i in seq_along(go_group_list)) {
			for (j in seq_along(go_group_list[[i]])) {
				write_tsv(go_group_list[[i]][[j]], path = paste0("suppl_files/go_groups_", names(go_group_list[[i]][j]), "_", names(go_group_list)[i], ".tsv"))
			}
		}
		print("running GO Enrichment Analysis")
		go_enrich_list <- list()
		for (i in seq_along(de_signif_list)) {
			go_enrich_list[[names(de_signif_list)[i]]] <- tk_go_enrich(geneset = de_signif_list[[i]], universe = de[[entrez_id_col]])
		}
		# Write out multi-tab excel
		for (i in seq_along(go_enrich_list)) {		
			wb <- openxlsx::write.xlsx(go_enrich_list[[i]], file = paste0("go_enrich_", names(go_enrich_list)[i], ".xlsx"))
			for (j in seq_along(go_enrich_list[[i]])) {
				setColWidths(wb, sheet = j, cols = 1, widths = 13)
				setColWidths(wb, sheet = j, cols = 2, widths = 40)
				saveWorkbook(wb, paste0("go_enrich_", names(go_enrich_list)[i], ".xlsx"), overwrite = TRUE)
			}
		}
		# Write out text files
		for (i in seq_along(go_enrich_list)) {
			for (j in seq_along(go_enrich_list[[i]])) {
				write_tsv(go_enrich_list[[i]][[j]], path = paste0("suppl_files/go_enrich_", names(go_enrich_list[[i]][j]), "_", names(go_enrich_list)[i], ".tsv"))
			}
		}
	}

	prerank <- tk_create_prerank(de, method = "lfc", logFC_col = logFC_col, padj_col = padj_col, sort = TRUE)


	# KEGG PATHWAYS ANALYSIS
	if (run_kegg == TRUE) {
		print("running KEGG")
		# Run function
		# CHOOSE BOTH UP AND DOWNREGULATED GENES
		kegg_results <- tk_kegg_enrich(geneset = de_signif_list[["all"]], prerank)
		#Excel list
		kegg_enrich_all <- kegg_results[["kegg_enrich_all"]]
		class(kegg_enrich_all$database_link) <- "hyperlink"
		kegg_enrich_signif <- kegg_results[["kegg_enrich_signif"]]
		class(kegg_enrich_signif$database_link) <- "hyperlink"
        kegg_enrich_list <- list(kegg_enrich_signif = kegg_enrich_signif, kegg_enrich_all = kegg_enrich_all)
		# Write out  excel
		wb <- openxlsx::write.xlsx(kegg_enrich_list, file = paste0("kegg_enrich.xlsx"))
		setColWidths(wb, sheet = 1, cols = 1, widths = 11)
		setColWidths(wb, sheet = 1, cols = 2, widths = 40)
		setColWidths(wb, sheet = 2, cols = 1, widths = 11)
		setColWidths(wb, sheet = 2, cols = 2, widths = 40)
		saveWorkbook(wb, paste0("kegg_enrich.xlsx"), overwrite = TRUE)
		# Write out text files
		for (i in seq_along(kegg_enrich_list)) {
			write_tsv(kegg_enrich_list[[i]], path = paste0("suppl_files/", names(kegg_enrich_list)[i], ".tsv"))
		}
	}
	
	
	
	# GSEA
	if (run_gsea == TRUE) {
		print("running GSEA")
		enricher_gsea_output <- tk_enricher_gsea(sorted_prerank_df = prerank, msigdb_long = msigdb_long, plot_title = plot_title, qvalue_threshold = 0.05, database_types = database_types)
		#universe = sorted_prerank_df[[entrez_id_col]]
	}
}














