library(data.table)
library(edgeR)
library(limma)
library(stringr)
library(Biobase)
library(biomaRt)
library(tools)
library(ggplot2)
library(Hmisc)
library(ggfortify)
library(ggrepel)

# Function for normalizing counts, making expressionset and DE analysis

# function for normalizing counts
count_norm <- function(raw_counts, gene_id_column, perc_low = 0.3){
  
  # identify lowly expressed genes (expression of a gene is less than 1 in % of the samples)
  gene_low <- rowSums(raw_counts <= 1) >= ceiling(perc_low * ncol(raw_counts))
  
  # remove lowly expressed genes
  if(is.null(perc_low)){raw_counts_low_rm <- raw_counts}else{raw_counts_low_rm <- raw_counts[!gene_low, ]}
  
  # create a DGElist from the the raw counts data frame
  dge <- DGEList(counts=raw_counts_low_rm)
  
  # Calculate normalization factors to scale the raw library sizes
  dge <- calcNormFactors(dge, method = "TMM")
  
  # normalize the raw counts
  v <- voom(dge, plot=TRUE)
  
  # show mean-variance plot
  v
  
  return(as.data.frame(v$E))
  
}

count_norm_dgelist <- function(raw_counts, gene_id_column, meta, perc_low = 0.3){
  
  # identify lowly expressed genes (expression of a gene is less than 1 in % of the samples)
  gene_low <- rowSums(raw_counts <= 1) >= ceiling(perc_low * ncol(raw_counts))
  
  # remove lowly expressed genes
  if(is.null(perc_low)){raw_counts_low_rm <- raw_counts}else{raw_counts_low_rm <- raw_counts[!gene_low, ]}
  
  # make annotation table
  if (org == "mouse"){
		ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
		annot <- getBM(c("ensembl_gene_id", "mgi_symbol"), mart=ensembl)
		# make gene id - gene symbol annotation table
		annot <- data.frame(Gene_ID = rownames(raw_counts_low_rm), Gene_Symbol = annot$mgi_symbol[match(rownames(raw_counts_low_rm), annot$ensembl_gene_id)])
		rownames(annot) <- rownames(raw_counts_low_rm)
		}else if (org == "human") {
			ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
			annot <- getBM(c("ensembl_gene_id", "hgnc_symbol"), mart=ensembl)
			# make gene id - gene symbol annotation table
			annot <- data.frame(Gene_ID = rownames(raw_counts_low_rm), Gene_Symbol = annot$hgnc_symbol[match(rownames(raw_counts_low_rm), annot$ensembl_gene_id)])
			rownames(annot) <- rownames(raw_counts_low_rm)}
  
  # create a DGElist from the the raw counts data frame
  
  
  if(identical(rownames(raw_counts_low_rm), rownames(annot)) & identical(colnames(raw_counts_low_rm), rownames(meta))){dge <- DGEList(counts=raw_counts_low_rm, genes = annot, samples = meta)}else{
	raw_counts_low_rm <- raw_counts_low_rm[, rownames(meta)]
	raw_counts_low_rm <- raw_counts_low_rm[rownames(annot), ]
	dge <- DGEList(counts=raw_counts_low_rm, genes = annot, samples = meta)}
  
  # Calculate normalization factors to scale the raw library sizes
  dge <- calcNormFactors(dge, method = "TMM")
  
  # normalize the raw counts
  v <- voom(dge, plot=TRUE)
  
  # show mean-variance plot
  show(v)
  
  return(v)
  
}

# Function for making expressionset
make_exprset <- function(expression_matrix, meta_data, sample_id_column, org = "mouse"){
  
  # convert expression_matrix to a matrix
  expression_matrix <- as.matrix(expression_matrix)
  
  # make annotation table
  if (org == "mouse"){
		ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
		annot <- getBM(c("ensembl_gene_id", "mgi_symbol"), mart=ensembl)
		# make gene id - gene symbol annotation table
		annot <- data.frame(Gene_ID = rownames(expression_matrix), Gene_Symbol = annot$mgi_symbol[match(rownames(expression_matrix), annot$ensembl_gene_id)])
		rownames(annot) <- rownames(expression_matrix)
		}else if (org == "human") {
			ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
			annot <- getBM(c("ensembl_gene_id", "hgnc_symbol"), mart=ensembl)
			# make gene id - gene symbol annotation table
			annot <- data.frame(Gene_ID = rownames(expression_matrix), Gene_Symbol = annot$hgnc_symbol[match(rownames(expression_matrix), annot$ensembl_gene_id)])
			rownames(annot) <- rownames(expression_matrix)}
  
  
  
  # checking if rownames and colnames of meta and expr are same
  if(identical(colnames(expression_matrix), rownames(meta_data))){
    # creating meta info for the expressionset
    phenoData <- new("AnnotatedDataFrame", data = meta_data)
    # creating feature info for the expressionset
    featureData <- new("AnnotatedDataFrame", data = annot)
    # constructing an expressionset
    exprset <- ExpressionSet(assayData = expression_matrix, phenoData = phenoData,featureData = featureData)
	}else{
		# match sample order between expr and meta
		expression_matrix <- expression_matrix[, rownames(meta_data)]
		# creating meta info for the expressionset
		phenoData <- new("AnnotatedDataFrame", data = meta_data)
		# creating feature info for the expressionset
		featureData <- new("AnnotatedDataFrame", data = annot)
		exprset <- ExpressionSet(assayData = expression_matrix, phenoData = phenoData,featureData = featureData)}
}

# Function for DE analysis
de_analysis_exprset_subset <- function(expressionset, contrast, treatment_column, p_value = 1, fc = 1, cov_vec = NULL){ 
  # subset expressionset
  expressionset <- expressionset[, expressionset[[treatment_column]]  == str_split(contrast[[1]], "-")[[1]][1] | expressionset[[treatment_column]]  == str_split(contrast[[1]], "-")[[1]][2]]
  
  print("Sample Size: ")
  print(table(pData(expressionset)[[treatment_column]]))
  
  # make design and contrast matrix
  if(!is.null(cov_vec)){
    design <- model.matrix(as.formula(paste(c("~0", "expressionset[[treatment_column]]", paste0("expressionset", "[[", "'", cov_vec, "'", "]]")), collapse=" + ")))
    colnames(design) <- c(c(str_split(colnames(design)[1:2], "]]")[[1]][2], str_split(colnames(design)[1:2], "]]")[[2]][2]) , paste0("cov_", seq(ncol(design))[-(1:2)]))
    ctrst <- makeContrasts(contrasts = contrast, levels = colnames(design))
  }else{
		design <- model.matrix(~ 0 + factor(expressionset[[treatment_column]]))
		colnames(design) <- c(str_split(colnames(design)[1:2], "]]\\)")[[1]][2], str_split(colnames(design)[1:2], "]]\\)")[[2]][2])
		ctrst <- makeContrasts(contrasts = contrast, levels = design)
  }
  
  
  # de analyis
  fit <- lmFit(expressionset, design)
  fit_2 <- contrasts.fit(fit, ctrst)
  fit_2 <- eBayes(fit_2)
  
  # saving the results
  res <- topTable(fit_2, n = Inf, p.value=p_value, lfc=log2(fc))
  
  # print the number of DEGs identified
  print(paste(nrow(res), "DEGs were detected"))
  
  return(res)
}


de_analysis_dgelist_subset <- function(exprlist, contrast, treatment_column, p_value = 1, fc = 1, cov_vec = NULL){ 
  # subset exprlist
  exprlist <- exprlist[, exprlist$targets[[treatment_column]]  == str_split(contrast[[1]], "-")[[1]][1] | exprlist$targets[[treatment_column]]  == str_split(contrast[[1]], "-")[[1]][2]]
  
  print("Sample Size: ")
  print(table(exprlist$targets[[treatment_column]]))
  
  # make design and contrast matrix
  if(!is.null(cov_vec)){
    design <- model.matrix(as.formula(paste(c("~0", "exprlist$targets[[treatment_column]]", paste0("exprlist$targets", "[[", "'", cov_vec, "'", "]]")), collapse=" + ")))
    colnames(design) <- c(c(str_split(colnames(design)[1:2], "]]")[[1]][2], str_split(colnames(design)[1:2], "]]")[[2]][2]) , paste0("cov_", seq(ncol(design))[-(1:2)]))
    ctrst <- makeContrasts(contrasts = contrast, levels = colnames(design))
  }else{
		design <- model.matrix(~ 0 + factor(exprlist$targets[[treatment_column]]))
		colnames(design) <- c(str_split(colnames(design)[1:2], "]]\\)")[[1]][2], str_split(colnames(design)[1:2], "]]\\)")[[2]][2])
		ctrst <- makeContrasts(contrasts = contrast, levels = design)
  }
  
  
  # de analyis
  fit <- lmFit(exprlist, design)
  fit_2 <- contrasts.fit(fit, ctrst)
  fit_2 <- eBayes(fit_2)
  
  # saving the results
  res <- topTable(fit_2, n = Inf, p.value=p_value, lfc=log2(fc))
  
  # print the number of DEGs identified
  print(paste(nrow(res), "DEGs were detected"))
  
  return(res)
}


de_analysis_dgelist_full <- function(exprlist, treatment_column, p_value = 1, fc = 1, cov_vec = NULL){ 
  
  # make design and contrast matrix
  if(!is.null(cov_vec)){
    design <- model.matrix(as.formula(paste(c("~0", "exprlist$targets[[treatment_column]]", paste0("exprlist$targets", "[[", "'", cov_vec, "'", "]]")), collapse=" + ")))
	
    colnames(design)[1:length(unique(exprlist$targets[[treatment_column]]))] <- gsub(".*]]", "", colnames(design)[1:length(unique(exprlist$targets[[treatment_column]]))])
	
	colnames(design)[(length(unique(exprlist$targets[[treatment_column]]))+1):ncol(design)] <- paste0("Cov", (length(unique(exprlist$targets[[treatment_column]]))+1):ncol(design))
	
	contrast_pairwise <- apply(combn(unique(exprlist$targets[[treatment_column]]), 2), 2, paste, collapse = "-")
	
	eval(parse(text = paste0("ctrst <- makeContrasts(", str_c(c(unname(sapply(contrast_pairwise, function(a) paste(paste0("'", a, "'"), "=", paste0("'", a, "'")))), "levels = design"), collapse=','), ")")))
	
  }else{
		design <- model.matrix(~ 0 + exprlist$targets[[treatment_column]])
		
		colnames(design) <- gsub(".*]]", "", colnames(design))
		
		contrast_pairwise <- apply(combn(unique(exprlist$targets[[treatment_column]]), 2), 2, paste, collapse = "-")
	
		eval(parse(text = paste0("ctrst <- makeContrasts(", str_c(c(unname(sapply(contrast_pairwise, function(a) paste(paste0("'", a, "'"), "=", paste0("'", a, "'")))), "levels = design"), collapse=','), ")")))
  }
  
  
  # de analyis
  fit <- lmFit(exprlist, design)
  fit_2 <- contrasts.fit(fit, ctrst)
  fit_2 <- eBayes(fit_2)
  
  res <- lapply(colnames(ctrst), function(x) {topTable(fit_2, coef = x, n = Inf, p.value=p_value, lfc=log2(fc), sort.by = "p")})
  
  # remove comparisons with 0 DEGs
  deg_non_0 <- sapply(res, function(x) nrow(x) != 0)
  
  res_deg0_rm <- res[deg_non_0]
  
  ctrst_deg0_rm <- colnames(ctrst)[deg_non_0]
  
  res_deg0_rm <- mapply(function(x, y) {x$Contrast <- rep(y, nrow(x)); return(x)}, res_deg0_rm, ctrst_deg0_rm, SIMPLIFY = F)
  
  res_deg0_rm <- rbindlist(res_deg0_rm)

  
  # print the number of DEGs identified
  print(sapply(unique(res_deg0_rm$Contrast), function(x) nrow(res_deg0_rm[Contrast == x])))
  
  return(res_deg0_rm)
}


# Plot Functional Enrichment Results
plot_functional_enrich_2group <- function(enrich_file, Contrast, save_dir, out_format = "png", dpi = 300, width = 25, height = 15){

	if(length(enrich_file) == 1){
    # read deg functional enrichment results
    deg_func <- fread(enrich_file)
    
    # format deg_func
    deg_func <- deg_func[!System %in% c("TranscriptionFactor", "MIR")]
    
    deg_func$neglogP <- -log10(deg_func$Corrected_P)
    
    deg_func_up <- deg_func[grepl("up", deg_func$Module)][1:15, ]
    deg_func_down <- deg_func[grepl("down", deg_func$Module)][1:15, ]
    
    deg_func_down$neglogP <- -1 * deg_func_down$neglogP
    
    deg_func_top <- rbind(deg_func_up, deg_func_down)
    
    deg_func_top <- deg_func_top[order(-deg_func_top$neglogP), ][!duplicated(deg_func_top$Gene.Category)]
    
    # plot top functions
    plot_deg_func_top <- ggplot(deg_func_top, aes(x=reorder(`Gene.Category`, neglogP), y=neglogP, label=neglogP)) + geom_bar(stat='identity', aes(fill=Module), width=.5) + scale_fill_manual(values = c("#00ba38", "#f8766d")) + coord_flip() + theme_bw() + theme(legend.position = "none", plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) + labs(x = "Functions", y = "-log10 P Value")
	
	# save plot
	ggsave(paste0(list.files(save_dir, "Ontology"), ".", out_format), plot_deg_func_top, device = out_format, dpi = dpi, width = width, height = height, unit = "cm")
	
}else{print("Enrichment file does not exist! Plotting terminated!")}
}


plot_functional_enrich_multi <- function(enrich_file, out_format = "png", dpi = 300, width = 25, height = 15){

	if(length(enrich_file) == 1){
    # read deg functional enrichment results
    deg_func <- fread(enrich_file)
    
    # format deg_func
    deg_func <- deg_func[!System %in% c("TranscriptionFactor", "MIR")]
    
    deg_func$neglogP <- -log10(deg_func$Corrected_P)
    
    deg_func_sublist <- lapply(unique(gsub("_up|_down", "", deg_func[[1]])), function(x) deg_func[grepl(x, deg_func$Module)])
	names(deg_func_sublist) <- unique(gsub("_up|_down", "", deg_func[[1]]))
	
	deg_func_up <- lapply(deg_func_sublist, function(x) x[grepl("up", x$Module)][1:15, ])
	deg_func_down <- lapply(deg_func_sublist, function(x) x[grepl("down", x$Module)][1:15, ])
	
	deg_func_down <- lapply(deg_func_down, function(x) {x$neglogP <- -1 * x$neglogP; return(x)})
	
	deg_func_top <- mapply(rbind, deg_func_up, deg_func_down, SIMPLIFY  = F)
	
	deg_func_top <- lapply(deg_func_top, function(x) x[order(-x$neglogP), ][!duplicated(x$Gene.Category)])
	
	plot <- lapply(deg_func_top, function(x) ggplot(x, aes(x=reorder(`Gene.Category`, neglogP), y=neglogP, label=neglogP)) + geom_bar(stat='identity', aes(fill=Module), width=.5) + scale_fill_manual(values = c("#00ba38", "#f8766d")) + coord_flip() + theme_bw() + theme(legend.position = "none", plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) + labs(x = "Functions", y = "-log10 P Value"))
	
	mapply(ggsave, paste0("Plot_DEG_Ontology_", names(deg_func_sublist), ".", out_format), plot, MoreArgs = list(device = out_format, dpi = dpi, width = width, height = height, unit = "cm")) 
	

	
}else{print("Enrichment file does not exist! Plotting terminated!")}
}


# Plot gene expression
plot_gene_exprset <- function(exprset, gene_symbol, treatment_col, save_dir, out_format, dpi = 300, width = 25, height =30, plot_type = "bar"){
		
	# subset expressionset
	geneid <- fData(exprset)[[1]][match(gene_symbol, fData(exprset)[[2]])][1]
	
	if(!is.na(geneid)){
	exprset_sub <- exprset[geneid, ]
	
	if(identical(rownames(t(exprs(exprset_sub))), rownames(pData(exprset_sub)))){exprset_sub_melt <- cbind(pData(exprset_sub), t(exprs(exprset_sub)))}
	
	if(plot_type == "violin"){
		# violin plot
		violin_plot <- ggplot(exprset_sub_melt, aes_string(x=treatment_col, y=geneid)) + geom_violin(trim=FALSE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + stat_summary(geom="pointrange", fun.data = "median_hilow", color="red", size = 1.2) + xlab("") + ylab("Expression") + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text=element_text(size=25, face = "bold"), axis.text.x = element_text(angle = 90), axis.title=element_text(size=25,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
	
		# save plot
		ggsave(paste0(save_dir, "/", gene_symbol, "_violin.", out_format), violin_plot, device = out_format, dpi = dpi, width = width, height = height, unit = "cm")
	}else if(plot_type == "bar"){
			# aggregate data
			exprset_sub_melt_bar <- data.frame(treatment = unique(exprset_sub_melt[[treatment_col]]), mean = sapply(unique(exprset_sub_melt[[treatment_col]]), function(x) mean(exprset_sub_melt[exprset_sub_melt[[treatment_col]] == x, ][[geneid]])), se = sapply(unique(exprset_sub_melt[[treatment_col]]), function(x) sd(exprset_sub_melt[exprset_sub_melt[[treatment_col]] == x, ][[geneid]])/sqrt(length(x))))
			
			# bar plot
			bar_plot <- ggplot(exprset_sub_melt_bar, aes_string(x="treatment", y="mean")) + geom_bar(stat="identity", width = 0.5, fill="tomato2") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9)) + xlab("") + ylab("Expression") + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text=element_text(size=25, face = "bold"), axis.text.x = element_text(angle = 90), axis.title=element_text(size=25,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
	
			# save plot
			ggsave(paste0(save_dir, "/", gene_symbol, "_bar.", out_format), bar_plot, device = out_format, dpi = dpi, width = width, height = height, unit = "cm")
	}}else{print(paste(gene_symbol, "not found in the dataset!"))}
}


plot_gene_dgelist <- function(exprlist, gene_symbol, treatment_col, save_dir, out_format, dpi = 300, width = 25, height =30, plot_type = "violin"){
		
	# subset expressionset
	geneid <- exprlist$genes[[1]][match(gene_symbol, exprlist$genes[[2]])][1]
	
	if(!is.na(geneid)){
	exprlist_sub <- exprlist[geneid, ]
	
	if(identical(rownames(t(exprlist_sub$E)), rownames(exprlist_sub$targets))){exprlist_sub_melt <- cbind(exprlist_sub$targets, t(exprlist_sub$E))}
	
	if(plot_type == "violin"){
		# violin plot
		violin_plot <- ggplot(exprlist_sub_melt, aes_string(x=treatment_col, y=geneid)) + geom_violin(trim=FALSE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + stat_summary(geom="pointrange", fun.data = "median_hilow", color="red", size = 1.2) + xlab("") + ylab("Expression") + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text=element_text(size=25, face = "bold"), axis.text.x = element_text(angle = 90), axis.title=element_text(size=25,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
	
		# save plot
		ggsave(paste0(save_dir, "/", gene_symbol, "_violin.", out_format), violin_plot, device = out_format, dpi = dpi, width = width, height = height, unit = "cm")
	}else if(plot_type == "bar"){
			# aggregate data
			exprlist_sub_melt_bar <- data.frame(treatment = unique(exprlist_sub_melt[[treatment_col]]), mean = sapply(unique(exprlist_sub_melt[[treatment_col]]), function(x) mean(exprlist_sub_melt[exprlist_sub_melt[[treatment_col]] == x, ][[geneid]])), se = sapply(unique(exprlist_sub_melt[[treatment_col]]), function(x) sd(exprlist_sub_melt[exprlist_sub_melt[[treatment_col]] == x, ][[geneid]])/sqrt(length(x))))
			
			# bar plot
			bar_plot <- ggplot(exprlist_sub_melt_bar, aes_string(x="treatment", y="mean")) + geom_bar(stat="identity", width = 0.5, fill="tomato2") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9)) + xlab("") + ylab("Expression") + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text=element_text(size=25, face = "bold"), axis.text.x = element_text(angle = 90), axis.title=element_text(size=25,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
	
			# save plot
			ggsave(paste0(save_dir, "/", gene_symbol, "_bar.", out_format), bar_plot, device = out_format, dpi = dpi, width = width, height = height, unit = "cm")
	}}else{print(paste(gene_symbol, "not found in the dataset!"))}
}

# plot sample PCA
plot_pca <- function(count_table, meta_table, shape_column, color_column, sample_label = TRUE, out_format, dpi = 300, width = 24, height =18, save_dir){

	count_table_t <- as.data.frame(t(count_table))
	
	count_table_t <- count_table_t[rownames(meta_table), ]
	
	if(identical(rownames(count_table_t), rownames(meta_table))){
		pca_res <- prcomp(count_table_t, scale. = TRUE)
		if(sample_label){
		pca_plot <- autoplot(pca_res, label = F, shape = shape_column, data = meta_table, colour = color_column) + geom_text_repel(aes(label = rownames(meta_table))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())}else{pca_plot <- autoplot(pca_res, label = F, shape = shape_column, data = meta_table, colour = color_column) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())}
    }else{print("Sample order does not match between count table and meta table")}
	
	show(pca_plot)
	
	# save plot
	ggsave(paste0(save_dir, "/", "Sample_PCA.", out_format), pca_plot, device = out_format, dpi = dpi, width = width, height = height, unit = "cm")
}