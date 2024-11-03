# Define the packages to be loaded
pkg = c('GUniFrac', 'Rcpp', 'vegan', 'ggplot2', 'matrixStats', 'Matrix', 
        'ape', 'parallel', 'stats', 'corncob', 'utils', 'statmod', 
        'rmutil', 'dirmult', 'MASS', 'ggrepel', 'foreach', 'modeest', 
        'foreach', 'dplyr', 'tidyverse', 'LinDA', 'metagenomeSeq', 
        'microeco', 'magrittr', 'edgeR', 'ANCOMBC', 'fastANCOM', 
        'dacomp', 'curatedMetagenomicData', 'phyloseq', 'limma', 
        'DEFormats', 'apeglm', 'Maaslin2', 'DESeq2', 'mia', 
        'ALDEx2', 'dplyr')  

# Load the specified packages and suppress startup messages
# sapply is used to apply the require function to each package;
# require loads the R package functions and returns TRUE if the package is already loaded
suppressPackageStartupMessages(sapply(pkg, require, character = TRUE))

#####################################
# The input 'data' is a list that should contain:
# - data: the dataset
# - group: the grouping variable
# - gold_name: the true differential components
#####################################
####Wilcox####
res_wilcox <- function(data=sim,p.adj="fdr",p.cutoff=0.05){
   p_wilcox <- apply(data$data,1,function(x) wilcox.test(x ~ group,
      data=cbind(x,data[["group"]]))$p.value)
   p_wilcox[is.nan(p_wilcox)] <- 1
   
   adjusted_p_wilcox <- p.adjust(p_wilcox, method = p.adj)
   gold_names_wilcox <- names(adjusted_p_wilcox)[adjusted_p_wilcox <= p.cutoff]
   
   return(list(p_value = adjusted_p_wilcox, gold_names=gold_names_wilcox))
}

#####################################
####LinDA####
res_linda <- function(data=sim,p.adj = "BH",p.cutoff=0.05){
   linda.obj <- LinDA::linda(data$data, data$group , formula = '~group', alpha = 0.05,
      prev.cut = 0, lib.cut = 1000, winsor.quan = 0.97,p.adj.method = p.adj)$output[[1]]
   
   adjusted_p_linda <- linda.obj$padj
   gold_names_linda <- rownames(linda.obj)[adjusted_p_linda <= p.cutoff]
   
   return(list(p_value = adjusted_p_linda, gold_names=gold_names_linda))
}

#####################################
####metagenomeSeq##### 
res_metagenomeseq <- function(data=sim,p.adj = "fdr",p.cutoff=0.05){
  
  metadata = AnnotatedDataFrame(data$group)
  taxa1 <- as.data.frame(data$data[,c(1,3)])
  taxa <- AnnotatedDataFrame(taxa1)
  otu <- data.frame(data$data)
  obj = newMRexperiment(otu,phenoData=metadata,featureData=taxa)
  
  p = cumNormStatFast(obj)  #CSS归一化
  obj <- cumNorm(obj,p = p)
  pd <- pData(obj)
  mod <- model.matrix(~1+group, data=pd)
  mSeq_res = fitFeatureModel(obj,mod)  
  
  adjusted_p_metagenomeseq <- as.data.frame(p.adjust(mSeq_res@pvalues,method = p.adj))
  
  gold_names_metagenomeseq <- na.omit(rownames(adjusted_p_metagenomeseq)[adjusted_p_metagenomeseq <= p.cutoff])
  return(list(p_value = adjusted_p_metagenomeseq, gold_names=gold_names_metagenomeseq))
}
#####################################
####lefse####
res_lefse <- function(data=sim,p.adj = "fdr",p.cutoff=0.05){
   dataset <- microtable$new(sample_table = data$group,
      otu_table = as.data.frame(data$data), 
      tax_table = tax_tab)
   lefse <- trans_diff$new(dataset = dataset, 
      method = "lefse", 
      group = "group", 
      p_adjust_method = p.adj,
      alpha = 1, 
      lefse_subgroup = NULL,
      taxa_level = "Kingdom")
   
   adjusted_p_lefse <- p.adjust(lefse$res_diff$P.unadj, method = p.adj)
   gold_names_lefse <- rownames(lefse$res_diff)[adjusted_p_lefse <= p.cutoff]
   
   return(list(p_value = adjusted_p_lefse, gold_names=gold_names_lefse))
}

#####################################
####edgeR####
res_edger <- function(data=sim,p.adj = "fdr",p.cutoff=0.05){
  group=as.factor(as.matrix(data$group))
  design <- model.matrix(~group) 
  
  dge <- DGEList(counts=data$data, group=group)
  keep.exprs <- filterByExpr(dge) 
  dge <- dge[keep.exprs,,keep.lib.sizes=FALSE] 
  dge <- calcNormFactors(dge, method = "TMM") 
  dge <- estimateDisp(dge, design, robust=T) 
  
  fit <- glmFit(dge, design, robust=T)  
  lrt <- glmLRT(fit, coef=2)
  tempDEG <- topTags(lrt, n = Inf)$table
  
  adjusted_p_edger <- p.adjust(tempDEG$PValue,method = p.adj)
  gold_names_edger <- row.names(tempDEG)[adjusted_p_edger <= p.cutoff]
  
  return(list(p_value = adjusted_p_edger, gold_names=gold_names_edger))
}
#####################################
####ALDEx2####
## type = "we" indicates a t-test, while type = "wi" indicates a Wilcoxon rank-sum test
res_aldex2 <- function(data=sim,p.adj = "BH",p.cutoff=0.05, type="we" ){
   x <- aldex.clr(data$data,data$group$group, mc.samples=128, verbose=FALSE)
   x.tt <- aldex.ttest(x, hist.plot = F, paired.test=FALSE,verbose = F)
   p_aldex2 = switch(type,
      "we" =  x.tt$we.ep ,
      "wi" =  x.tt$wi.ep )
   
   adjusted_p_aldex2 <- p.adjust(p_aldex2, method=p.adj)
   gold_names_aldex2 <- (rownames(x.tt))[adjusted_p_aldex2 <= p.cutoff]
   
   return(list(p_value=adjusted_p_aldex2, gold_names=gold_names_aldex2))
}

#####################################
####fastANCOM####
res_fastancom <- function(data=sim,p.adj = "BH",p.cutoff=0.05){
   fit <- fastANCOM(Y=as.matrix(t(data$data)), 
      x=as.matrix(data$group), 
      sig = p.cutoff, pseudo = 0.5)
   
   adjusted_p_fastancom <-p.adjust(fit$results$final$log2FC.pval,method = p.adj)
   gold_names_fastancom <- rownames(fit$results$final)[adjusted_p_fastancom <= p.cutoff]
   
   return(list(p_value = adjusted_p_fastancom, gold_names=gold_names_fastancom))
}

#####################################
####DACOMP####
res_dacomp <- function(data=sim,p.adj = "BH",p.cutoff=0.05){
   ##Choosing the minimal number of reference taxa so that at least 50 reads are available under the reference for all samples
   result.selected.references = dacomp.select_references(X = as.matrix(t(data$data)), 
      minimal_TA = 50, verbose = T)
   dacomp = dacomp.test(X = as.matrix(t(data$data)),
      y = as.vector(as.matrix(data$group)),
      ind_reference_taxa = result.selected.references, 
      test = DACOMP.TEST.NAME.WILCOXON,nr_perm = 1000,
      verbose = FALSE,disable_DSFDR=TRUE)
   
   adjusted_p_dacomp <- p.adjust(dacomp$p.values.test,method = p.adj)
   adjusted_p_dacomp[is.na(adjusted_p_dacomp)] <- 1
   
   gold_names_dacomp <- rownames(data$data)[adjusted_p_dacomp <= p.cutoff]
   
   return(list(p_value = adjusted_p_dacomp, gold_names=gold_names_dacomp))
}

#####################################
####ZicoSeq####
res_zicoseq <- function(data=sim,p.adj = "BH",p.cutoff=0.05){
   data$data[apply(data$data,1,sd)==0,1] <- 1
   zico.obj <- ZicoSeq(meta.dat = as.data.frame(data$group), feature.dat = data$data,
      grp.name = 'group', feature.dat.type = "count",
      # Filter to remove rare taxa
      prev.filter = 0.2, mean.abund.filter = 0, max.abund.filter = 0.002, min.prop = 0,
      # Winsorization to replace outliers
      is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
      # Posterior sampling to impute zeros
      is.post.sample = TRUE, post.sample.no = 25,
      # Multiple link functions to capture diverse taxon-covariate relation
      link.func = list(function (x) x^0.25, function (x) x^0.5, function (x) x^0.75),
      stats.combine.func = max,
      # Permutation-based multiple testing correction
      perm.no = 99, strata = NULL,
      # Reference-based multiple stage normalization
      ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
      # Family-wise error rate control
      is.fwer = FALSE,
      verbose = TRUE, return.feature.dat = FALSE)
   
   adjusted_p_zicoseq <- p.adjust(zico.obj$p.raw,method = p.adj) 
   gold_names_zicoseq <- names(adjusted_p_zicoseq)[adjusted_p_zicoseq <= p.cutoff]
   
   return(list(p_value = adjusted_p_zicoseq, gold_names=gold_names_zicoseq))
}

#####################################
####limma-voom####
res_limma_voom <- function(data = sim, p.adj = "BH", p.cutoff = 0.05) {
  physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE),     # Create a phyloseq object with the OTU table and sample data
                     sample_data(data.frame(sim$group)))
  # Convert phyloseq object to DESeq2 format
  dds <- phyloseq_to_deseq2(physeq, ~ group)  # Convert to DESeq2 and DGEList objects
  dge <- as.DGEList(dds)  # Convert to DGEList format
  dge <- calcNormFactors(dge, method = "TMM")  # Calculate TMM normalization factors
  mm <- model.matrix(~ group, dge$samples)  # Create model matrix for group
  y <- voom(dge, mm)  # Obtain voom weights for the data
  fit <- lmFit(y, mm)  # Fit linear model with limma
  # Perform empirical Bayes moderation on the linear model statistics
  fit <- eBayes(fit)
  
  adjusted_p_limma <- p.adjust(fit$p.value[, 2], method = p.adj)  
  gold_names_limma <- names(adjusted_p_limma)[adjusted_p_limma <= p.cutoff]
  
  return(list(p_value = adjusted_p_limma, gold_names = gold_names_limma))
}

#####################################
####DESeq2####
res_deseq2 <- function(data=sim,p.adj = "BH",p.cutoff=0.05){
      physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
      sample_data(data.frame(sim$group)))
   dds <- phyloseq_to_deseq2(physeq, ~ group)      #convert to DESeq2 and DGEList objects
   dds1 <- DESeq(dds, test = "Wald", fitType = "local", sfType = "poscounts")
   res <- lfcShrink(dds1, coef=2, type="apeglm")    
   
   adjusted_p_deseq2 <- p.adjust(res$pvalue,method = p.adj)
   gold_names_deseq2 <- rownames(res)[adjusted_p_deseq2 <= p.cutoff]
   
   return(list(p_value = adjusted_p_deseq2, gold_names=gold_names_deseq2))
}

#####################################
####Corncob####
res_corncob <- function(data=sim,p.adj = "fdr",p.cutoff=0.05){
   physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
      sample_data(data.frame(sim$group)))
   corn_da <- differentialTest(formula = ~ group,
      phi.formula = ~ 1,
      formula_null = ~ 1,
      phi.formula_null = ~ 1,
      data = physeq,
      test = "Wald", boot = FALSE,
      fdr_cutoff = 0.05)
   
   adjusted_p_corncob <- na.omit(p.adjust(corn_da$p,method = p.adj))
   gold_names_corncob <- names(adjusted_p_corncob)[adjusted_p_corncob <= p.cutoff]
   
   return(list(p_value = adjusted_p_corncob, gold_names=gold_names_corncob))
}
#####################################
####maaslin2####
res_maaslin2 <- function(data=sim,p.adj = "fdr",p.cutoff=0.05){
   physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
      sample_data(data.frame(sim$group)))
   res_mas <- Maaslin2(input_data = data.frame(otu_table(physeq)),
      input_metadata = data.frame(sample_data(physeq)),
      output = "./Maaslin2_default_output",
      min_abundance = 0.0,
      min_prevalence = 0.0,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      max_significance = 0.05,
      fixed_effects = "group",
      correction = "BH",
      standardize = FALSE,
      cores = 1)
   
   adjusted_p_maaslin2 <- p.adjust(res_mas$results$pval,method = p.adj) 
   gold_names_maaslin2 <- str_replace(res_mas$results$feature,"X","")[adjusted_p_maaslin2 <= p.cutoff]
   
   return(list(p_value = adjusted_p_maaslin2, gold_names=gold_names_maaslin2))
}




