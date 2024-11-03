#### OTU = 100, 200 #####
setwd("/public/linzhf/simulation_13DAA")  # Set working directory
source('all_function.R')  # Load custom functions
source('simulatemseq_new.R')  # Load simulation functions

## Simulation parameters setup   
notu <- c(100, 200)  # Number of OTUs to simulate

# Create a save path based on the number of OTUs
save_path <- paste0("/public/linzhf/simulation_13DAA/simulation_", notu, "/")
if (dir.exists(save_path)) file.remove(save_path)  # Remove existing directory
if (!dir.exists(save_path)) dir.create(save_path)  # Create new directory

# Sample sizes
nsample <- c(20, 50, 100, 200, 300, 500)

# Types of OTUs
otutype <- c("abundant", "rare", "mix")

# Effect sizes for simulation
effectsize <- c(1, 2, 5, 10)

# Number of simulation iterations
B <- 100

# Load throat OTU table data
data(throat.otu.tab)  # Load the data
comm <- t(throat.otu.tab)  # Transpose the data
comm <- comm[rowMeans(comm != 0) > 0.05, ]  # Filter out low-abundance OTUs

# Simulate binary covariate, 5% signal density, abundant differential OTUs, unbalanced change
# This setting simulates strong compositional effects
p.adj <- "fdr"  # Method for p-value adjustment
p.cutoff <- 0.05  # Significance threshold

set.seed(123)
for(i in nsample) {
  for(j in notu){
    for(k in otutype){
      for(l in effectsize){
        res_sim <- list()
        for(p in 1:B){
          sim.obj <- SimulateMSeq_new(
            ref.otu.tab = comm, 
            nSam = i, 
            nOTU = j,
            # True signal setting
            diff.otu.pct = 0.1, 
            diff.otu.direct = c("unbalanced"),
            diff.otu.mode = k,
            covariate.type = c("binary"), 
            grp.ratio = 1,
            covariate.eff.mean = log(l), 
            covariate.eff.sd = 0,
            # Confounder signal setting
            confounder.type = c("none"), 
            conf.cov.cor = 0.6,
            conf.diff.otu.pct = 0.1, 
            conf.nondiff.otu.pct = 0.1,
            confounder.eff.mean = 1.0, 
            confounder.eff.sd = 0,
            # Depth setting
            depth.mu = 10000, 
            depth.theta = 5, 
            depth.conf.factor = 0
          )
          #####################################
          sim <- list(data=sim.obj[["otu.tab.sim"]],
                      group = sim.obj[["covariate"]] %>% as.data.frame() %>% setNames("group"),
                      gold_name = sim.obj$otu.names[sim.obj$diff.otu.ind]
          )
          
          colnames(sim$data) <- paste0("sample",1:dim(sim$data)[2])
          rownames(sim$group)<- paste0("sample",1:dim(sim$data)[2])
          
          physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
                             sample_data(data.frame(sim$group)))
          
          tax_tab <- matrix(1,nrow = nrow(sim$data), ncol = 7)
          rownames(tax_tab) = rownames(sim$data)
          colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species")
          tax_tab <- as.data.frame(tax_tab)
          tax_tab$Kingdom <- rownames(sim$data)
          
          rl_name <- paste(i,j,k,l,p,sep = "_")
          
          tryCatch({
            start_time <- Sys.time() 
            res_sim[paste(rl_name,"wilcox",sep="_")] <- list(res_wilcox(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"wilcox",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error wilcox in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"linda",sep="_")] <- list(res_linda(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"linda",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error linda in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"metagenomeseq",sep="_")] <- list(res_metagenomeseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"metagenomeseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error metagenomeseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"lefse",sep="_")] <- list(res_lefse(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"lefse",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error lefse in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"rf",sep="_")] <- list(res_rf(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"rf",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error rf in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"edger",sep="_")] <- list(res_edger(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"edger",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error edger in ",rl_name))} )
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"ancombc2",sep="_")] <- list(res_ancombc2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"ancombc2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error ancombc2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_we",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_we",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_we in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_wi",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff, type="wi"))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_wi",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_wi in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"fastancom",sep="_")] <- list(res_fastancom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"fastancom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error fastancom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"dacomp",sep="_")] <- list(res_dacomp(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"dacomp",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error dacomp in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"zicoseq",sep="_")] <- list(res_zicoseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"zicoseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error zicoseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"limma_voom",sep="_")] <- list(res_limma_voom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"limma_voom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error limma_voom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"deseq2",sep="_")] <- list(res_deseq2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"deseq2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error deseq2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"corncob",sep="_")] <- list(res_corncob(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"corncob",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error corncob in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"maaslin2",sep="_")] <- list(res_maaslin2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"maaslin2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error maaslin2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name, "new", sep = "_")] <- list(res_new(data = sim))
            end_time <- Sys.time()
            time_diff <- difftime(end_time, start_time, units = "secs")
            res_sim[[paste(rl_name, "new", sep = "_")]]$time <- as.numeric(time_diff)
          },
          error = function(e){cat(paste0("error new in ",rl_name))})
          
          cat(paste0(rl_name," completed!!! \n"))
          
          save(res_sim,file=paste0(save_path,paste(i,j,k,l,p,sep = "_"),".RData"))
        } 
      }
    }
  }
}


#######################
####OTU = 300####
setwd("/public/linzhf/simulation_13DAA")  # Set working directory
source('all_function.R')  # Load custom functions
source('simulatemseq_new.R')  # Load simulation functions

## Simulation parameters setup   
notu <- 300  # Number of OTUs to simulate

# Create a save path based on the number of OTUs
save_path <- paste0("/public/linzhf/simulation_13DAA/simulation_", notu, "/")
if (dir.exists(save_path)) file.remove(save_path)  # Remove existing directory
if (!dir.exists(save_path)) dir.create(save_path)  # Create new directory

# Sample sizes
nsample <- c(20, 50, 100, 200, 300, 500)

# Types of OTUs
otutype <- c("abundant", "rare", "mix")

# Effect sizes for simulation
effectsize <- c(1, 2, 5, 10)

# Number of simulation iterations
B <- 100

# Load throat OTU table data
data(throat.otu.tab)  # Load the data

throat.otu.tab.tmp <- bind_cols(throat.otu.tab,throat.otu.tab)
names(throat.otu.tab.tmp) <- str_replace(names(throat.otu.tab.tmp),"[//.]+","_")
comm <- t(throat.otu.tab.tmp)

comm <- t(throat.otu.tab)  # Transpose the data
comm <- comm[rowMeans(comm != 0) > 0.05, ]  # Filter out low-abundance OTUs
comm <- comm[1:310,]

# Simulate binary covariate, 5% signal density, abundant differential OTUs, unbalanced change
# This setting simulates strong compositional effects
p.adj <- "fdr"  # Method for p-value adjustment
p.cutoff <- 0.05  # Significance threshold

set.seed(123)
for(i in nsample) {
  for(j in notu){
    for(k in otutype){
      for(l in effectsize){
        res_sim <- list()
        for(p in 1:B){
          sim.obj <- SimulateMSeq_new(
            ref.otu.tab = comm, 
            nSam = i, 
            nOTU = j,
            # True signal setting
            diff.otu.pct = 0.1, 
            diff.otu.direct = c("unbalanced"),
            diff.otu.mode = k,
            covariate.type = c("binary"), 
            grp.ratio = 1,
            covariate.eff.mean = log(l), 
            covariate.eff.sd = 0,
            # Confounder signal setting
            confounder.type = c("none"), 
            conf.cov.cor = 0.6,
            conf.diff.otu.pct = 0.1, 
            conf.nondiff.otu.pct = 0.1,
            confounder.eff.mean = 1.0, 
            confounder.eff.sd = 0,
            # Depth setting
            depth.mu = 10000, 
            depth.theta = 5, 
            depth.conf.factor = 0
          )
          #####################################
          sim <- list(data=sim.obj[["otu.tab.sim"]],
                      group = sim.obj[["covariate"]] %>% as.data.frame() %>% setNames("group"),
                      gold_name = sim.obj$otu.names[sim.obj$diff.otu.ind]
          )
          
          colnames(sim$data) <- paste0("sample",1:dim(sim$data)[2])
          rownames(sim$group)<- paste0("sample",1:dim(sim$data)[2])
          
          physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
                             sample_data(data.frame(sim$group)))
          
          tax_tab <- matrix(1,nrow = nrow(sim$data), ncol = 7)
          rownames(tax_tab) = rownames(sim$data)
          colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species")
          tax_tab <- as.data.frame(tax_tab)
          tax_tab$Kingdom <- rownames(sim$data)
          
          rl_name <- paste(i,j,k,l,p,sep = "_")
          
          tryCatch({
            start_time <- Sys.time() 
            res_sim[paste(rl_name,"wilcox",sep="_")] <- list(res_wilcox(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"wilcox",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error wilcox in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"linda",sep="_")] <- list(res_linda(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"linda",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error linda in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"metagenomeseq",sep="_")] <- list(res_metagenomeseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"metagenomeseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error metagenomeseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"lefse",sep="_")] <- list(res_lefse(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"lefse",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error lefse in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"rf",sep="_")] <- list(res_rf(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"rf",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error rf in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"edger",sep="_")] <- list(res_edger(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"edger",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error edger in ",rl_name))} )
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"ancombc2",sep="_")] <- list(res_ancombc2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"ancombc2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error ancombc2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_we",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_we",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_we in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_wi",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff, type="wi"))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_wi",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_wi in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"fastancom",sep="_")] <- list(res_fastancom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"fastancom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error fastancom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"dacomp",sep="_")] <- list(res_dacomp(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"dacomp",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error dacomp in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"zicoseq",sep="_")] <- list(res_zicoseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"zicoseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error zicoseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"limma_voom",sep="_")] <- list(res_limma_voom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"limma_voom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error limma_voom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"deseq2",sep="_")] <- list(res_deseq2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"deseq2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error deseq2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"corncob",sep="_")] <- list(res_corncob(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"corncob",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error corncob in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"maaslin2",sep="_")] <- list(res_maaslin2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"maaslin2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error maaslin2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name, "new", sep = "_")] <- list(res_new(data = sim))
            end_time <- Sys.time()
            time_diff <- difftime(end_time, start_time, units = "secs")
            res_sim[[paste(rl_name, "new", sep = "_")]]$time <- as.numeric(time_diff)
          },
          error = function(e){cat(paste0("error new in ",rl_name))})
          
          cat(paste0(rl_name," completed!!! \n"))
          
          save(res_sim,file=paste0(save_path,paste(i,j,k,l,p,sep = "_"),".RData"))
        } 
      }
    }
  }
}

#######################
####OTU = 400####
setwd("/public/linzhf/simulation_13DAA")  # Set working directory
source('all_function.R')  # Load custom functions
source('simulatemseq_new.R')  # Load simulation functions

## Simulation parameters setup   
notu <- 400  # Number of OTUs to simulate

# Create a save path based on the number of OTUs
save_path <- paste0("/public/linzhf/simulation_13DAA/simulation_", notu, "/")
if (dir.exists(save_path)) file.remove(save_path)  # Remove existing directory
if (!dir.exists(save_path)) dir.create(save_path)  # Create new directory

# Sample sizes
nsample <- c(20, 50, 100, 200, 300, 500)

# Types of OTUs
otutype <- c("abundant", "rare", "mix")

# Effect sizes for simulation
effectsize <- c(1, 2, 5, 10)

# Number of simulation iterations
B <- 100

# Load throat OTU table data
data(throat.otu.tab)  # Load the data

throat.otu.tab.tmp <- bind_cols(throat.otu.tab,throat.otu.tab)
names(throat.otu.tab.tmp) <- str_replace(names(throat.otu.tab.tmp),"[//.]+","_")
comm <- t(throat.otu.tab.tmp)
comm <- comm[1:410,]

comm <- t(throat.otu.tab)  # Transpose the data
comm <- comm[rowMeans(comm != 0) > 0.05, ]  # Filter out low-abundance OTUs

# Simulate binary covariate, 5% signal density, abundant differential OTUs, unbalanced change
# This setting simulates strong compositional effects
p.adj <- "fdr"  # Method for p-value adjustment
p.cutoff <- 0.05  # Significance threshold

set.seed(123)
for(i in nsample) {
  for(j in notu){
    for(k in otutype){
      for(l in effectsize){
        res_sim <- list()
        for(p in 1:B){
          sim.obj <- SimulateMSeq_new(
            ref.otu.tab = comm, 
            nSam = i, 
            nOTU = j,
            # True signal setting
            diff.otu.pct = 0.1, 
            diff.otu.direct = c("unbalanced"),
            diff.otu.mode = k,
            covariate.type = c("binary"), 
            grp.ratio = 1,
            covariate.eff.mean = log(l), 
            covariate.eff.sd = 0,
            # Confounder signal setting
            confounder.type = c("none"), 
            conf.cov.cor = 0.6,
            conf.diff.otu.pct = 0.1, 
            conf.nondiff.otu.pct = 0.1,
            confounder.eff.mean = 1.0, 
            confounder.eff.sd = 0,
            # Depth setting
            depth.mu = 10000, 
            depth.theta = 5, 
            depth.conf.factor = 0
          )
          #####################################
          sim <- list(data=sim.obj[["otu.tab.sim"]],
                      group = sim.obj[["covariate"]] %>% as.data.frame() %>% setNames("group"),
                      gold_name = sim.obj$otu.names[sim.obj$diff.otu.ind]
          )
          
          colnames(sim$data) <- paste0("sample",1:dim(sim$data)[2])
          rownames(sim$group)<- paste0("sample",1:dim(sim$data)[2])
          
          physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
                             sample_data(data.frame(sim$group)))
          
          tax_tab <- matrix(1,nrow = nrow(sim$data), ncol = 7)
          rownames(tax_tab) = rownames(sim$data)
          colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species")
          tax_tab <- as.data.frame(tax_tab)
          tax_tab$Kingdom <- rownames(sim$data)
          
          rl_name <- paste(i,j,k,l,p,sep = "_")
          
          tryCatch({
            start_time <- Sys.time() 
            res_sim[paste(rl_name,"wilcox",sep="_")] <- list(res_wilcox(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"wilcox",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error wilcox in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"linda",sep="_")] <- list(res_linda(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"linda",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error linda in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"metagenomeseq",sep="_")] <- list(res_metagenomeseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"metagenomeseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error metagenomeseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"lefse",sep="_")] <- list(res_lefse(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"lefse",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error lefse in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"rf",sep="_")] <- list(res_rf(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"rf",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error rf in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"edger",sep="_")] <- list(res_edger(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"edger",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error edger in ",rl_name))} )
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"ancombc2",sep="_")] <- list(res_ancombc2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"ancombc2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error ancombc2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_we",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_we",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_we in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_wi",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff, type="wi"))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_wi",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_wi in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"fastancom",sep="_")] <- list(res_fastancom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"fastancom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error fastancom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"dacomp",sep="_")] <- list(res_dacomp(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"dacomp",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error dacomp in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"zicoseq",sep="_")] <- list(res_zicoseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"zicoseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error zicoseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"limma_voom",sep="_")] <- list(res_limma_voom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"limma_voom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error limma_voom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"deseq2",sep="_")] <- list(res_deseq2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"deseq2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error deseq2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"corncob",sep="_")] <- list(res_corncob(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"corncob",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error corncob in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"maaslin2",sep="_")] <- list(res_maaslin2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"maaslin2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error maaslin2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name, "new", sep = "_")] <- list(res_new(data = sim))
            end_time <- Sys.time()
            time_diff <- difftime(end_time, start_time, units = "secs")
            res_sim[[paste(rl_name, "new", sep = "_")]]$time <- as.numeric(time_diff)
          },
          error = function(e){cat(paste0("error new in ",rl_name))})
          
          cat(paste0(rl_name," completed!!! \n"))
          
          save(res_sim,file=paste0(save_path,paste(i,j,k,l,p,sep = "_"),".RData"))
        } 
      }
    }
  }
}

#######################
####OTU = 500####
setwd("/public/linzhf/simulation_13DAA")  # Set working directory
source('all_function.R')  # Load custom functions
source('simulatemseq_new.R')  # Load simulation functions

## Simulation parameters setup   
notu <- 500  # Number of OTUs to simulate

# Create a save path based on the number of OTUs
save_path <- paste0("/public/linzhf/simulation_13DAA/simulation_", notu, "/")
if (dir.exists(save_path)) file.remove(save_path)  # Remove existing directory
if (!dir.exists(save_path)) dir.create(save_path)  # Create new directory

# Sample sizes
nsample <- c(20, 50, 100, 200, 300, 500)

# Types of OTUs
otutype <- c("abundant", "rare", "mix")

# Effect sizes for simulation
effectsize <- c(1, 2, 5, 10)

# Number of simulation iterations
B <- 100

# Load throat OTU table data
data(throat.otu.tab)  # Load the data

throat.otu.tab.tmp <- bind_cols(throat.otu.tab,throat.otu.tab)
names(throat.otu.tab.tmp) <- str_replace(names(throat.otu.tab.tmp),"[//.]+","_")
comm <- t(throat.otu.tab.tmp)
comm <- comm[1:510,]

comm <- t(throat.otu.tab)  # Transpose the data
comm <- comm[rowMeans(comm != 0) > 0.05, ]  # Filter out low-abundance OTUs

# Simulate binary covariate, 5% signal density, abundant differential OTUs, unbalanced change
# This setting simulates strong compositional effects
p.adj <- "fdr"  # Method for p-value adjustment
p.cutoff <- 0.05  # Significance threshold

set.seed(123)
for(i in nsample) {
  for(j in notu){
    for(k in otutype){
      for(l in effectsize){
        res_sim <- list()
        for(p in 1:B){
          sim.obj <- SimulateMSeq_new(
            ref.otu.tab = comm, 
            nSam = i, 
            nOTU = j,
            # True signal setting
            diff.otu.pct = 0.1, 
            diff.otu.direct = c("unbalanced"),
            diff.otu.mode = k,
            covariate.type = c("binary"), 
            grp.ratio = 1,
            covariate.eff.mean = log(l), 
            covariate.eff.sd = 0,
            # Confounder signal setting
            confounder.type = c("none"), 
            conf.cov.cor = 0.6,
            conf.diff.otu.pct = 0.1, 
            conf.nondiff.otu.pct = 0.1,
            confounder.eff.mean = 1.0, 
            confounder.eff.sd = 0,
            # Depth setting
            depth.mu = 10000, 
            depth.theta = 5, 
            depth.conf.factor = 0
          )
          #####################################
          sim <- list(data=sim.obj[["otu.tab.sim"]],
                      group = sim.obj[["covariate"]] %>% as.data.frame() %>% setNames("group"),
                      gold_name = sim.obj$otu.names[sim.obj$diff.otu.ind]
          )
          
          colnames(sim$data) <- paste0("sample",1:dim(sim$data)[2])
          rownames(sim$group)<- paste0("sample",1:dim(sim$data)[2])
          
          physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
                             sample_data(data.frame(sim$group)))
          
          tax_tab <- matrix(1,nrow = nrow(sim$data), ncol = 7)
          rownames(tax_tab) = rownames(sim$data)
          colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species")
          tax_tab <- as.data.frame(tax_tab)
          tax_tab$Kingdom <- rownames(sim$data)
          
          rl_name <- paste(i,j,k,l,p,sep = "_")
          
          tryCatch({
            start_time <- Sys.time() 
            res_sim[paste(rl_name,"wilcox",sep="_")] <- list(res_wilcox(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"wilcox",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error wilcox in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"linda",sep="_")] <- list(res_linda(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"linda",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error linda in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"metagenomeseq",sep="_")] <- list(res_metagenomeseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"metagenomeseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error metagenomeseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"lefse",sep="_")] <- list(res_lefse(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"lefse",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error lefse in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"rf",sep="_")] <- list(res_rf(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"rf",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error rf in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"edger",sep="_")] <- list(res_edger(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"edger",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error edger in ",rl_name))} )
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"ancombc2",sep="_")] <- list(res_ancombc2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"ancombc2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error ancombc2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_we",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_we",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_we in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_wi",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff, type="wi"))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_wi",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_wi in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"fastancom",sep="_")] <- list(res_fastancom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"fastancom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error fastancom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"dacomp",sep="_")] <- list(res_dacomp(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"dacomp",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error dacomp in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"zicoseq",sep="_")] <- list(res_zicoseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"zicoseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error zicoseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"limma_voom",sep="_")] <- list(res_limma_voom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"limma_voom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error limma_voom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"deseq2",sep="_")] <- list(res_deseq2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"deseq2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error deseq2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"corncob",sep="_")] <- list(res_corncob(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"corncob",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error corncob in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"maaslin2",sep="_")] <- list(res_maaslin2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"maaslin2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error maaslin2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name, "new", sep = "_")] <- list(res_new(data = sim))
            end_time <- Sys.time()
            time_diff <- difftime(end_time, start_time, units = "secs")
            res_sim[[paste(rl_name, "new", sep = "_")]]$time <- as.numeric(time_diff)
          },
          error = function(e){cat(paste0("error new in ",rl_name))})
          
          cat(paste0(rl_name," completed!!! \n"))
          
          save(res_sim,file=paste0(save_path,paste(i,j,k,l,p,sep = "_"),".RData"))
        } 
      }
    }
  }
}

#######################
####OTU = 800####
setwd("/public/linzhf/simulation_13DAA")  # Set working directory
source('all_function.R')  # Load custom functions
source('simulatemseq_new.R')  # Load simulation functions

## Simulation parameters setup   
notu <- 800  # Number of OTUs to simulate

# Create a save path based on the number of OTUs
save_path <- paste0("/public/linzhf/simulation_13DAA/simulation_", notu, "/")
if (dir.exists(save_path)) file.remove(save_path)  # Remove existing directory
if (!dir.exists(save_path)) dir.create(save_path)  # Create new directory

# Sample sizes
nsample <- c(20, 50, 100, 200, 300, 500)

# Types of OTUs
otutype <- c("abundant", "rare", "mix")

# Effect sizes for simulation
effectsize <- c(1, 2, 5, 10)

# Number of simulation iterations
B <- 100

# Load throat OTU table data
data(throat.otu.tab)  # Load the data

throat.otu.tab.tmp <- bind_cols(throat.otu.tab,throat.otu.tab)
names(throat.otu.tab.tmp) <- str_replace(names(throat.otu.tab.tmp),"[//.]+","_")
comm <- t(throat.otu.tab.tmp)
comm <- comm[1:810,]

comm <- t(throat.otu.tab)  # Transpose the data
comm <- comm[rowMeans(comm != 0) > 0.05, ]  # Filter out low-abundance OTUs

# Simulate binary covariate, 5% signal density, abundant differential OTUs, unbalanced change
# This setting simulates strong compositional effects
p.adj <- "fdr"  # Method for p-value adjustment
p.cutoff <- 0.05  # Significance threshold

set.seed(123)
for(i in nsample) {
  for(j in notu){
    for(k in otutype){
      for(l in effectsize){
        res_sim <- list()
        for(p in 1:B){
          sim.obj <- SimulateMSeq_new(
            ref.otu.tab = comm, 
            nSam = i, 
            nOTU = j,
            # True signal setting
            diff.otu.pct = 0.1, 
            diff.otu.direct = c("unbalanced"),
            diff.otu.mode = k,
            covariate.type = c("binary"), 
            grp.ratio = 1,
            covariate.eff.mean = log(l), 
            covariate.eff.sd = 0,
            # Confounder signal setting
            confounder.type = c("none"), 
            conf.cov.cor = 0.6,
            conf.diff.otu.pct = 0.1, 
            conf.nondiff.otu.pct = 0.1,
            confounder.eff.mean = 1.0, 
            confounder.eff.sd = 0,
            # Depth setting
            depth.mu = 10000, 
            depth.theta = 5, 
            depth.conf.factor = 0
          )
          #####################################
          sim <- list(data=sim.obj[["otu.tab.sim"]],
                      group = sim.obj[["covariate"]] %>% as.data.frame() %>% setNames("group"),
                      gold_name = sim.obj$otu.names[sim.obj$diff.otu.ind]
          )
          
          colnames(sim$data) <- paste0("sample",1:dim(sim$data)[2])
          rownames(sim$group)<- paste0("sample",1:dim(sim$data)[2])
          
          physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
                             sample_data(data.frame(sim$group)))
          
          tax_tab <- matrix(1,nrow = nrow(sim$data), ncol = 7)
          rownames(tax_tab) = rownames(sim$data)
          colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species")
          tax_tab <- as.data.frame(tax_tab)
          tax_tab$Kingdom <- rownames(sim$data)
          
          rl_name <- paste(i,j,k,l,p,sep = "_")
          
          tryCatch({
            start_time <- Sys.time() 
            res_sim[paste(rl_name,"wilcox",sep="_")] <- list(res_wilcox(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"wilcox",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error wilcox in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"linda",sep="_")] <- list(res_linda(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"linda",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error linda in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"metagenomeseq",sep="_")] <- list(res_metagenomeseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"metagenomeseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error metagenomeseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"lefse",sep="_")] <- list(res_lefse(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"lefse",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error lefse in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"rf",sep="_")] <- list(res_rf(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"rf",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error rf in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"edger",sep="_")] <- list(res_edger(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"edger",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error edger in ",rl_name))} )
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"ancombc2",sep="_")] <- list(res_ancombc2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"ancombc2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error ancombc2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_we",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_we",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_we in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_wi",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff, type="wi"))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_wi",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_wi in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"fastancom",sep="_")] <- list(res_fastancom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"fastancom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error fastancom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"dacomp",sep="_")] <- list(res_dacomp(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"dacomp",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error dacomp in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"zicoseq",sep="_")] <- list(res_zicoseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"zicoseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error zicoseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"limma_voom",sep="_")] <- list(res_limma_voom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"limma_voom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error limma_voom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"deseq2",sep="_")] <- list(res_deseq2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"deseq2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error deseq2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"corncob",sep="_")] <- list(res_corncob(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"corncob",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error corncob in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"maaslin2",sep="_")] <- list(res_maaslin2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"maaslin2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error maaslin2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name, "new", sep = "_")] <- list(res_new(data = sim))
            end_time <- Sys.time()
            time_diff <- difftime(end_time, start_time, units = "secs")
            res_sim[[paste(rl_name, "new", sep = "_")]]$time <- as.numeric(time_diff)
          },
          error = function(e){cat(paste0("error new in ",rl_name))})
          
          cat(paste0(rl_name," completed!!! \n"))
          
          save(res_sim,file=paste0(save_path,paste(i,j,k,l,p,sep = "_"),".RData"))
        } 
      }
    }
  }
}

#######################
####OTU = 1000####
setwd("/public/linzhf/simulation_13DAA")  # Set working directory
source('all_function.R')  # Load custom functions
source('simulatemseq_new.R')  # Load simulation functions

## Simulation parameters setup   
notu <- 1000  # Number of OTUs to simulate

# Create a save path based on the number of OTUs
save_path <- paste0("/public/linzhf/simulation_13DAA/simulation_", notu, "/")
if (dir.exists(save_path)) file.remove(save_path)  # Remove existing directory
if (!dir.exists(save_path)) dir.create(save_path)  # Create new directory

# Sample sizes
nsample <- c(20, 50, 100, 200, 300, 500)

# Types of OTUs
otutype <- c("abundant", "rare", "mix")

# Effect sizes for simulation
effectsize <- c(1, 2, 5, 10)

# Number of simulation iterations
B <- 100

# Load throat OTU table data
data(throat.otu.tab)  # Load the data

throat.otu.tab.tmp <- bind_cols(throat.otu.tab,throat.otu.tab)
names(throat.otu.tab.tmp) <- str_replace(names(throat.otu.tab.tmp),"[//.]+","_")
comm <- t(throat.otu.tab.tmp)
comm <- comm[1:1010,]

comm <- t(throat.otu.tab)  # Transpose the data
comm <- comm[rowMeans(comm != 0) > 0.05, ]  # Filter out low-abundance OTUs

# Simulate binary covariate, 5% signal density, abundant differential OTUs, unbalanced change
# This setting simulates strong compositional effects
p.adj <- "fdr"  # Method for p-value adjustment
p.cutoff <- 0.05  # Significance threshold

set.seed(123)
for(i in nsample) {
  for(j in notu){
    for(k in otutype){
      for(l in effectsize){
        res_sim <- list()
        for(p in 1:B){
          sim.obj <- SimulateMSeq_new(
            ref.otu.tab = comm, 
            nSam = i, 
            nOTU = j,
            # True signal setting
            diff.otu.pct = 0.1, 
            diff.otu.direct = c("unbalanced"),
            diff.otu.mode = k,
            covariate.type = c("binary"), 
            grp.ratio = 1,
            covariate.eff.mean = log(l), 
            covariate.eff.sd = 0,
            # Confounder signal setting
            confounder.type = c("none"), 
            conf.cov.cor = 0.6,
            conf.diff.otu.pct = 0.1, 
            conf.nondiff.otu.pct = 0.1,
            confounder.eff.mean = 1.0, 
            confounder.eff.sd = 0,
            # Depth setting
            depth.mu = 10000, 
            depth.theta = 5, 
            depth.conf.factor = 0
          )
          #####################################
          sim <- list(data=sim.obj[["otu.tab.sim"]],
                      group = sim.obj[["covariate"]] %>% as.data.frame() %>% setNames("group"),
                      gold_name = sim.obj$otu.names[sim.obj$diff.otu.ind]
          )
          
          colnames(sim$data) <- paste0("sample",1:dim(sim$data)[2])
          rownames(sim$group)<- paste0("sample",1:dim(sim$data)[2])
          
          physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
                             sample_data(data.frame(sim$group)))
          
          tax_tab <- matrix(1,nrow = nrow(sim$data), ncol = 7)
          rownames(tax_tab) = rownames(sim$data)
          colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species")
          tax_tab <- as.data.frame(tax_tab)
          tax_tab$Kingdom <- rownames(sim$data)
          
          rl_name <- paste(i,j,k,l,p,sep = "_")
          
          tryCatch({
            start_time <- Sys.time() 
            res_sim[paste(rl_name,"wilcox",sep="_")] <- list(res_wilcox(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"wilcox",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error wilcox in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"linda",sep="_")] <- list(res_linda(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"linda",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error linda in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"metagenomeseq",sep="_")] <- list(res_metagenomeseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"metagenomeseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error metagenomeseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"lefse",sep="_")] <- list(res_lefse(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"lefse",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error lefse in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"rf",sep="_")] <- list(res_rf(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"rf",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error rf in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"edger",sep="_")] <- list(res_edger(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"edger",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error edger in ",rl_name))} )
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"ancombc2",sep="_")] <- list(res_ancombc2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"ancombc2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error ancombc2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_we",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_we",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_we in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_wi",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff, type="wi"))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_wi",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_wi in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"fastancom",sep="_")] <- list(res_fastancom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"fastancom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error fastancom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"dacomp",sep="_")] <- list(res_dacomp(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"dacomp",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error dacomp in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"zicoseq",sep="_")] <- list(res_zicoseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"zicoseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error zicoseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"limma_voom",sep="_")] <- list(res_limma_voom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"limma_voom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error limma_voom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"deseq2",sep="_")] <- list(res_deseq2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"deseq2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error deseq2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"corncob",sep="_")] <- list(res_corncob(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"corncob",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error corncob in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"maaslin2",sep="_")] <- list(res_maaslin2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"maaslin2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error maaslin2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name, "new", sep = "_")] <- list(res_new(data = sim))
            end_time <- Sys.time()
            time_diff <- difftime(end_time, start_time, units = "secs")
            res_sim[[paste(rl_name, "new", sep = "_")]]$time <- as.numeric(time_diff)
          },
          error = function(e){cat(paste0("error new in ",rl_name))})
          
          cat(paste0(rl_name," completed!!! \n"))
          
          save(res_sim,file=paste0(save_path,paste(i,j,k,l,p,sep = "_"),".RData"))
        } 
      }
    }
  }
}


#######################
####OTU = 2000####
setwd("/public/linzhf/simulation_13DAA")  # Set working directory
source('all_function.R')  # Load custom functions
source('simulatemseq_new.R')  # Load simulation functions

## Simulation parameters setup   
notu <- 2000  # Number of OTUs to simulate

# Create a save path based on the number of OTUs
save_path <- paste0("/public/linzhf/simulation_13DAA/simulation_", notu, "/")
if (dir.exists(save_path)) file.remove(save_path)  # Remove existing directory
if (!dir.exists(save_path)) dir.create(save_path)  # Create new directory

# Sample sizes
nsample <- c(20, 50, 100, 200, 300, 500)

# Types of OTUs
otutype <- c("abundant", "rare", "mix")

# Effect sizes for simulation
effectsize <- c(1, 2, 5, 10)

# Number of simulation iterations
B <- 100

# Load throat OTU table data
data(throat.otu.tab)  # Load the data

throat.otu.tab.tmp <- bind_cols(throat.otu.tab,throat.otu.tab,throat.otu.tab)
names(throat.otu.tab.tmp) <- str_replace(names(throat.otu.tab.tmp),"[//.]+","_")
comm <- t(throat.otu.tab.tmp)
comm <- comm[1:2010,]

comm <- t(throat.otu.tab)  # Transpose the data
comm <- comm[rowMeans(comm != 0) > 0.05, ]  # Filter out low-abundance OTUs

# Simulate binary covariate, 5% signal density, abundant differential OTUs, unbalanced change
# This setting simulates strong compositional effects
p.adj <- "fdr"  # Method for p-value adjustment
p.cutoff <- 0.05  # Significance threshold

set.seed(123)
for(i in nsample) {
  for(j in notu){
    for(k in otutype){
      for(l in effectsize){
        res_sim <- list()
        for(p in 1:B){
          sim.obj <- SimulateMSeq_new(
            ref.otu.tab = comm, 
            nSam = i, 
            nOTU = j,
            # True signal setting
            diff.otu.pct = 0.1, 
            diff.otu.direct = c("unbalanced"),
            diff.otu.mode = k,
            covariate.type = c("binary"), 
            grp.ratio = 1,
            covariate.eff.mean = log(l), 
            covariate.eff.sd = 0,
            # Confounder signal setting
            confounder.type = c("none"), 
            conf.cov.cor = 0.6,
            conf.diff.otu.pct = 0.1, 
            conf.nondiff.otu.pct = 0.1,
            confounder.eff.mean = 1.0, 
            confounder.eff.sd = 0,
            # Depth setting
            depth.mu = 10000, 
            depth.theta = 5, 
            depth.conf.factor = 0
          )
          #####################################
          sim <- list(data=sim.obj[["otu.tab.sim"]],
                      group = sim.obj[["covariate"]] %>% as.data.frame() %>% setNames("group"),
                      gold_name = sim.obj$otu.names[sim.obj$diff.otu.ind]
          )
          
          colnames(sim$data) <- paste0("sample",1:dim(sim$data)[2])
          rownames(sim$group)<- paste0("sample",1:dim(sim$data)[2])
          
          physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
                             sample_data(data.frame(sim$group)))
          
          tax_tab <- matrix(1,nrow = nrow(sim$data), ncol = 7)
          rownames(tax_tab) = rownames(sim$data)
          colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species")
          tax_tab <- as.data.frame(tax_tab)
          tax_tab$Kingdom <- rownames(sim$data)
          
          rl_name <- paste(i,j,k,l,p,sep = "_")
          
          tryCatch({
            start_time <- Sys.time() 
            res_sim[paste(rl_name,"wilcox",sep="_")] <- list(res_wilcox(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"wilcox",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error wilcox in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"linda",sep="_")] <- list(res_linda(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"linda",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error linda in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"metagenomeseq",sep="_")] <- list(res_metagenomeseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"metagenomeseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error metagenomeseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"lefse",sep="_")] <- list(res_lefse(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"lefse",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error lefse in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"rf",sep="_")] <- list(res_rf(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"rf",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error rf in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"edger",sep="_")] <- list(res_edger(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"edger",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error edger in ",rl_name))} )
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"ancombc2",sep="_")] <- list(res_ancombc2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"ancombc2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error ancombc2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_we",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_we",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_we in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"aldex2_wi",sep="_")] <- list(res_aldex2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff, type="wi"))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"aldex2_wi",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error aldex2_wi in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"fastancom",sep="_")] <- list(res_fastancom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"fastancom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error fastancom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"dacomp",sep="_")] <- list(res_dacomp(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"dacomp",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error dacomp in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"zicoseq",sep="_")] <- list(res_zicoseq(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"zicoseq",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error zicoseq in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"limma_voom",sep="_")] <- list(res_limma_voom(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"limma_voom",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error limma_voom in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"deseq2",sep="_")] <- list(res_deseq2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"deseq2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error deseq2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"corncob",sep="_")] <- list(res_corncob(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"corncob",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error corncob in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name,"maaslin2",sep="_")] <- list(res_maaslin2(data=sim, p.adj = p.adj, p.cutoff = p.cutoff))
            end_time <- Sys.time()
            res_sim[[paste(rl_name,"maaslin2",sep="_")]]$time <- as.numeric(end_time - start_time) 
          },
          error = function(e){cat(paste0("error maaslin2 in ",rl_name))})
          
          tryCatch({
            start_time <- Sys.time()
            res_sim[paste(rl_name, "new", sep = "_")] <- list(res_new(data = sim))
            end_time <- Sys.time()
            time_diff <- difftime(end_time, start_time, units = "secs")
            res_sim[[paste(rl_name, "new", sep = "_")]]$time <- as.numeric(time_diff)
          },
          error = function(e){cat(paste0("error new in ",rl_name))})
          
          cat(paste0(rl_name," completed!!! \n"))
          
          save(res_sim,file=paste0(save_path,paste(i,j,k,l,p,sep = "_"),".RData"))
        } 
      }
    }
  }
}

