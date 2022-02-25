suppressPackageStartupMessages({
    library(CellChat, quietly = T)
    library(liana, quietly = T)
    library(Seurat, quietly = T)
    library(data.table, quietly = T)
    library(dplyr, quietly = T)
    library(magrittr, quietly = T)
})
options(stringsAsFactors = FALSE)

rev_path = '/data3/hratch/tc2c_analyses_1/natcomm_revisions/'
expr_files = '/data2/eric/Tensor-Revisions/COVID-19-BALF-log1p.h5ad'

raw_counts_path<-'/data2/hratch/immune_CCI/covid/expression_data/covid_data/'
# normalized_counts_path<-paste0(rev_path, 'interim/tc2c_external_inputs/liana/liana_inputs2/')
n.cores<-20
seed=888

# rewrites CellChat::subsetCommunicationProbability to not remove 0 scores
# basically, we want all scores regardless of p-values and if they're 0
format.cellchat.communication<-function(cellchat.obj, thresh = NULL){
    prob <- cellchat.obj@net$prob
    pval <- cellchat.obj@net$pval
    prob[pval >= thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source","target","interaction_name")
    net.pval <- reshape2::melt(pval, value.name = "pval")
    net$pval <- net.pval$pval
    
    pairLR <- dplyr::select(cellchat.obj@LR$LRsig, c("interaction_name_2", "pathway_name", "ligand",  "receptor" ,"annotation","evidence"))
    idx <- match(net$interaction_name, rownames(pairLR))
    net <- cbind(net, pairLR[idx,])
    return(net)
    }

# adapts liana::call_cellchat to avoid CellChat::subsetCommunicationProbability, which excludes 0 scores
# default params from Daniel Dimitrov for consistency with LIANA
run_cellchat<-function(so, seed, nboot=1000, expr_prop=0.1, de_thresh=1, thresh=NULL){
    labels <- Seurat::Idents(so)
    meta <- data.frame(group = labels, row.names = names(labels))
    cellchat.obj<-createCellChat(object = GetAssayData(so, assay = 'RNA', slot = "data"), 
                                     meta = meta, 
                                   group.by = "group")
    cellchat.obj <- CellChat::addMeta(cellchat.obj, meta = meta)
    cellchat.obj <- CellChat::setIdent(cellchat.obj, ident.use = "group")
    
    ccDB<-CellChatDB.human
    resource<-c('CellChatDB')
    resource %<>% select_resource
    ccDB<-liana::cellchat_formatDB(ccDB, op_resource=resource$CellChatDB, exclude_anns=c())
    cellchat.obj@DB <- ccDB# human organism
    
    cellchat.obj <- subsetData(cellchat.obj) # subset the expression data of signaling genes, assign to @data.signalling 
    cellchat.obj <- identifyOverExpressedGenes(cellchat.obj, thresh.pc = expr_prop,
                                                thresh.p = de_thresh)
    cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
    cellchat.obj <- projectData(cellchat.obj, PPI.human)
    cellchat.obj <- computeCommunProb(cellchat.obj, raw.use = F, type = 'triMean', trim = NULL, 
                                      seed.use = seed, population.size=T, do.fast=T, nboot=nboot) 
    cellchat.obj <- filterCommunication(cellchat.obj, min.cells = 1) # as in LIANA
    
    cm<-format.cellchat.communication(cellchat.obj, thresh=thresh) # replace CellChat::subsetCommunicationProbability
    cm <- cm %>%
        dplyr::select(source,
               target,
               ligand,
               receptor,
               prob,
               pval) %>%
        as_tibble()
    return(cm)
    
}

# get samples
samples<-c()
for (fn in list.files(raw_counts_path)){
    samples<-c(samples, strsplit(fn, '_')[[1]][[2]])
}
samples<-unique(samples)


# get scores
for (sample in samples){
    print(sample)
    raw.counts<-read.csv(paste0(raw_counts_path, 'DGE_', sample, '_External_Tool.csv'), row.names = 1)
    md<-read.csv(paste0(raw_counts_path, 'Meta_', sample, '_External_Tool.csv'), row.names = 1)
    
    so<-CreateSeuratObject(counts = raw.counts, meta.data = md)
    so<-Seurat::NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1e6)
    Idents(so)<-'celltype'
    
    cellchat_score<-run_cellchat(so, seed = seed)
    
    # write the communication scores
    fwrite(cellchat_score, 
          paste0(rev_path, 'interim/tc2c_external_inputs/liana/liana_outputs/', sample, '_communication_scores_', 
                 'cellchat', '.csv')
          )
    }
}