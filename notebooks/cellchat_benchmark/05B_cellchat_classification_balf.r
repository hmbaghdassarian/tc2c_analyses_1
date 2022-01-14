suppressPackageStartupMessages({
    library(CellChat, quietly = T)
    library(patchwork, quietly = T)
    library(RhpcBLASctl, quietly = T)
    library(Matrix, quietly = T)
    library(data.table, quietly = T)
})
options(stringsAsFactors = FALSE)
source('embedding_seed_utils.r') # for non-stochastic embeddings
# RhpcBLASctl::blas_set_num_threads(25) # no multithreading

# paths
expression_data_path = '/data2/hratch/immune_CCI/covid/expression_data/covid_data/'#'/data2/eric/CCC-Benchmark/data/External/'
external_expression_path = F # set to T if using commented out path , small parsing differences

output_path = '/data2/hratch/immune_CCI/covid/balf_classification/'
rev_path = '/data3/hratch/tc2c_analyses_1/natcomm_revisions/'
input_data_path = '/data2/hratch/immune_CCI/covid/inputs/'

# parameters
type_<-'functional'#'structural'
group<-0
version<-1 #0<--didn't have seeds properly set, used for classifier (see markdown for details)
seed<-888L
set.seed(seed)

if (external_expression_path){
    cell_grouper<-'cell_type'
}else{
    cell_grouper<-'celltype'
}

humandb<-CellChatDB.human

fns = list()
for (fn in list.files(expression_data_path)){
    sn = strsplit(fn, '_')[[1]]
    sample.name = sn[[2]]
    type = sn[[1]]
    fns[[sample.name]][[type]] = fn
}
sample.names<-names(fns)

# load the UMI counts
read_sample_csv<-function(sn){
    counts<-as.data.frame(fread(paste0(expression_data_path, fns[[sn]]$DGE)))
    rownames(counts)<-counts$Gene
    counts<-counts[ , !(colnames(counts) %in% c('Gene'))]
    if(external_expression_path){colnames(counts)<-paste0(sn, '.', colnames(counts))} # consistency with metadata
    return(as.matrix(counts))
}

read_meta<-function(sn){
    meta<-read.csv(paste0(expression_data_path, fns[[sn]]$Meta))
    rownames(meta) = meta$Cell
    return(meta)
}


if (!group){                             
    counts<-lapply(setNames(sample.names, sample.names), function(sn) read_sample_csv(sn))
    meta<-lapply(setNames(sample.names, sample.names), function(sn) read_meta(sn))
                                  
}else{ # group by context
    stop('Should not group by context for classification')        
}

md.cell<-do.call("rbind", meta)
if (!external_expression_path){rownames(md.cell)<-md.cell$ID}
if (type_ == 'functional'){# functional requires same cell types to work                
    # filter for intersection of cell types across samples/contexts - to make comparable with Tensor-cell2cell
    cell.types<-Reduce(intersect, lapply(counts, function(df) unique(md.cell[colnames(df), cell_grouper])))
    cell.ids<-rownames(md.cell[md.cell[[cell_grouper]] %in% cell.types, ])   
    for (n in names(counts)){
        df<-counts[[n]]
        counts[[n]]<-df[,colnames(df) %in% cell.ids]
    }
}       

suppressWarnings({
    suppressMessages({
        # create cellchat object for each sample or sample.name
        covid.list<-list()
        for (sample.name in names(counts)){
            print(sample.name)
            # loop through each sample.name and create a cell type future
            expr<-CellChat::normalizeData(counts[[sample.name]])
            cellchat<-createCellChat(object = as(expr, "dgCMatrix"), meta = md.cell[colnames(expr),], 
                                           group.by = cell_grouper)
            cellchat@DB <- humandb # human organism

            cellchat <- subsetData(cellchat) # subset the expression data of signaling genes, assign to @data.signalling 
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat) # generate @ LR slot used by computeCommunProb
            cellchat <- projectData(cellchat, PPI.human) # shallow sequencing depth

            cellchat <- computeCommunProb(cellchat, raw.use = F, type = 'triMean', trim = NULL, seed.use = seed, 
                                         population.size = F) 

            # The functional similarity analysis requires the same cell population composition between two datasets.
            cellchat <- filterCommunication(cellchat, min.cells = 10)
            cellchat <- computeCommunProbPathway(cellchat)
            covid.list[[sample.name]]<-cellchat
        }
        print('merge')
        # merge and analyze
        cellchat <- mergeCellChat(covid.list, add.names = names(covid.list))
        cellchat <- computeNetSimilarityPairwise(cellchat, type = type_)
        cellchat.embed <- netEmbedding2(cellchat, type = type_, seed.use = seed)
    })
})
print('complete')

# for natcomm_revisions/01_cellchat_covid_balf.ipynb
saveRDS(cellchat.embed, paste0(rev_path, 'interim/covid_balf_cellchat_', type_, '_',
                               version, '.rds'))