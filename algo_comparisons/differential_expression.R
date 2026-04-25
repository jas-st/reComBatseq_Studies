library(reComBatseq)
suppressMessages(library(limma))
suppressMessages(library(DESeq2))

deseq2.dea <- function(counts, colData, fileprefix) {
    dds <- DESeqDataSetFromMatrix(
        countData=counts, 
        colData=colData, 
        design=~Disease
    )
    dds$Disease <- relevel(dds$Disease, ref = "Healthy")
    dds <- DESeq(dds)
    resLFC <- lfcShrink(
        dds, 
        coef="Disease_Psoriasis_vs_Healthy", 
        type="apeglm"
    )
    write.table(
        resLFC, 
        file = paste0(fileprefix, '.dea.lfcshrink.tsv'),
        sep = '\t', 
        quote = F
    )
    write.table(
        results(dds), 
        file = paste0(fileprefix, '.dea.tsv'), 
        sep = '\t', 
        quote = F
    )
}


limma.dea <- function(data, group, fileprefix) {
    group <- as.factor(group)
    group <- relevel(group, ref = "Healthy")
    mm <- model.matrix(~group)
    fit <- lmFit(data, mm)
    efit <- eBayes(fit)
    write.table(
        topTable(efit, number = Inf), 
        file = paste0(fileprefix, '.dea.tsv'), 
        sep = '\t', 
        quote = F
    )
}


accessions <- c('SRP238713', 'SRP065812', 'SRP035988', 'SRP026042', 'ERP110816')

message('Reading raw data')
raw.data <- read.table(
    'psoriasis/psoriasis.all.exp.tsv',
    sep = '\t',
    quote = '',
    header = T,
    row.names = 'X'
)
message('Reading recombatseq data')
rcs.data <- read.table(
    'psoriasis/psoriasis.all.exp.recombatseq.tsv',
    sep = '\t',
    quote = '',
    header = T
)
message('Reading recombat data')
rc.data <- read.table(
    'psoriasis/psoriasis.all.exp.recombat.tsv',
    sep = '\t',
    quote = '',
    header = T,
    row.names = 'X'
)
message('Reading metadata')
meta.data <- read.table(
    'psoriasis/psoriasis.all.meta.tsv',
    sep = '\t',
    quote = '',
    header = T,
    row.names = 'X'
)
message('Performing DEA for full dataset')
deseq2.dea(t(raw.data), meta.data, 'psoriasis/psoriasis.all')
deseq2.dea(rcs.data, meta.data, 'psoriasis/psoriasis.all.recombatseq')
limma.dea(t(rc.data), meta.data$Disease, 'psoriasis/psoriasis.all.recombat')

for(acc in accessions) {
    message('Performing DEA for ', acc)
    boolidx <- meta.data$sra_study_acc == acc
    meta <- meta.data[boolidx,]
    
    raw <- raw.data[boolidx,]
    deseq2.dea(t(raw), meta, paste0('psoriasis/psoriasis.', acc))
    
    rcs <- rcs.data[,boolidx]
    deseq2.dea(rcs, meta, paste0('psoriasis/psoriasis.', acc, '.recombatseq'))
    
    rc <- rc.data[boolidx,]
    limma.dea(t(rc), meta$Disease, paste0('psoriasis/psoriasis.', acc, '.recombat'))
}