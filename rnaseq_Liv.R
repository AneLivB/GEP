###### RNA-seq data by Rebecca Nagel ######
###### Script by Bernice Sepers ######
###### Arctocephalus gazella 2018-2019 and 2019-2020 focal individuals ######


##-------------- create tx2gene file to be able to link gene_ID to transcript ID --------------##

library(tximport) # package for getting Kallisto results into R
library(readr)
library(data.table)
genes <- fread("/vol/cluster-data/bsepers/ref_genome/BIE3448_Antarctic_Fur_Seal.annotation.gff")
head(genes)
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
genes <- genes[type == "gene"]
unique(genes$chr)
unique(genes$source)
unique(genes$type)
unique(genes$score)
unique(genes$phase)
unique(genes$attributes)
genes$score <- NULL
genes$phase <- NULL
genes[, c("transcript_ID", "gene") := tstrsplit(attributes, ";", fixed=TRUE)]
head(genes$transcript_ID)
genes$transcript_ID<-gsub("ID=","",genes$transcript_ID)
head(genes$gene)
genes[, c("gene_ID", "protein") := tstrsplit(gene, ": ", fixed=TRUE)]
genes$gene<-NULL
head(genes$gene_ID)
genes$gene_ID<-gsub("Note=Similar to ","",genes$gene_ID)
head(genes$protein)
# it is very difficult to remove the different species because there are about 130 and there are multiple " (" in one cell
genes[, c("protein_ID", "species", "info", "info2") := tstrsplit(protein, " (", fixed=TRUE)]
find.list<-unique(grep(" OX%", genes$species, value = TRUE))
genes$protein_ID<-NULL
genes$species<-NULL
genes$info<-NULL
genes$info2<-NULL
find.string <- paste(unlist(find.list), collapse = "|\\(")
head(genes$protein)
tail(genes$protein)
genes$protein<-gsub(find.string, replacement = "", genes$protein)
genes[, c("gene", "species", "info") := tstrsplit(gene_ID, " (", fixed=TRUE)]
find.list<-unique(grep(" OX%", genes$species, value = TRUE))
genes$gene<-NULL
genes$species<-NULL
genes$info<-NULL
find.string <- paste(unlist(find.list), collapse = "|\\(")
head(genes$gene_ID)
tail(genes$gene_ID)
genes$gene_ID<-gsub(find.string, replacement = "", genes$gene_ID)
head(genes)
genes$gene_ID<-tolower(genes$gene_ID)
gene<-as.data.frame(genes$gene_ID)
gene<-unique(gene)
gene$n<-nchar(gene$`genes$gene_ID`)
dup<-as.data.frame(duplicated(gene$lower))
rm(gene,dup)
library(dplyr)
tx2g<-as.data.frame(select(genes, c('transcript_ID','gene_ID')))
rm(genes)
head(tx2g)
tx2g$target_id <- paste(tx2g$transcript_ID, "RA", sep="-") #add -RA to end of each transcript ID to match names in kallisto abundance files
tx2gene<-as.data.frame(select(tx2g, c('target_id','gene_ID')))
rm(tx2g)
length(unique(tx2gene$target_id)) #23408
length(unique(tx2gene$gene_ID)) #15569


##--------------  import Kallisto transcript counts into R using Tximport --------------##

## collect the paths to the kallisto abundance files (using a sample sheet)
targets <- read_tsv("/vol/cluster-data/bsepers/rnaseq/sample_sheet_rna.txt") #all RNA files
nrow(targets) #189
targets[targets == "NA"] <- NA
targets<-subset(targets, targets$sample != "F24_FWB_mum_start_2018B") #exclude outlier based on MDS plot (very low library size)
path <- file.path("/vol/cluster-data/bsepers/rnaseq/pseudoalign",targets$sample,"abundance.tsv")
all(file.exists(path)) #true, all files exist
#add sample ID to the row to make sure that the sample ID information is maintained when data is imported
sampleLabels <- targets$sample
head(path)
names(path)<-sampleLabels

txi_gene <- tximport(path,
                     type = "kallisto",
                     tx2gene = tx2gene,
                     dropInfReps=TRUE, #otherwise it will try to import .h5 files
                     txOut = FALSE) # TRUE = transcript or FALSE = gene level data

counts<-as.data.frame(txi_gene$counts) #kallisto counts in one dataframe
rm(tx2gene,txi_gene,targets,path,find.list,find.string,sampleLabels)


##-------------- FILTERING AND NORMALIZATION --------------##

library(edgeR) # loads limma as a dependency
library(DESeq2) # needed for rowCounts
d0 <- DGEList(counts) #create DGEList object
head(d0$counts)
head(d0$samples)
tail(d0$samples)
dim(d0) # 15569 genes, 188 samples

#add sample info
d0$samples$Sample <- row.names(d0$samples)
targets <- read.csv("sample_sheet_rna_Liv.csv", header = TRUE, sep = ";") 
targets[targets == "NA"] <- NA
colnames(targets)[2] <- "Sample"
targets<-subset(targets,targets$Sample != "F24_FWB_mum_start_2018B") # exclude outlier
sampleinfo <- left_join(d0$samples, targets, by="Sample")
d0$samples <- sampleinfo # add sample info 
d0$samples$group <- NULL # uninformative column
rownames(d0$samples) <- d0$samples$Sample
colnames(d0$samples)[13] <- "group" #replace time point by group
rm(sampleinfo,targets,counts)

#filter low-expressed genes, only include genes with cpm > 0.2 in at least 20 samples
#0.1 is not enough to remove the dip, 0.2 is
voom(d0, plot=TRUE) #there is a dip left side in the plot, indicative of a mean-variance relationship in lowly expressed genes
keep <- rowSums(cpm(d0$counts) >= 0.2) >= 20
d <- d0[keep,]
dim(d0) # 15569 genes, 188 samples
dim(d) # 12639 genes, 188 samples
voom(d, plot=TRUE)
rm(keep,d0)

#keep genes covered in at least 10 samples in at least one group in each comparison (colony, time point, lhs, year)
filter <- cpm(d$counts) >= 0.2
head(filter)
filt <- matrix(c(rowCounts(filter[, c(colnames(filter)[grepl("SSB",colnames(filter))==T])]),
                 rowCounts(filter[, c(colnames(filter)[grepl("FWB",colnames(filter))==T])]),
                 rowCounts(filter[, c(colnames(filter)[grepl("start",colnames(filter))==T])]),
                 rowCounts(filter[, c(colnames(filter)[grepl("end",colnames(filter))==T])]),
                 rowCounts(filter[, c(colnames(filter)[grepl("mum",colnames(filter))==T])]),
                 rowCounts(filter[, c(colnames(filter)[grepl("pup",colnames(filter))==T])]),
                 rowCounts(filter[, c(colnames(filter)[grepl("2018",colnames(filter))==T])]),
                 rowCounts(filter[, c(colnames(filter)[grepl("2019",colnames(filter))==T])])), ncol=8)
colnames(filt) <- c("SSB", "FWB", "start", "end",
                    "mum", "pup", "y2018", "y2019")
rownames(filt) <- rownames(d)                                                      
filt<-as.data.frame(filt)
filt$keep_col <- case_when(filt$SSB >= 10 | filt$FWB >= 10 ~ "keep",
                           filt$SSB < 10 & filt$FWB < 10 ~ "not_keep")
filt$keep_tp <- case_when(filt$start >= 10 | filt$end >= 10 ~ "keep",
                            filt$start < 10 & filt$end < 10 ~ "not_keep")
filt$keep_lhs <- case_when(filt$mum >= 10 | filt$pup >= 10 ~ "keep",
                            filt$mum < 10 & filt$pup < 10 ~ "not_keep")
filt$keep_year <- case_when(filt$y2018 >= 10 | filt$y2019 >= 10 ~ "keep",
                            filt$y2018 < 10 & filt$y2019 < 10 ~ "not_keep")
filtwhich.n <- subset(filt, filt$keep_col == "keep" & filt$keep_tp == "keep" &
                        filt$keep_lhs == "keep" & filt$keep_year == "keep") 
d_n <- d[rownames(d) %in% rownames(filtwhich.n),] # get gene names to keep & select these
dim(d_n) 
rm(filter, filt, filtwhich.n, d)
# unfiltered:                                                    15569 genes, 188 samples
# cutoff cpm > 0.2 in n>=10:                                     12917 genes, 188 samples
# cutoff cpm > 0.2 in n>=20:                                     12639 genes, 188 samples
# cutoff cpm > 0.2 in n>=10 in at least 1 group each comparison: 12679 genes, 188 samples
# cutoff cpm > 0.2 in n>=20 & cpm > 0.2 in n>=10 in at least 1 group per comparison: 12639 genes, 188 samples

#Calculate normalization factors for use downstream
d_n <- calcNormFactors(d_n, method="TMM")


##-------------- SAMPLE CLUSTERING --------------##

#Derive information from the sample names
snames <- rownames(d_n$samples) # Sample names
snames
ID_tp <- snames
ID_tp <- gsub("_SSB_pup_start_2018B", "_1", ID_tp)
ID_tp <- gsub("_SSB_pup_start_2019B", "_1", ID_tp)
ID_tp <- gsub("_FWB_pup_start_2018B", "_1", ID_tp)
ID_tp <- gsub("_FWB_pup_start_2019B", "_1", ID_tp)
ID_tp <- gsub("_FWB_pup_start_2018A", "_1", ID_tp)
ID_tp <- gsub("_FWB_pup_start_2018", "_1", ID_tp)
ID_tp <- gsub("_SSB_pup_end_2018B", "_2", ID_tp)
ID_tp <- gsub("_SSB_pup_end_2019B", "_2", ID_tp)
ID_tp <- gsub("_FWB_pup_end_2018B", "_2", ID_tp)
ID_tp <- gsub("_FWB_pup_end_2019B", "_2", ID_tp)
ID_tp <- gsub("_SSB_pup_end_2018", "_2", ID_tp)
ID_tp <- gsub("_SSB_mum_start_2018B", "_1", ID_tp)
ID_tp <- gsub("_SSB_mum_start_2019B", "_1", ID_tp)
ID_tp <- gsub("_FWB_mum_start_2018B", "_1", ID_tp)
ID_tp <- gsub("_FWB_mum_start_2019B", "_1", ID_tp)
ID_tp <- gsub("_FWB_mum_start_2018A", "_1", ID_tp)
ID_tp <- gsub("_FWB_mum_start_2018", "_1", ID_tp)
ID_tp <- gsub("_SSB_mum_end_2018B", "_2", ID_tp)
ID_tp <- gsub("_SSB_mum_end_2019B", "_2", ID_tp)
ID_tp <- gsub("_FWB_mum_end_2018B", "_2", ID_tp)
ID_tp <- gsub("_FWB_mum_end_2019B", "_2", ID_tp)
ID_tp <- gsub("_SSB_mum_end_2018", "_2", ID_tp)
rm(ID_tp)
ID <- substr(ID_tp, 1, nchar(ID_tp) -2)
ID <- as.factor(ID)
year <- substr(snames, nchar(snames) - 4, nchar(snames))
year[grep(x=year,pattern="2018")]<-"2018"
year[grep(x=year,pattern="2019")]<-"2019"
year <- as.factor(year)
timepoint <- substr(snames, nchar(snames) - 10, nchar(snames) - 5)
timepoint[grep(x=timepoint,pattern="end")]<-"end"
timepoint[grep(x=timepoint,pattern="start")]<-"start"
timepoint <- as.factor(timepoint)
colony <- substr(snames, 4, nchar(snames) - 13)
colony[grep(x=colony,pattern="FWB")]<-"FWB"
colony[grep(x=colony,pattern="SSB")]<-"SSB"
colony <- as.factor(colony)
batch <- d_n$samples$rna_batch
batch <- as.factor(batch)
lhs <- substr(snames, 8, nchar(snames) - 9)
lhs[grep(x=lhs,pattern="mum")]<-"mum"
lhs[grep(x=lhs,pattern="pup")]<-"pup"
lhs <- as.factor(lhs)

#multidimensional scaling plots
plotMDS(d_n, col = as.numeric(timepoint)) # clear grouping by timepoint within lhs
plotMDS(d_n, col = as.numeric(batch))
plotMDS(d_n, col = as.numeric(year)) # grouping by year within mums
plotMDS(d_n, col = as.numeric(colony))
plotMDS(d_n, col = as.numeric(lhs)) # clear grouping by lhs

library("tibble")
d_n$samples <- add_column(d_n$samples, "ID_tp" = paste(d_n$samples$ID, d_n$samples$group, sep = "_"), .before = "group")
plotMDS(cpm(d_n, log=TRUE), labels = d_n$samples$ID_tp, 
        col = as.numeric(d_n$samples$sample_collection_Dec_days))
plotMDS(cpm(d_n, log=TRUE), labels = d_n$samples$ID_tp, 
        col = as.numeric(d_n$samples$group)) 
plotMDS(cpm(d_n, log=TRUE), labels = d_n$samples$ID_tp, 
        col = as.numeric(d_n$samples$delta_days_samples+1)) 
plotMDS(cpm(d_n, log=TRUE), labels = d_n$samples$ID_tp, 
        col = as.numeric(year))
plotMDS(cpm(d_n, log=TRUE), labels = d_n$samples$ID_tp, 
        col = as.numeric(colony))
plotMDS(cpm(d_n, log=TRUE), labels = d_n$samples$ID_tp, 
        col = as.numeric(batch))
plotMDS(cpm(d_n, log=TRUE), labels = d_n$samples$ID_tp, 
        col = as.numeric(lhs))

#test if there is an association between the loadings and the grouping factors
mds<-plotMDS(cpm(d_n, log=TRUE))
library(lme4)
library(lmerTest)
x.id<-lm(mds$x ~ as.factor(ID))
anova(x.id) #significant difference between individuals
y.id<-lm(mds$y ~ as.factor(ID))
anova(y.id)#significant difference between individuals
x.col<-lmer(mds$x ~ as.factor(colony) + (1|ID))
anova(x.col) #NS
y.col<-lmer(mds$y ~ as.factor(colony) + (1|ID))
anova(y.col) #NS
x.bat<-lmer(mds$x ~ as.factor(batch) + (1|ID))
anova(x.bat) #NS
y.bat<-lmer(mds$y ~ as.factor(batch) + (1|ID))
anova(x.bat) #NS
x.year<-lmer(mds$x ~ as.factor(year) + (1|ID))
anova(x.year) #NS
y.year<-lmer(mds$y ~ as.factor(year) + (1|ID))
anova(y.year) #significant difference between years
x.tp<-lmer(mds$x ~ as.factor(timepoint) + (1|ID))
anova(x.tp)
y.tp<-lmer(mds$y ~ as.factor(timepoint) + (1|ID))
anova(y.tp) # significant difference between timepoints
x.days<-lmer(mds$x ~ as.numeric(d_n$samples$delta_days_samples+1) + (1|ID))
anova(x.days) #significant
y.days<-lmer(mds$y ~ as.numeric(d_n$samples$delta_days_samples+1) + (1|ID))
anova(y.days) #significant
x.tp_days<-lmer(mds$x ~ as.factor(timepoint) + as.numeric(d_n$samples$delta_days_samples+1) + (1|ID))
anova(x.tp_days) # both significant
y.tp_days<-lmer(mds$y ~ as.factor(timepoint) +as.numeric(d_n$samples$delta_days_samples+1) + (1|ID))
anova(y.tp_days) 
x.dec<-lmer(mds$x ~ as.numeric(d_n$samples$sample_collection_Dec_days+1) + (1|ID))
anova(x.dec) #significant
y.dec<-lmer(mds$y ~ as.numeric(d_n$samples$sample_collection_Dec_days+1) + (1|ID))
anova(y.dec) #significant, but it is highly correlated with age/delta days
# it is really hard to differentiate between timepoint effects, (within)seasonal effects, age effects and days between sampling effects
rm(list=ls(pattern="x."))
rm(list=ls(pattern="y."))
rm(mds,batch,colony,ID,lhs,snames,timepoint)

saveRDS(d_n, file = "rnaseq_filtered_Liv.Rdata")


##-------------- ANALYSIS --------------##

d_n<-readRDS("rnaseq_filtered_Liv.Rdata")
d_n$samples$group <- as.factor(d_n$samples$group) # group = timepoint (start or end)
d_n$samples$ID <- as.factor(d_n$samples$ID)
design <- model.matrix(~ -1 + group, data = d_n$samples) #setting up model contrasts is more straight forward in the absence of an intercept for timepoint
contr.matrix <- makeContrasts(
  end_vs_start = group2-group1,
  levels = colnames(design)) #create matrix for contrasts
# Use voom to remove variance dependency on mean, see https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from d itself. 
#  If filtering of lowly-expressed genes is insufficient, a drop in variance levels can be observed at the low end of the expression scale due to very small counts.
#voom() converts the read counts to log2-cpm, with associated weights, ready for linear modelling
vobj_tmp <- voom(d_n, design, plot=TRUE)
# Fit a random effect
dupcor <- duplicateCorrelation(vobj_tmp, design, block = d_n$samples$ID)
dupcor$consensus #the estimated intra individual correlation
#the intra ind correlation will change the voom weights slightly
# run voom again considering the duplicateCorrelation results in order to compute more accurate precision weights
vobj = voom(d_n, design, plot=TRUE, 
            block=d_n$samples$ID, correlation=dupcor$consensus)
#update the correlation for the new voom weights
dupcor <- duplicateCorrelation(vobj, design, block = d_n$samples$ID)
dupcor$consensus # the estimated intra individual correlation

#analyzing repeated measures data using duplicateCorrelation.
#The model forces the magnitude of the random effect to be the same across all genes.
# Estimate linear mixed model with a single variance component
# Fit the model for each gene
# But this step uses only the genome-wide average for the random effect
vfit <- lmFit(vobj, design, block=d_n$samples$ID, correlation=dupcor$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
# Fit Empirical Bayes for moderated t-statistics, empirical Bayes moderation is carried out by borrowing information across all the genes to obtain more precise estimates of gene-wise variability
efit <- eBayes(vfit)
topTable(efit, n=20)
plotSA(efit, main="Final model: Mean-variance trend")
plotMD(efit)
abline(h=0,col="darkgrey")

#adjusted p-value cutoff that is set at 5% by default.
efit$F
head(efit$F.p.value) #pval of the moderated F-statistics (equivalent to one-way ANOVA), it combines all contrast, so in this case it is the same as the t-statistic
head(efit$p.value) #same
summary(decideTests(efit, adjust.method="fdr"))
#set log-fold-changes (log-FCs) to be above a minimum value. The treat method (McCarthy and Smyth 2009) used to calculate p-values from empirical Bayes moderated t-statistics with a minimum log-FC requirement. 
#no need to use the ebayes output for this
lfc<-0.14
tfit <- treat(vfit, lfc=lfc) #(means that the expression of that gene is increased by a multiplicative factor of 2^0.14 ≈ 1.101905).
topTreat(tfit)
nrow(tfit)
dt <- decideTests(tfit, adjust.method="fdr")
summary(dt) 
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1])

nrow(tfit)
results<-as.data.frame(topTreat(tfit, n=12594))
results_pups2<-subset(results_pups, results_pups$logFC >= 0.14)
results_pups3<-subset(results_pups, results_pups$logFC <= -0.14)
result<-rbind(results_pups2,results_pups3)
result$adj.P.r <- round(result$adj.P.Val, digits = 2)
result2<-subset(result, result$adj.P.r < 0.05)
nrow(subset(result2, result2$logFC > 0))
nrow(subset(result2, result2$logFC < 0))