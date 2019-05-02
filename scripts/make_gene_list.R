library(mTEC.10x.pipeline)
library(Single.mTEC.Transcriptomes)
library(mTEC.10x.data)
library(dropseqr)
library(dplyr)
library(DropletUtils)
library(grid)
library(Biobase)
library(GEOquery)
library(limma)

data("tras")
data("aireDependentSansom")
data("fantom")
data("geneNames")
data("biotypes")
proteinCoding <- names( which( biotype == "protein_coding" ) )

stage_color_df <- data.frame("cTEC" = "#CC6600", "Immature" = "#009933",
                            "Intermediate" = "#0066CC", "Mature" = "#660099",
                            "Late_mature" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#666666")

stage_color <- t(stage_color_df)[ , 1]


working_dir <- "/home/kwells4/mTEC_dev/RankL_ablation/"
full_gene_list <- Seurat::Read10X(data.dir = paste0(working_dir, 
                      "isoControlBeg/isoControlBeg_count/outs/filtered_gene_bc_matrices/mm10/"))

##################
# Set gene lists #
##################

#######################################################
# All genes defined as all original genes used in 10x #
#######################################################
full_gene_list <- rownames(full_gene_list)

############################
# Protein coding gene list #
############################
proteinCoding <- unique(proteinCoding)
isProteinCodingEns <- gene_to_ensembl$V1 %in% proteinCoding
isProteinCodingAll <- gene_to_ensembl[isProteinCodingEns, ]
isProteinCodingGene <- unique(isProteinCodingAll$V2)

####################################################
# TRAs as defined in part one of single mTEC study #
####################################################
tras <- unique(tras$'gene.ids')
isTRA_ens <- gene_to_ensembl$V1 %in% tras
isTRA_all <- gene_to_ensembl[isTRA_ens, ]
isTRA_gene <- unique(isTRA_all$V2)
isTRA_gene <- intersect(isTRA_gene, isProteinCodingGene)

############################################################
# Aire dependent genes downloaded from Sansom study (2014) #
############################################################

aireDependent <- aireDependentSansom
aireDependent_ens <- gene_to_ensembl$V1 %in% aireDependent
aireDependent_all <- gene_to_ensembl[aireDependent_ens, ]
aireDependent_gene <- unique(aireDependent_all$V2)
aire_tra <- intersect(aireDependent_gene, isTRA_gene)
aire_tra <- intersect(aire_tra, isProteinCodingGene)

#########################
# Genes not in TRA list #
#########################

non_tra_gene <- isProteinCodingGene[!isProteinCodingGene %in% isTRA_gene]

#######################################################################
# 30 Aire dependent genes from review (DOI: 10.1016/j.it.2017.07.010) #
#######################################################################

aire_dependent_short <- c("Amelx", "Dppa3", "Ins1", "Klk6", "Mpo", "Ngp",
                          "Sprr2f", "Scgb1a1", "Tac2", "Camk2b", "Gsta2",
                          "Igf2", "Kng1", "Mt4", "Npy", "Slpi", "S100a8",
                          "Tff2", "C4bp", "Hbb-y", "Krt2", "Ly6g", "Mup4",
                          "Procr", "Spt1", "Sprr1b", "Zp3")

########################################################################
# 30 Fezf2 dependent genes from review (DOI: 10.1016/j.it.2017.07.010) #
########################################################################

fezf2_dependent <- c("Anxa10", "Apoc2", "Car8", "Cyp24a1", "Fabp9", "Krt10",
                     "Lypd1", "Myo15b", "Resp18", "Apoa4", "Bhmt", "Col17a1",
                     "F2", "Gc", "Kcnj5", "Maoa", "Plagl1", "Uox", "Apob",
                     "Calca", "Crisp1", "Fabp7", "Itih3", "Lgals7", "Muc3",
                     "Pld1", "Zp2")

############################
# TRA based on FANTOM data #
############################
mtecCombinedSub <- mtecCombined
mtecCombinedSub@assay$ablation_DE <- NULL
mtecSub <- Seurat::SubsetData(mtecCombinedSub, ident.use = mtec_subset)
count_table <- mtecSub@data

# Figure out exactly what these functions are doing
meansFANTOM <- sapply( split(seq_len(ncol(dxdFANTOM)),
                      colData( dxdFANTOM )$tissue), function(x){
                      rowMeans(
                      counts(dxdFANTOM, normalized=TRUE)[,x, drop=FALSE] )
                })
meansFANTOM <- sapply(
                       split( seq_len(
                       nrow(meansFANTOM) ),
                       sapply( strsplit( rownames( meansFANTOM ), "," ), "[[", 1 )),
                              function(x){
                              colMeans( meansFANTOM[x,,drop=FALSE] )
                     })

meansFANTOM <- t( meansFANTOM )
cat( sprintf("The total number of tissues used from the FANTOM dataset was:%s\n",
length( unique( colnames(meansFANTOM) ) ) ) )
## The total number of tissues used from the FANTOM dataset was:91
matchedIndexes <- match( geneNames, rownames(meansFANTOM))
stopifnot(
         rownames(meansFANTOM)[matchedIndexes[!is.na(matchedIndexes)]] ==
         geneNames[!is.na(matchedIndexes)] )
rownames(meansFANTOM)[matchedIndexes[!is.na(matchedIndexes)]] <-
               names( geneNames[!is.na(matchedIndexes)] )
meansFANTOM <- meansFANTOM[grep("ENS", rownames(meansFANTOM)),]

gene_df <- as.data.frame(geneNames)

meansFantom_gene <- merge(gene_df, meansFANTOM, by = "row.names")
rownames(meansFantom_gene) <- meansFantom_gene$geneNames
meansFantom_gene$Row.names <- NULL 
meansFantom_gene$geneNames <- NULL

# This is numbers of tissues for ALL genes
numbersOfTissues <- rowSums( meansFantom_gene > 5 )

count_table <- as.matrix(count_table)
count_table_expr <- count_table[rowSums(count_table) > 0, ]

numbersOfTissues <- numbersOfTissues[names(numbersOfTissues) %in%
                                     rownames(count_table)]
aireDependent_gene <- aireDependent_gene[aireDependent_gene %in%
                                    rownames(count_table)]

aireDependent_gene <- intersect(aireDependent_gene, isProteinCodingGene)

table( numbersOfTissues < 5 )

TRA_fantom <- names(numbersOfTissues[numbersOfTissues < 5])
TRA_fantom <- intersect(TRA_fantom, isProteinCodingGene)

##################
# Non fantom TRA #
##################

non_fantom_tra_gene <- isProteinCodingGene[!isProteinCodingGene %in% TRA_fantom]

# This is number of tissues for only DE genes from their data set
# Why did they choose to only look at DE genes here?
#numbersOfTissues <- numbersOfTissues[names(numbersOfTissues) %in% deGenesNone]
#aireDependent <- aireDependent[aireDependent %in% deGenesNone]

#################
# Fezf2 targets #
#################

# This code came directly from: https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE69105
# It is from Geo2r. Probably worth rechecking at some point
# load series and platform data from GEO

gset <- getGEO("GSE69105", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1261", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "0000011111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0) ||
          (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
# Remove cutoff of 250
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = 30000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

# This code is my own.
# Cutoff at 1.5 and pval 0.05. as in this paper
# http://dx.doi.org/10.1016/j.it.2017.07.010 1
# But many of their genes seem to insead have a cutoff of 1.
# This is consistent with a 2 fold change used to define Aire genes,
# So I will go with that
tT <- tT[tT$logFC < -1,]
tT <- tT[tT$adj.P.Val < 0.05, ]
fezf2_gene <- unique(tT$Gene.symbol)
fezf2_gene <- unlist(strsplit(fezf2_gene, "///"))
fezf2_gene <- unique(fezf2_gene)
fezf2_gene <- intersect(fezf2_gene, isProteinCodingGene)

###################
# Aire fantom TRA #
###################
aire_tra_fantom <- intersect(aireDependent_gene, TRA_fantom)

#############
# Fezf2 TRA #
#############
fezf2_tra <- intersect(fezf2_gene, isTRA_gene)

####################
# Fezf2 TRA fantom #
####################
fezf2_tra_fantom <- intersect(fezf2_gene, TRA_fantom)

########################
# All not listed genes #
########################
all_gene_list <- c(TRA_fantom, aireDependent_gene, fezf2_gene)
other_protein <- isProteinCodingGene[!isProteinCodingGene %in% all_gene_list]


###### 

# Add all gene lists you want to analyze here
gene_lists <- list(protein_coding = isProteinCodingGene, non_tra = non_fantom_tra_gene,
                   all_other_genes = other_protein,
                   tra_brennecke = isTRA_gene, tra_fantom = TRA_fantom,
                   aire_genes = aireDependent_gene, aire_tra = aire_tra_fantom,
                   fezf2_genes = fezf2_gene, fezf2_tra = fezf2_tra_fantom)


################################################################################
save(gene_lists, file = "/home/kwells4/mTEC_dev/data/gene_lists.rda")