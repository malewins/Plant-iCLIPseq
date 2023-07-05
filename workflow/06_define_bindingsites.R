################################################################################
###                                                                          ###
###                 Binding site definition (Timing 2 h)                     ###
###                                                                          ###
################################################################################
#
# Supplementary R-Script for the binding site definition. Step 233 - 248.


# 233.	Load the BindingSiteFinder package and additional packages
library("BindingSiteFinder") 
library("rtracklayer") 
library("GenomicFeatures") 
library("tidyr") 
library("dplyr") 

# 234.	Read the PureCLIP-called peaks and the crosslink signal (bigWig) from each repli-cate into the working memory.
peaksInitial = "05_called_peaks/AtGRP7-GFP.bed"
peaksInitial = import(con = peaksInitial, format = "BED", 
                      extraCols=c("additionalScores" = "character"))
clipFiles = "04_mapped/01_deduped_BAMs/beds/01_bigwigs"
clipFiles = list.files(clipFiles, pattern = ".bw$", full.names = TRUE)
clipFilesP = clipFiles[grep(clipFiles, pattern = "plus")]
clipFilesM = clipFiles[grep(clipFiles, pattern = "minus")]

# 235.	Create the initial BindingSiteFinder object from both resources and provide meta information for each sample.
colData = data.frame(
    name = "AtGRP7",
    id = c(1:3),
    condition = factor(c("WT"), levels = c("WT")),
    clPlus = clipFilesP,
    clMinus = clipFilesM)
bdsInitial = BSFDataSetFromBigWig(ranges = peaksInitial, meta = colData)

# 236.	Globally remove peaks with the lowest 2% PureCLIP score by applying the PureCLIP global filter function, to reduce noise
bds = pureClipGlobalFilter(bdsInitial, cutoff = 0.02)

# 237.	Read gene annotations from the AtRTD3 and Ensembl into memory and link them both by using the unique geneID. 
anno = rtracklayer::import(
    con = "annotation/atRTD3_TS_21Feb22_transfix.gtf",
    format = "gtf")
ensAnno = rtracklayer::import(
    con = "annotation/Arabidopsis_thaliana.TAIR10.56.gff3",
    format = "gff3")
mcols(anno)$name = ensAnno$Name[match(anno$gene_id, ensAnno$gene_id)]
mcols(anno)$geneType = ensAnno$biotype[match(anno$gene_id, ensAnno$gene_id)]

# 238.	Use the AtRTD3 annotation to create a TxDB database that can be queried for gene and transcript regions.
annoDb = makeTxDbFromGFF(
    file = "annotation/atRTD3_TS_21Feb22_transfix.gtf", 
    format = "gtf", circ_seqs = character()) 

# 239.	Extract positions from all genes in the annotation.
gns = genes(annoDb)
idx = match(gns$gene_id, anno$gene_id)
elementMetadata(gns) = cbind(elementMetadata(gns),
                             elementMetadata(anno)[idx,]) 

# 240.	Use the estimateBsWidth() function to simultaneously estimate the parameters for binding site size and gene-wise filtering. 
bds = estimateBsWidth(object = bds, anno.genes = gns, est.maxBsWidth = 19,
                      geneResolution = "medium", bsResolution = "medium",
                      est.subsetChromosome = "Chr1")

# 241.	Apply the gene-wise filter to selectively retain the top 70% of all PureCLIP sites for each gene, while eliminating the lowest 30%.
bds = pureClipGeneWiseFilter(bds, anno.genes = gns, cutoff = 0.3, 
                             overlaps = "keepSingle")

# 242.	Utilizing the estimated binding site width of 5 nt and compute a set of equally sized binding sites using the makeBindingSites() function. 
bds = makeBindingSites(bds, bsSize = 5, minWidth = 3)

# 243.	Apply the replicate reproducibility filter. 
bds = reproducibilityFilter(bds, cutoff = 0.1, nReps = 2)

# 244.	Assign binding sites to their host genes based on overlaps. 
bds = assignToGenes(bds, anno.genes = gns,
                    overlaps = "frequency",
                    match.geneID = "gene_id",
                    match.geneName = "name",
                    match.geneType = "geneType")

# 245.	Extract specific transcript regions 
cdseq = cds(annoDb)
intrns = unlist(intronsByTranscript(annoDb))
utrs3 = unlist(threeUTRsByTranscript(annoDb))
utrs5 = unlist(fiveUTRsByTranscript(annoDb))
trl = GRangesList(CDS = cdseq, Intron = intrns, UTR3 = utrs3, UTR5 = utrs5)

# 246.	Assign binding sites to their hosting transcript regions. 
bds = assignToTranscriptRegions(bds, anno.transcriptRegionList = trl,
                                overlaps = "frequency")

# 247.	Re-assign PureCLIP scores to each binding site by using the highest score from all overlapping PureCLIP sites.
bds = annotateWithScore(bds, match.ranges = peaksInitial, match.option = "max")

# 248.	Export a bed file containing the transcriptome-wide binding sites of the processed sample. 
exportToBED(bds, con = "06_binding_sites/AtGRP7-GFP.bsites.bed")

