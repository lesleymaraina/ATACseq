library(Rsamtools)
library(Biostrings)
bamfile <- 'LY1_ATAC_chr18.bam'
bam <- scanBam(bamfile)
#look at fields in bam object
names(bam[[1]])
# we will work with the isize field: insert size, insertion size between adapters and is equivalent to fragment size
# look at distribution of isize
hist(bamm[[1]]$isize)
hist(bam[[1]]$isize)
fragsize <- bam[[1]]$isize [ bam[[1]]$isize > 1 & bam[[1]]$isize < 500]
hist(fragsize)
head(fragsize)
table(fragsize)
#count number of times you see each insert size: table(fragsize)
plot(table(fragsize))
#ATAC seq contains a lot of mitochondrial DNA
# remove pcr duplicates (collapse using peak calling) and mitochondrial DNA
# call peaks
library(ChIPpeakAnno)
source("http://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)
peaks <- toGRanges("ATAC_LY1_peaks_chr18.bed")
mean(peaks)
sum(peaks)
sum(peaks[[1]])
head(peaks)
length(peaks)
mean(width(peaks))
# width: gives mean of all geomean objects
# annotate chip seq
anno <- annotatePeakInBatch(peaks, AnnotationData=TSS.human.GRCh37,
output="overlapping", maxgap=1000L)
data(TSS.human.GRCh37)
anno <- annotatePeakInBatch(peaks, AnnotationData=TSS.human.GRCh37,
output="overlapping", maxgap=1000L)
?annotatePeakInBatcm\h
?annotatePeakInBatch
anno
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
anno <- addGeneIDs(annotatedPeak=anno,
orgAnn="org.Hs.eg.db",
IDs2Add="symbol")
anno
# characterize overlap, the most overlap occurs within the inside classification
table(as.data.frame(anno)$insideFeature)
pie(table(as.data.frame(anno)$insideFeature))
# overlap ATACseq with Chipseq -> strong overlap b/w ATACseq and ChipSeq peaks --> tells whether transcription factor is active in cell
# ChiP seq peak and open chromatin: gives a clue of txn binding
# ATACseq: gives data on open chromatin
PU1peaks <- toGRanges("PU1_LY1_peaks_chr18.bed")
# ChiPseq: can tell whether txn factor is bindng to DNA
PU1peaks <- toGRanges("PU1_LY1_peaks_chr18.bed")
length(PU1peaks)
annoPU1 <- annotatePeakInBatch(peaks, AnnotationData=PU1peaks)
annoPU1
annoPU1 <- annotatePeakInBatch(peaks, AnnotationData=PU1peaks, output="overlap")
annoPU1
ol <- findOverlapsOfPeaks(peaks, PU1peaks, maxgap=0)
peaklist <-ol$peaklist
names(peaklist)
length(peaklist)
length(peaklist[["peaks///PU1peaks"]])
length(peaklist[["peaks"]])
length(peaklist[["PU1peaks"]])
(length(peaklist[["peaks///PU1peaks"]]))/(length(peaklist[["peaks"]]))
(length(peaklist[["peaks///PU1peaks"]]))/(length(peaklist[["peaks"] + length(peaklist[["PU1peaks"]])]))
makeVennDiagram(ol, totalTest=1e+5)
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
library(Rsamtools)
library(Biostrings)
source("http://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")
bamfile <- 'LY1_ATAC_chr18.bam'
bam <- scanBam(bamfile)
fragsize <- bam[[1]]$isize [ bam[[1]]$isize > 1 & bam[[1]]$isize < 500]
head(fragsize)
plot(table(fragsize))
peaks <- toGRanges("ATAC_LY1_peaks_chr18.bed")
length(peaks)
mean(width(peaks))
data(TSS.human.GRCh37)
library(ChIPpeakAnno)
anno <- annotatePeakInBatch(peaks, AnnotationData=TSS.human.GRCh37,
output="overlapping", maxgap=1000L)
anno
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
anno <- addGeneIDs(annotatedPeak=anno,
orgAnn="org.Hs.eg.db",
IDs2Add="symbol")
anno
table(as.data.frame(anno)$insideFeature)
pie(table(as.data.frame(anno)$insideFeature))
PU1peaks <- toGRanges("PU1_LY1_peaks_chr18.bed")
length(PU1peaks)
annoPU1 <- annotatePeakInBatch(peaks, AnnotationData=PU1peaks)
annoPU1 <- annotatePeakInBatch(peaks, AnnotationData=PU1peaks, output="overlap")
ol <- findOverlapsOfPeaks(peaks, PU1peaks, maxgap=0)
peaklist <-ol$peaklist
names(peaklist)
makeVennDiagram(ol, totalTest=1e+5)
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
peaksWithSequences <- getAllPeakSequence(peaks, upstream=0,
downstream=0, genome=Hsapiens)
peaksWithSequences
s <- DNAStringSet(peaksWithSequences$sequence)
biocLite("MotifDb")
biocLite("seqLogo")
library(seqLogo)
library(MotifDb)
pfm.pu1 <- query(MotifDb,"Spib")[[1]]
pfm.pu1
seqLogo(pfm.pu1)
query(MotifDb, "Spib")[[3]]
sequo(query(MotifDb, "Spib")[[3]])
seqLogoo(query(MotifDb, "Spib")[[3]])
seqLogo(query(MotifDb, "Spib")[[3]])
?seqLogo
seqLogo(query(MotifDb, "Spib")[[3]],ic.scale=F)
seqLogo(query(MotifDb, "Spib")[[1]],ic.scale=F)
seqLogoo(query(MotifDb, "Spib")[[3]])
seqLogo(query(MotifDb, "Spib")[[3]])
pcm.pu1 <- round(100 * pfm.pu1)
pcm.pu1
mmpu1 <- matchPWM(pcm.pu1, unlist(s), "90%")
length(mmpu1)
mmpu1
# Motif Matching
length(mmpu1)
length(peaks)
peaks_shifted = shift(peaks, width(peaks))
peaksWithSequences_shifted <- getAllPeakSequence(peaks_shifted = shift(peaks, width(peaks)), upstream=0,
downstream=0, genome=Hsapiens)
peaksWithSequences_shifted <- getAllPeakSequence(peaks_shifted, upstream=0,
downstream=0, genome=Hsapiens)
ss <- DNAStringSet(peaksWithSequences_shifted$sequence)
mmpu1_s <- matchPWM(pcm.pu1, unlist(ss), "90%")
length(mmpu1_s)
#chi square test
prop.test(c(2406,2220),c(6383,6383))
prop.test(c(2406,2220),c(6383,6383), alternative="greater")
binom.test(2406, 6383 ,2220/6383)
binom.test(2406, 6383 ,2220/6383, alternative="greater")
