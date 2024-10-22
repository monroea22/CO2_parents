# CO<sub>2</sub> Parents 

### running on server, with bash scripts
##### check quality of fastq files
```{bash}
fastqc ./
```

##### combine into one multi view
```{bash}
multiqc ./
```

## trim files with trimmomatic
##### sliding window 4:20, min length of 40, phred 33

```{bash}
time java -jar /home/monroeaa/Programs/Trimmomatic-0.36/trim.jar PE -threads 2 rawreads/72-1-H5_L001_R1.fastq.gz rawreads/72-1-H5_L001_R2.fastq.gz -baseout 72h5.fq.gz ILLUMINACLIP:trim_adapt.fasta:2:30:10 HEADCROP:14 SLIDINGWINDOW:4:20 MINLEN:40

```

## hisat2
##### building index for alignment with apoly gtf file
```{bash}
hisat2_extract_splice_sites.py apoly.gtf > splicesites.tsv
hisat2_extract_exons.py apoly.gtf > exons.tsv
hisat2-build -p 10 --ss splicesites.tsv --exon exons.tsv apoly_primary_contigs.fasta.masked apoly

```
##### make txt file of sample names to run loop
```{bash}
ls co2parents/pairtrim | cut -d _ -f 1 | sort | uniq > samples.txt
```

##### running hisat2
```{bash}
 for f in $(<samples.txt); do hisat2 -p 8 --dta -x ~/index/apoly -1 ~/co2parents/pairtrim/${f}_1.fq.gz -2 ~/co2parents/pairtrim/${f}_2.fq.gz -S ${f}.sam; done
```

## running feature counts
##### running in same folder with bam files from hisat2
```{bash}
featureCounts -T 8 -t exon -g gene_id -p -s2 -Q 20 -B -a apoly_primary_geneannotation_v1.gtf -o co2par_2.txt *.bam
```

### Move over to R for differential expression

## DESeq2

#### load libraries
```{r}
library(DESeq2)
library(vegan)
library(rgl)
library(ape)
library(ggplot2)
library(pheatmap)
library(here)
```
### Importing counts table from feature counts
#### clean up columns and names to prep for DEseq2
```{r}
countdata <- read.table("co2par_2.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\.bam$", "", colnames(countdata))
colnames(countdata) <- gsub("X", "", colnames(countdata))
countdata <- as.matrix(countdata)
```
#### importing condition data including parent id and clutch id
##### double check read count table for accuracy
```{r}
coldata<-read.table('co2parentscond.txt', row.names=1)  ## cond and parents seperate
coldata<-coldata[,c("condition","parents","clutch")]
all(rownames(coldata) == colnames(countdata))
totalCounts=colSums(countdata)
totalCounts
mean(totalCounts)
```

#### create DEseq object from counts table and condition data and run LRT (likelihood ratio test) 
```{r}
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds$condition <- relevel(dds$condition, ref="control")
dds.c<- DESeq(dds, test="LRT",full=~condition, reduced=~1) ## genes due to condition
res.c<-results(dds.c)
summary(res.c)
res.c<-res.c[order(res.c$padj),]
table(res.c$padj<0.05) 
```
##### LRT based on parental identity
```{r}
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~parents)
dds.c<- DESeq(dds, test="LRT",full=~parents, reduced=~1) ## genes due to parental identity
res.c<-results(dds.c)
summary(res.c)
res.c<-res.c[order(res.c$padj),]
table(res.c$padj<0.05)
```

### PCA looking at impact of condition and parent assignment on DEG profiles
```{r}
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition+parents)
vst <- vst(dds, blind=FALSE)
pcData<-plotPCA(vst, intgroup=c("condition", "parents"), returnData=TRUE, ntop=20000)
percentVar <- round(100 * attr(pcData, "percentVar"))
pdf("pca4.pdf")
ggplot(pcData, aes(PC1, PC2, color= condition, shape=parents, group=condition)) +
  geom_point(size=4)  +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("#ca0020", "#f4a542","#5A4D4C","#0571B0" )) +
  coord_fixed() +
  scale_shape_manual(values=c(16,18))+
  theme_classic()+
  stat_ellipse(type="norm", linetype=2)
dev.off()
```
![<src="figures/pca4.pdf"/>](https://github.com/monroea22/CO2_parents/blob/main/figures/pca4.pdf)


### PCA with just control and transgenerational individuals (removing acute conditions)
```{r}
vsta <- vst[,vst$condition %in% c("control","transgen")]
pcData<-plotPCA(vsta, intgroup=c("condition", "parents"), returnData=TRUE, ntop=5000)
percentVar <- round(100 * attr(pcData, "percentVar"))
pdf("pca2.pdf")
ggplot(pcData, aes(PC1, PC2, color=parents, shape=condition, group=parents)) +
  geom_point(size=6)  +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c( "#507aa6", '#4f9e85')) +
  coord_fixed() +
  theme_classic() +
  stat_ellipse(type="norm", linetype=2)
dev.off()
```
![](/Users/monroea/Downloads/pca2.pdf)

### PCA with just acute conditions (control vs CO2)
```{r}
vstb <- vst[,vst$condition %in% c("acuteCO2","acuteCtrl")]
pcData<-plotPCA(vstb, intgroup=c("condition", "parents"), returnData=TRUE, ntop=5000)
percentVar <- round(100 * attr(pcData, "percentVar"))
pdf("pca3.pdf")
ggplot(pcData, aes(PC1, PC2, color=parents, shape=condition, group=parents)) +
  geom_point(size=6)  +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("#507aa6", '#4f9e85' )) +
  coord_fixed() +
  theme_classic() +
  stat_ellipse(type="norm", linetype=2) 
dev.off()
```
![](/Users/monroea/Downloads/pca3.pdf)


### LRT test looking at impact of clutch on DEGs
```{r}
coldata<-read.table('co2parentscondcomb.txt', row.names=1)
coldata<-coldata[,c("condition","clutch")]
all(rownames(coldata) == colnames(countdata))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~clutch)
dds.c<- DESeq(dds, test="LRT",full=~clutch, reduced=~1)
res.c<-results(dds.c)
summary(res.c)
res.c<-res.c[order(res.c$padj),]
table(res.c$padj<0.05)
```

## Running DESeq2 for comparisons
```{r}
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~clutch+condition)
dds<-dds[,-48]  ##removig 67b2 outlier
dds.c<- DESeq(dds)
res.c<-results(dds.c)
summary(res.c)
res.c<-res.c[order(res.c$padj),]
table(res.c$padj<0.05)
```
#### variance stabilized data of deseq object
```{r}
vsd=getVarianceStabilizedData(dds.c)
vals=cbind(res.c$pvalue, res.c$padj)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, 'co2par_vsdpvals_1outrem.csv', quote=F)
```

#### reading in genome annotation for functional comparisons
```{r}
xx<-read.csv('Apoly_2018_annotation3.csv', header = T)
row.names(xx)=xx[,1]
xx[,1]<-NULL
head(xx)
```
### comparison within treatments, between parent types


#### control
```{r}
res<-results(dds.c, contrast=c('condition','controlNTT', 'controlTNT'))
mcols(res,use.names=TRUE)
table(res$padj<0.05)
write.csv(res, file="ctrlNTTvsTNT.csv")## for GO_MWU
res<-read.csv("ctrlNTTvsTNT.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
head(res)
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "ctrlNTTvsTNT_sig.csv", quote=F)
```

#### trangenerational

```{r}
res<-results(dds.c, contrast=c('condition','transgenNTT', 'transgenTNT'))
mcols(res,use.names=TRUE)
table(res$padj<0.05)
write.csv(res, file="transgenNTTvsTNT.csv")## for GO_MWU
res<-read.csv("transgenNTTvsTNT.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
head(res)
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "transgenNTTvsTNT_sig.csv", quote=F)
```

#### acute control
```{r}
res<-results(dds.c, contrast=c('condition','acuteCtrlNTT', 'acuteCtrlTNT'))
mcols(res,use.names=TRUE)
table(res$padj<0.05)
write.csv(res, file="acuteCtrlNTTvsTNT.csv")## for GO_MWU
res<-read.csv("acuteCtrlNTTvsTNT.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
head(res)
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "acuteCtrlNTTvsTNT_sig.csv", quote=F)
```

#### acute CO2
```{r}
res<-results(dds.c, contrast=c('condition','acuteCO2NTT', 'acuteCO2TNT'))
mcols(res,use.names=TRUE)
table(res$padj<0.05)
write.csv(res, file="acuteCO2NTTvsTNT.csv")## for GO_MWU
res<-read.csv("acuteCO2NTTvsTNT.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
head(res)
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "acuteCO2NTTvsTNT_sig.csv", quote=F)
```

### Using the LFC shrink option to contrast all comparisons

#### within treatments, between parent types

control

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','controlNTT', 'controlTNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="ctrlNTTvsTNT_lfc.csv")## for GO_MWU
res<-read.csv("ctrlNTTvsTNT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
head(res)
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "ctrlNTTvsTNT_sigLFC.csv", quote=F)
```
transgenerational

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','transgenNTT', 'transgenTNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="transNTTvsTNT_lfc.csv")## for GO_MWU
res<-read.csv("transNTTvsTNT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
head(res)
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "transNTTvsTNT_sigLFC.csv", quote=F)
```

acute control

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','acuteCtrlNTT', 'acuteCtrlTNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="acuteCtrlNTTvsTNT_lfc.csv")## for GO_MWU
res<-read.csv("acuteCtrlNTTvsTNT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "acuteCtrlNTTvsTNT_sigLFC.csv", quote=F)
```

acute CO<sub>2</sub>


```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','acuteCO2NTT', 'acuteCO2TNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="acuteCO2NTTvsTNT_lfc.csv")## for GO_MWU
res<-read.csv("acuteCO2NTTvsTNT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "acuteCO2NTTvsTNT_sigLFC.csv", quote=F)
```
#### within parents, between treatments

NT♂T♀  Parents, control vs trangenerational

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','controlNTT', 'transgenNTT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="ctrlvstransNTT_lfc.csv")

res<-read.csv("ctrlvstransNTT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "ctrlvstransNTT_sigLFC.csv", quote=F)
```
T♂NT♀  parents, control vs transgenerational

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','controlTNT', 'transgenTNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="ctrlvstransTNT_lfc.csv")

res<-read.csv("ctrlvstransTNT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "ctrlvstransTNT_sigLFC.csv", quote=F)
```
NT♂T♀  acute control vs acute CO<sub>2</sub>

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','acuteCtrlNTT', 'acuteCO2NTT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="acuteCtrlvsacuteCO2NTT_lfc.csv")
sig<-resLFC[which(resLFC$padj<0.05),]
write.csv(sig, file="acuteCtrlvsacuteCO2NTT_sig.csv")
```
T♂NT♀  acute control vs acute CO<sub>2</sub>

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','acuteCtrlTNT', 'acuteCO2TNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="acuteCtrlvsacuteCO2TNT_lfc.csv")
sig<-resLFC[which(resLFC$padj<0.05),]
write.csv(sig, file="acuteCtrlvsacuteCO2TNT_sig.csv")
```
NT♂T♀  control vs acute control

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','controlNTT', 'acuteCtrlNTT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="ctrlvsacuteCtrlNTT_lfc.csv")
res<-read.csv("ctrlvsacuteCtrlNTT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "ctrlvsacuteCtrlNTT_sigLFC.csv", quote=F)
```
T♂NT♀  control vs acute control

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','controlTNT', 'acuteCtrlTNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="ctrlvsacuteCtrlTNT_lfc.csv")
res<-read.csv("ctrlvsacuteCtrlTNT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "ctrlvsacuteCtrlTNT_sigLFC.csv", quote=F)
```
T♂NT♀  control vs acute CO<sub>2</sub>

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','controlTNT', 'acuteCO2TNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="ctrlvsacuteCO2TNT_lfc.csv")
res<-read.csv("ctrlvsacuteCO2TNT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "ctrlvsacuteCO2TNT_sigLFC.csv", quote=F)
```
NT♂T♀  control vs acute CO<sub>2</sub>

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','controlNTT', 'acuteCO2NTT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="ctrlvsacuteCO2NTT_lfc.csv")
res<-read.csv("ctrlvsacuteCO2NTT_lfc.csv", header=T)
row.names(res)=res[,1]
res[,1]<-NULL
mervals<-merge(xx, res, by='row.names')
head(mervals)
mervals<-mervals[which(mervals$padj<0.05),]
nrow(mervals)
write.csv(mervals, "ctrlvsacuteCO2NTT_sigLFC.csv", quote=F)
```
T♂NT♀  transgenerational vs acute CO<sub>2</sub>

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','transgenTNT', 'acuteCO2TNT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="transvsacuteCO2TNT_lfc.csv")
sig<-resLFC[which(resLFC$padj<0.05),]
write.csv(sig, file="transvsacuteCO2TNT_sig.csv")
```
NT♂T♀  transgenerational vs acute CO<sub>2</sub>

```{r}
resLFC <- lfcShrink(dds.c, contrast=c('condition','transgenNTT', 'acuteCO2NTT'))
table(resLFC$padj<0.05)
write.csv(resLFC, file="transvsacuteCO2NTT_lfc.csv")
sig<-resLFC[which(resLFC$padj<0.05),]
write.csv(sig, file="transvsacuteCO2NTT_sig.csv")
```

## bubble graph GO terms
```{r}
data<-read.table("goTerm_co2par2.txt", header=T)
ggplot(data, aes(x = condition, y = go)) +        
	geom_point(aes(color = Delta_rank, size = FDR), alpha = 0.8) + ## alpha is opacity of bubbles
	scale_color_gradient2(low = '#00897B', mid = 'white', high = '#8E24AA') + ##green to purple
	scale_size(range = c(10, 4)) + ##bubble sizes, decreasing with increased FDR
	theme_classic() ## removes background
```
![](/Users/monroea/Downloads/bubbleGraph_goTerms.pdf)
