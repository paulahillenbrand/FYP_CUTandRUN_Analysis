***R CODE*** 
### Set up RStudio
## libraries 

library(DESeq2)
library(dplyr)
library(ggpubr) 
library(corrplot)
library(stringr)
library(ggplot2)
library(viridis)

## Set Project Path
projPath = "/Users/paula/Downloads/cutandrun"
sampleList = c("AP_1_1", "AP_1_2", "AP_2_1", "AP_2_2","K27_1_1", "K27_1_2", "Pre_1_1", "Pre_1_2", "Pre_2_1", "Pre_2_2")
histList = c("AP", "K27", "Pre")
Pre-Processing 
Quality Control
fastqc {samp_name}.fq.gz
multiqc /cutandrun/rawdata
multiqc /cutandrun/trimmed

### Correlation Heatmap

reprod = c()
fragCount = NULL
for(hist in sampleList){
  if(is.null(fragCount)){
    fragCount = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
  }else{
    fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}
M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 
corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", title = "Sample Correlation Heatmap", mar=c(0,0,2,0), addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 1, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))

### Alignment results 

alignResult = c()
for(hist in sampleList){
  alignRes = read.table(paste0(projPath, "/alignment/sam/alignment_txt/", hist, "_bowtie2.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  alignResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           AlignmentRate_HSymV2.1 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
alignResult$Histone = factor(alignResult$Histone, levels = histList)
alignResult %>% mutate(AlignmentRate_HSymV2.1 = paste0(AlignmentRate_HSymV2.1, "%"))

## plot 
fig3A = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 30) +
    ylab("Sequencing Depth per Million") +
    xlab("") + 
    ggtitle("A. Sequencing Depth")

fig3B = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_HSymV2.1, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 30) +
    ylab("% of Mapped Fragments") +
    xlab("") +
    ggtitle("B. Alignment Rate (HSymV2.1)")
ggarrange(fig3A, fig3B, ncol = 2, common.legend = TRUE, legend="bottom")


### Duplication Rate

dupResult = c()
for(hist in sampleList){
  dupRes = read.table(paste0(projPath, "/alignment/sam/dupMark/", hist, "_picard.dupMark.txt"), header = TRUE, fill = TRUE)

  histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric, DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}
dupResult$Histone = factor(dupResult$Histone, levels = histList)

## Plot

fig4A = dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Duplication Rate (*100%)") +
    xlab("") +
    ylim(0,100) +
  ggtitle("A. Duplication Rate")

fig4B = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Unique Fragments per Million") +
    xlab("") +
    ylim(0,20)+
  ggtitle("C. Unique Fragments")
ggarrange(fig4A, fig4Bnrow=1, common.legend = TRUE, legend="bottom")

### Fragment Length Distribution

samtools view -F 0x04 alignment/sam/{samp_name}_bowtie2.sorted.dupMarked.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >alignment/sam/{samp_name}_fragmentLen.txt

fragLen = c()
for(hist in sampleList){

  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(paste0(projPath, "/alignment/sam/fragmentLen/", hist, "_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[1], Replicate = histInfo[2], sampleInfo = hist) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)

## Plot

fig5 = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("Fragment Length Distributions")

ggarrange(fig5)






