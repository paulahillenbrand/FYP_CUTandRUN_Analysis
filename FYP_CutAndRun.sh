*** LINUX COMMAND LINE ***

### Set Up Conda Environment
conda create --name cutandrun
conda activate cutandrun

# Packages from Analysis Pipeline
conda install bioconda::fastqc=0.11.9
conda install bioconda::bowtie2=2.3.4.3
conda install bioconda::samtools=1.10
conda install bioconda::bedtools=2.29.1
conda install bioconda::picard=2.18.29
conda install bioconda::seacr=1.3
conda install bioconda::deeptools=2.0 

# Additional Packages
conda install bioconda::multiqc=1.0
conda install bioconda::trim-galore=0.6.10
conda install bioconda::sambamba=1.0.1


### Adapter trimming
trim_galore \
--paired \
--o trimmed \
rawdata/{samp_name}1.fq.gz \
rawdata/{samp_name}2.fq.gz

### Alignment
bowtie2 \
--end-to-end \
--very-sensitive \
--no-mixed \
--no-discordant \
--phred33 \
-I 10 -X 700 \
-x genome/HSymV2.1 \
-1 trimmed/{samp_name}_1.fq.gz \
-2 trimmed/{samp_name}_2.fq.gz \
-S alignment/sam/{samp_name}_bowtie2.sam \
&> alignment/sam/{samp_name}_bowtie2.txt


### Mark Duplicates

picard SortSam \
I={samp_name}_bowtie2.sam \
O={samp_name}_bowtie2.sorted.sam \
SORT_ORDER=coordinate

picard MarkDuplicates \
I={samp_name}_bowtie2.sorted.sam \
O={samp_name}_bowtie2.sorted.dupMarked.sam \
METRICS_FILE={samp_name}_picard.dupMark.txt

### File Format Conversion

## sam to bam 
samtools view \
-Sb \
-F 0x04 \
alignment/sam/{samp_name}_bowtie2.sorted.dupMarked.sam \
-o alignment/bam/{samp_name}_bowtie2.bam

## bam to bed 
bedtools bamtobed -i {samp_name}.bed


### Mark Multimappers
sambamba view \
-f bam \
-F "[XS] == null and not unmapped" \
alignment/bam/{samp_name}_bowtie2.bam \
-o alignment/bam/{samp_name}_bowtie2_nomultimap.bam

### File Format Conversion
samtools sort \
/cutandrun/alignment/bam/no_multimappers/{samp_name}_bowtie2_nomultimap. bam \
-o /cutandrun/alignment/bam/no_multimappers/{samp_name}_bowtie2_nomultimap.sorted.bam

for file in /alignment/bam/no_multimappers/*.bam; do samtools index $file; done

## bedgraph
bamCoverage \
-b /cutandrun/alignment/bam/no_multimappers/{samp_name}_bowtie2_nomultimap.sorted.bam \
 -o /cutandrun/alignment/bedgraph/{samp_name}_bowtie2_nomultimap.sorted.norm.bedgraph \
–outFileFormat bedgraph\
--normalizeUsing RPKM \
--exactScaling \
--effectiveGenomeSize 185312851

## bigwig
bamCoverage \
-b /cutandrun/alignment/bam/no_multimappers/{samp_name}_bowtie2_nomultimap.sorted.bam \
 -o /cutandrun/alignment/bam/no_multimappers/normalised_bigwig/{samp_name}_bowtie2_nomultimap.sorted.norm.bw \
–outFileFormat bigwig\
--normalizeUsing RPKM \
--exactScaling \
--effectiveGenomeSize 185312851

### Peak Calling
seacr="/Users/paula/opt/anaconda3/pkgs/seacr-1.3-hdfd78af_2/bin/SEACR_1.3.sh" \

## in AP2 samples
bash $seacr \ 
/cutandrun/alignment/bedgraph/{AP_tech_rep1}_bowtie2_nomultimap.sorted.norm.bedgraph \
/cutandrun/alignment/bedgraph/{Pre_rep1}_bowtie2_nomultimap.sorted.norm.bedgraph \
non stringent \
/cutandrun/alignment/bedgraph/peakCalling/{AP_tech_rep1}_ctrls.peaks

## in H3K27me3 samples
bash $seacr \ 
/cutandrun/alignment/bedgraph/{K27_tech_rep1}_bowtie2_nomultimap.sorted.norm.bedgraph \
0.01 \
non stringent \
/cutandrun/alignment/bedgraph/peakCalling/{AP_tech_rep1}_ctrls.peaks

## Consensus Peaks
bedtools intersect \
-a /cutandrun/alignment/bedgraph/peaks/{AP_tech_rep_Pre1}_ctrl_1_2.peaks.stringent.bed \
-b /cutandrun/alignment/bedgraph/peaks/{AP_tech_rep_Pre2}_ctrl.peaks.stringent.bed \
> /cutandrun/alignment/bedgraph/peaks/{ AP_techrep}_intersect.bed

# count lines
wc -l { AP_techrep}_intersect.bed


### PyGenomeTracks visualisation of Alignment Results in the Genomic Context
## make suitable conda environment
conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks python=3.7

## load tracks.ini file (including tracks of interest) 
pyGenomeTracks --tracks tracks.ini --region {genomic_region} -o image.png 


### Heatmap of H3K27me3 Distribution over transcriptional units
computeMatrix scale-regions \
-S /alignment/bam/no_multimappers/normalised_bigwig/{K27_rep1}_bowtie2_nomultimap.sorted.norm.bw /alignment/bam/no_multimappers/normalised_bigwig/{K27_rep2}_bowtie2_nomultimap.sorted.norm.bw \
-R /Genome/HSymV2.1.gtf \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--skipZeros \
-o /Users/paula/Downloads/cutandrun/matrix_gene.mat.gz  \

plotHeatmap \
-m /matrix_gene.mat.gz \
-out /K27_gene.png \
--sortUsing sum
--samplesLabel "{samp_name}" "{samp_name}"

### Heatmap of peak shapes
Peak shape
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' \
/alignment/bedgraph/peaks/{samp_name}_seacr_ctrl.peaks.stringent.bed \
> /alignment/bedgraph/peaks/{samp_name}_seacr_ctrl.summitRegion.bed

computeMatrix reference-point \
-scoreFileName /alignment/bam/no_multimappers/normalised_bigwig/{samp_name}_bowtie2_nomultimap.sorted.norm.bw \
-regionsFileName /alignment/bedgraph/peaks/{samp_name}_seacr_ctrl.summitRegion.bed \
--skipZeros \
-o /alignment/bedgraph/peaks/{samp_name}_SEACR.mat.gz \
--beforeRegionStartLength 3000 \
--afterRegionStartLength 3000 \
--referencePoint center

plotHeatmap \
-m /alignment/bedgraph/peaks/{samp_name}_SEACR.mat.gz \
-out /alignment/bedgraph/peaks/{samp_name}_SEACR_heatmap.png \
--sortUsing sum \
--startLabel "Peak Start" \
--endLabel "Peak End" \
--regionsLabel "Peaks" \
--samplesLabel "{samp_name}"