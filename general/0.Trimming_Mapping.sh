# Submit as a sample-based positional script 
RUN=$1

# Define the working directories and output directory
wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/WGS/ILLUMINA/
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq

# File input variables defined by the RUN variable
file1=$(ls ${wd}/${RUN}__R1__*)
file2=$(ls ${wd}/${RUN}__R2__*)

# Create a temporary directory in /tmp for the current job
SCRATCH=/tmp/$SLURM_JOB_ID

# Perform adapter trimming using bbduk.sh on the input files
bbduk.sh t=6 -Xmx24g in=${file1} in2=${file2} out=${outdir}/${RUN}.trim.fastq.gz minlen=25 qtrim=rl trimq=2 ktrim=r k=23 mink=11 ref=~/modules/bbmap/adapters.fa hdist=1 tpe tbo

# Index the genome using BWA
bwa index ${genome}

# Define RG tags based on metadata files specific to the current RUN
ID="$(cat ${wd}/RGs/${RUN}.ID | tail -1)"
PU="$(cat ${wd}/RGs/${RUN}.PU | tail -1)"
SM="$(cat ${wd}/RGs/${RUN}.SM | tail -1)"
LB="$(cat ${wd}/RGs/${RUN}.LB | tail -1)"

# Perform alignment with BWA MEM, including RG information, then sort with samtools, write output to temp directory
bwa mem -M -p -t 10 -R "@RG\tID:${ID}\tSM:${SM}\tPL:ILLUMINA\tPU:${PU}\tLB:${LB}" ${genome} ${qcdata}/${RUN}.trim.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -

# Filter the alignment for properly paired reads, write the output to the final output directory
samtools view -b -f 1 -F 524 ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam

# Index the final bam file
samtools index -b ${outdir}/${RUN}.bam

# Make output directories for final bams and stats
mkdir tmp
mkdir $SCRATCH
mkdir ${bamdir}/merged
mkdir ${bamdir}/merged/stats
mkdir ${bamdir}/merged/stats/alignments
mkdir ${bamdir}/merged/stats/coverage

# Merge all the individual library bams for the current RUN into one bam, sort, and index the merged bam
samtools merge - ${bamdir}/${RUN}*.bam | samtools sort -@3 -o $SCRATCH/${RUN}.bam
samtools index $SCRATCH/${RUN}.bam

# Mark duplicates with sambamba, and index the result
sambamba markdup --tmpdir $SCRATCH/ $SCRATCH/${RUN}.bam $SCRATCH/${RUN}.dup.bam
samtools index $SCRATCH/${RUN}.dup.bam

# Remove reads extending beyond scaffold ends using GATK's CleanSam
gatk CleanSam -I $SCRATCH/${RUN}.dup.bam -O ${bamdir}/merged/${RUN}.bam -R ${genome} --CREATE_INDEX true

# Collect alignment stats with GATK's CollectAlignmentSummaryMetrics
gatk CollectAlignmentSummaryMetrics -R ${genome} -I ${bamdir}/merged/${RUN}.bam -O ${bamdir}/merged/stats/alignments/${RUN}.alignments.txt -LEVEL SAMPLE -LEVEL READ_GROUP

# Calculate coverage in 100kb windows, using the median value 
mosdepth --threads 3 --use-median --by 100000 --fast-mode --no-per-base ${bamdir}/merged/stats/coverage/${RUN} ${bamdir}/merged/${RUN}.bam
