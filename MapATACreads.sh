#!/bin/bash
#SBATCH --job-name=Lewislab_ATAC.%j.job
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zlewis@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=./MapATAC.%j.out
#SBATCH --error=./MapATAC.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

OUTDIR=../${OutputFolderName}
mkdir ${OUTDIR}


# #process reads using trimGalore
#
module load Trim_Galore/0.6.10-GCCcore-12.3.0
trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
#
FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#
 mkdir "${OUTDIR}/SortedBamFiles"
 mkdir "${OUTDIR}/BigWigs"
 mkdir "${OUTDIR}/Peaks"
 mkdir "${OUTDIR}/Beds"
#mkdir "$OUTDIR/HomerTagDirectories"
#mkdir "$OUTDIR/TdfFiles"
#
#Iterate over the files
for f in $FILES
do
#
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
# 		#dir=${f%/*}
		
	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	name=${file/%_S[1-12]*_L001_R1_001_val_1.fq.gz/}

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam file
 	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	shifted="${OUTDIR}/SortedBamFiles/${name}.shifted.bam"
  	deduped="${OUTDIR}/SortedBamFiles/${name}_deduped.bam"
  	bed="${OUTDIR}/Beds/${name}.bed"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#
module load  SAMtools/1.21-GCC-13.3.0
module load  BWA/0.7.18-GCCcore-13.3.0

#
bwa mem -M -v 3 -a -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"

#module load picard/2.27.4-Java-13.0.2

#java -jar picard.jar MarkDuplicates \
 #     I=${bam} \
  #    O=${deduped} \
   #   M=./${OUTDIR}/SortedBamFiles/marked_dup_metrics.txt \
    #  --REMOVE_DUPLICATES=TRUE
    
#perl ./shiftTn5_BAM_2_BED.pl "${bam}" > "${name}.bed"

############################
# #deeptools
module load deepTools/3.5.5-gfbf-2023a
alignmentSieve -p $THREADS --ATACshift --bam ${bam} -o ${name}.tmp.bam

# the bam file needs to be sorted again
samtools sort -@ $THREADS -O bam -o ${shifted} ${name}.tmp.bam
samtools index -@ $THREADS ${shifted}
rm ${name}.tmp.bam

#Plot all reads
bamCoverage -p $THREADS --Offset 1 1 -bs 1 --smoothLength 6 --normalizeUsing BPM  -of bigwig -b ${shifted} -o "${bigwig}.ATAC_bin_1smooth_1_Bulk.bw"

#plot mononucleosomes
#bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

done
