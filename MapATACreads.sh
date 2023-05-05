#!/bin/bash
#SBATCH --job-name=Lewislab_ATAC.%j.job
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zlewis@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapATAC.%j.out
#SBATCH --error=../ATAC.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

OUTDIR=../${OutputFolderName}
mkdir ${OUTDIR}


# #process reads using trimGalore
#
#  ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
#  trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
# #
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
  deduped="${OUTDIR}/SortedBamFiles/${name}_deduped.bam"
  bed="${OUTDIR}/Beds/${name}.bed"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#

ml SAMtools/1.9-GCC-8.3.0
ml BWA/0.7.17-GCC-8.3.0
#
bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"

#module load picard/2.27.4-Java-13.0.2

#java -jar picard.jar MarkDuplicates \
 #     I=${bam} \
  #    O=${deduped} \
   #   M=./${OUTDIR}/SortedBamFiles/marked_dup_metrics.txt \
    #  --REMOVE_DUPLICATES=TRUE
    
#perl ./shiftTn5_BAM_2_BED.pl "${bam}" > "${name}.bed"

############################
# # #deeptools
ml deepTools/3.3.1-intel-2019b-Python-3.7.4
alignmentSieve -p $THREADS --ATACshift --bam ${bam} -o ${name}.tmp.bam

# the bam file needs to be sorted again
samtools sort -@ $THREADS -O bam -o ${name}.shifted.bam ${name}.tmp.bam
samtools index -@ $THREADS ${name}.shifted.bam
rm ${name}.tmp.bam

#Plot all reads
bamCoverage -p $THREADS -bs 1 --normalizeUsing BPM --Offset 1 3 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes
#bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

done
