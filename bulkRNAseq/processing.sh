#!/usr/bin/bash

:<<"RNA_SEQ"

RNA_SEQ

__VERSION='1.1 - 07/2020'

## default values
C_THREADS=20          # compute threads, e.g. STAR aligner
ASSUME_STRANDED_RF=0  # (KAPA stranded mRNA / Illumina TruSeq)
ASSUME_STRANDED_FR=0
IS_TWO_COLOR_SEQ=0

if [ $# -eq 0 ]; then
  echo "Usage: $0 -1 FQ1 -2 FQ2 -s SAMPLE -g GENOME [ -a ANNO_FILE -c C_THREADS -n -x -y ]"
  echo ""
  echo "FQ1       : first fastq file"
  echo "FQ2       : second fastq file"
  echo "SAMPLE    : sample name, used for output naming"
  echo "GENOME    : Either mm9/mm10 or hg19/hg38. We are using GENCODE."
  echo "            So this is actually NCBIM37/GRCm38.p6 or GRCh37.p13/GRCh38.p12."
  echo "ANNO_FILE : annotation file used for STAR alignment; here for transcript assembly by stringtie"
  echo "            This is GTF/GFF3 formatted. def: genomes/STAR/current/\$GENOME/gencode.gtf"
  echo "C_THREADS : computing(alignment) threads to use."
  echo "n           Special flag for 2-color devices NextSeq and NovaSeq."
  echo ""
  echo "x         : Assume stranded library fr-firststrand  (KAPA stranded mRNA / Illumina TruSeq)"
  echo "y         : Assume stranded library fr-secondstrand (Others)"
  echo ""
  echo " Version : $__VERSION"
  exit
fi

## parameters
while getopts 1:2:s:g:a:c:nxy option
do
    case "${option}"
        in
        1) FQ1_INPUT=${OPTARG};;
        2) FQ2_INPUT=${OPTARG};;
        s) SAMPLE=${OPTARG};;
        g) GENOME=${OPTARG};;
        a) ANNO_INPUT=${OPTARG};;
        c) C_THREADS=${OPTARG};;
        n) IS_TWO_COLOR_SEQ=1;;
        x) ASSUME_STRANDED_RF=1;;
        y) ASSUME_STRANDED_FR=1;;
    esac
done

## safer way
set -e

## either option must be used, not both
if (( $ASSUME_STRANDED_RF && $ASSUME_STRANDED_FR )); then
    echo "[ERROR] You cannot specify -x and -y at the same time."
    exit
fi

## strandness of RNA library for stringtie
STRAND=''
if (( $ASSUME_STRANDED_RF )); then
    STRAND='--rf' # KAPA stranded mRNA
elif (( $ASSUME_STRANDED_FR )); then
    STRAND='--fr'
fi

## common stuff
working_dir=$(pwd)
date_flag=$(date +%Y%m%d)  #.%H%M%S)  # 20181025.144207
declare -a remove_me  # array for final cleanup

## avoid funny locales for core components QC, mapping, transcript assembly
export LC_ALL="POSIX"

## see https://github.com/alexdobin/STAR/issues/512
ulimit -n 4096

## programs to use
     fastqc=bin/fastqc
   cutadapt=bin/cutadapt
       star=bin/STAR
bamCoverage=bin/bamCoverage
   samtools=bin/samtools
  stringtie=bin/stringtie
    multiqc=bin/multiqc

## Effective Genome Sizes (for bamCoverage)
## https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
declare -A egs
# correct names
egs[grcm38]=2652783500

# common names
egs[mm10]=${egs[grcm38]}

## Genome Name Lookup
declare -A genome_index
# correct names
genome_index[grcm38]='genomes/STAR/current/GRCm38.p6'

# common names
genome_index[mm10]='genomes/STAR/current/mm10'

## Adaptors to trim -- todo: add parameter for different adaptors
declare -A adaptors
adaptors[illumina]='AGATCGGAAGAGC'
adaptors[smallrna]='TGGAATTCTCGG'
adaptors[nextera]='CTGTCTCTTATA'

## be verbose
function print_str {
  str=$1
  echo ""
  echo $(date +"### [%Y-%m-%d %H:%M:%S] $str")
}


## PARAMETER CHECK >>>>>>>>>>>>>>>
echo ""

## ----------------------------------------------------------------
## (c) Threads - some basic testing for positive integer
if [[ ! $C_THREADS =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] Not a valid number '$C_THREADS' for computing threads."
    exit
fi

# increase ulimit to allow for STAR to use more threads.
max_open_files=$(ulimit -n)
no_of_bins=50
max_threads=$(($max_open_files/$no_of_bins)) # bash:1024/50=20 float:20.48 - we want 20 :-)

if [ $C_THREADS -gt $max_threads ]; then
    echo "[ERROR] Too many threads used ($C_THREADS). Use max. $max_threads threads."
    exit
fi

## respects sth. like "cutadapt --threads N [..] | cutadapt --threads N [..]"
## as $C_THREADS corresponds to --memory=$C_THREADS for 'mxqsub'.
p_threads=$(($C_THREADS/2)) # processing threads, e.g. adaptor trimming

## ----------------------------------------------------------------
## (1,2) Fastq input data
if [ ! -r "$FQ1_INPUT" ] || [ ! -r "$FQ2_INPUT" ]; then
    echo "$FQ1_INPUT or $FQ2_INPUT cannot be read."
    exit
fi

## absolute path
FQ1=$(realpath -s $FQ1_INPUT)
FQ2=$(realpath -s $FQ2_INPUT)


## ----------------------------------------------------------------
## (g) Genome Reference
## returns empty string if $genome_index[$GENOME] is NOT set, "is_set_value" if set.
## --> bash parameter substitution
if [ -z ${genome_index[$GENOME]+"is_set_value"} ];then
    echo "[ERROR] Unknown reference genome '$GENOME'."
    exit
fi

## ----------------------------------------------------------------
## (a) Annotation used for transcript assembly.
## + if not provided, then we check index folder for "standard" gencode.gtf
## + to override use -a, but this must be readable
## + if -a set and file is readable, then go ahead
if [ -z $ANNO_INPUT ]; then
    if [ -r ${genome_index[$GENOME]}/gencode.gtf ]; then
        echo "[INFO] Using standard GTF file '${genome_index[$GENOME]}/gencode.gtf'."
        ANNO_FILE=${genome_index[$GENOME]}/gencode.gtf
    else
        "[ERROR] Cannot read annotation file for transcript assembly."
        exit
    fi
elif [ ! -r $ANNO_INPUT ]; then
    echo "[ERROR] Cannot read annotation file '$ANNO_INPUT' for transcript assembly."
    exit
else
    ANNO_FILE=$(realpath -s $ANNO_INPUT)
fi

## ----------------------------------------------------------------
## (s) Sample Name
if [ -z $SAMPLE ]; then
    echo "[ERROR] Please provide SAMPLE name."
    exit
fi

## PARAMETER CHECK <<<<<<<<<<<<<<<


## START WORKFLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ----------------------------------------------------------------
## (1) Prepare working/results directory
## ----------------------------------------------------------------
results_dir=${working_dir}/rnaseq_${SAMPLE}_${GENOME}_${date_flag}
fastqc_results_dir=$results_dir/00_FastQC
mkdir -p $fastqc_results_dir
cd $results_dir


## ----------------------------------------------------------------
## (2) pre-QC
## ----------------------------------------------------------------
print_str "FastQC : pre-QC"

$fastqc \
  --outdir $fastqc_results_dir \
  --threads $C_THREADS \
  --dir $fastqc_results_dir \
  $FQ1 $FQ2


## ----------------------------------------------------------------
## (3) Adaptor Trimming, incl. poly-A and quality
## ----------------------------------------------------------------
print_str "Running cutadapt - Trimming adaptors, polyAT and Quality"

min_length=25
min_qual=20

# output names, no absolute paths necessary
fq1_trimmed=${SAMPLE}_R1.trm.fq.gz
fq2_trimmed=${SAMPLE}_R2.trm.fq.gz

# If data has been generated on NextSeq or NovaSeq 2-color machines,
# we need to apply a different type of quality clipping (due to dark cycles)
trimming_type='--quality-cutoff'
if (( $IS_TWO_COLOR_SEQ )) ; then
    echo "[INFO] Applying 2-color-device trimming mode"
    trimming_type='--nextseq-trim'
fi

# quality clipping is done BEFORE adaptor clipping
$cutadapt \
  $trimming_type $min_qual \
  --overlap 5 \
  --minimum-length $min_length \
  --cores $(($p_threads+4)) \
  --adapter ${adaptors["illumina"]} \
  -A ${adaptors["illumina"]} \
  --interleaved \
  $FQ1 \
  $FQ2 \
  2>.log_a | \
$cutadapt \
  --interleaved \
  --overlap 20 \
  --minimum-length $min_length \
  --cores $(($p_threads+4)) \
  --adapter "A{100}" \
  --adapter "T{100}" \
  --output $fq1_trimmed \
  --paired-output $fq2_trimmed \
  - &>.log_b

cat .log_a .log_b > ${SAMPLE}_FastQ-Trimming.log

remove_me+=(".log_a" ".log_b" $fq1_trimmed $fq2_trimmed)


## ----------------------------------------------------------------
## (4) post-QC
## ----------------------------------------------------------------
print_str "FastQC : post-QC"

$fastqc \
  --outdir $fastqc_results_dir \
  --threads $C_THREADS \
  --dir $fastqc_results_dir \
  $fq1_trimmed $fq2_trimmed


## ----------------------------------------------------------------
## (5) Alignment - STAR
## ----------------------------------------------------------------
print_str "STAR Alignment : $SAMPLE on $GENOME (${genome_index[$GENOME]})"

bam_file=${SAMPLE}_${GENOME}.star.srt.bam
temp_dir=temp.STAR.${SAMPLE}_${GENOME}

$star \
  --genomeDir ${genome_index[$GENOME]} \
  --runMode alignReads \
  --runThreadN $C_THREADS \
  --readFilesCommand gzip -dc \
  --outFileNamePrefix ${SAMPLE}_${GENOME}_ \
  --outTmpDir $temp_dir \
  --outStd Log \
  --outSAMtype BAM SortedByCoordinate \
  --chimSegmentMin 20 \
  --outSAMstrandField intronMotif \
  --quantMode GeneCounts \
  --readFilesIn $fq1_trimmed $fq2_trimmed

## file name given by STAR
bam_star=${SAMPLE}_${GENOME}_Aligned.sortedByCoord.out.bam

## be sure that BAM is really completely written
while ! $samtools quickcheck $bam_star; do
    print_str "BAM file ($bam_star) not yet complete. Waiting ..."
    sleep 2m
done

print_str "Renaming Standard BAM to $bam_file"
mv $bam_star $bam_file

print_str "Creating BAM Index for $bam_file"
$samtools index -@ $((2*$p_threads)) $bam_file

## rename output files consistently
for postfix in 'Log.final.out' 'Log.out' 'Log.progress.out' 'ReadsPerGene.out.tab' 'SJ.out.tab'; do
    if [ -e ${SAMPLE}_${GENOME}_$postfix ]; then
        mv -v ${SAMPLE}_${GENOME}_$postfix ${SAMPLE}_${GENOME}.star.$postfix
    fi
done

remove_me+=($temp_dir)


## ----------------------------------------------------------------
## (6) Transcript Assembly - Stringtie
## ----------------------------------------------------------------
print_str "Transcript Assembly - Stringtie"

assembled_transcripts=${SAMPLE}_${GENOME}.stringtie.transcripts.gtf
       gene_abundance=${SAMPLE}_${GENOME}.stringtie.gene_abundance.txt
        coverage_file=${SAMPLE}_${GENOME}.stringtie.coverage.gff3

ballgown_tables_dir=10_Ballgown_Tables
mkdir -p $ballgown_tables_dir

## respect strandness
$stringtie \
  $bam_file \
  $STRAND \
  -G $ANNO_FILE \
  -o $assembled_transcripts \
  -C $coverage_file \
  -p $C_THREADS \
  -A $gene_abundance \
  -b $ballgown_tables_dir \
  -e \
  -v


## ----------------------------------------------------------------
## (7) Create BigWig Coverage File
## ----------------------------------------------------------------
## RPKM = Reads Per Kilobase per Million mapped reads
##  CPM = Counts Per Million mapped reads, same as CPM in RNA-seq
##  BPM = Bins Per Million mapped reads, same as TPM in RNA-seq
## RPGC = reads per genomic content (1x normalization)
print_str "bamCoverage (deeptools): compute coverage (bigwig) using Effective Genome Size '$GENOME' = ${egs[$GENOME]}"

export LC_ALL=de_DE.UTF-8

for n in 'None' 'RPKM'; do
    $bamCoverage \
      --bam $bam_file \
      --normalizeUsing $n \
      --outFileName ${bam_file}.cov.${n}.bw \
      --numberOfProcessors $C_THREADS \
      --effectiveGenomeSize ${egs[$GENOME]}

    ## strand-specific coverage plot
    for f in 'forward' 'reverse'; do
        $bamCoverage \
          --bam $bam_file \
          --normalizeUsing $n \
          --filterRNAstrand  $f \
          --outFileName ${bam_file}.cov.${n}.${f:0:3}.bw \
          --numberOfProcessors $C_THREADS \
          --effectiveGenomeSize ${egs[$GENOME]}
    done
done


## ----------------------------------------------------------------
## (8) Cleaning up
## ----------------------------------------------------------------
print_str "Cleaning up"

for item in "${remove_me[@]}"
do
    echo " removing '$item'"
    rm -fr $item
done
