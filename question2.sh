#!/bin/sh
##### variant calling Pipeline for Clinical reporting#####
echo "Welcome to Clinical Reporting Pipeline"
NUMCPUS=8

####The path for the following program must be set in the server ###
#### CHECK for Program paths ####
FAILURE=""
SAMTOOLS=$(which samtools)
if [ ! -x $SAMTOOLS ]
then
          FAILURE=$SAMTOOLS
fi

BCFTOOLS=$(which bcftools)
if [ ! -x $BCFTOOLS ]
then
	FAILURE=$BCFTOOLS
fi

BWA=$(which bwa)
if [ ! -x $BWA ]
then
	FAILURE=$BWA
fi

GATK="$(which gatk)"
if [ ! -x $GATK ]
then
	FAILURE=$GATK
fi

if [ "$FAILURE" ]
then
  echo "ERROR: $FAILURE program not found, please edit the configuration script."
  exit 1
else
echo "SUCCESS: All program found. proceeding for analysis."
fi

### BASE DIRECTORY FOR DATABASES ###
#PATH={Home Directory}
BASEDIRDB=$PATH/Human_genome/hg38
DBSNP="$PATH/dbsnp154_hg38_p13.vcf"
GENOMEIDX1="$PATH/hg38.p13_23chr.fna"

TEMPLOC="./tmp" #this will be relative to the output directory

### READ SAMPLE LIST ###
#BASEDIRDATA="$(pwd)"
reads1=(${BASEDIRDATA}/*_1P.fastq)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1P./_2P.}")

if [ $reads1 -a $reads2 ]
then
echo "***Following SAMPLES are submitted for Germline variant calling***"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"
    echo "$sample"
    done
else
 echo "ERROR: $reads1 $reads2 program not found, please edit the configuration script."
  exit 1
fi
### GATK haplotype caller ###
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"
    stime=`date +"%Y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $sample"
   
## fastp and trimmomatic software may be used for read quality control

    echo [$stime] "   * Alignment of reads to genome (BWA)"
   $BWA mem -t $NUMCPUS ${GENOMEIDX1} \
              ${BASEDIRDATA}/${reads1[$i]} \
             ${BASEDIRDATA}/${reads2[$i]} \
 		-o ${sample}.sam 
    
  echo [$stime] "   * Convertion of SAM to BAM and sorting"
 $GATK --java-options "-Xmx24G" SortSam \
      	--INPUT ${sample}.sam     \
     	--OUTPUT ${sample}_sort.bam --SORT_ORDER coordinate

##Alignment QC##
    
  echo [$stime] "   * Alignment summary"
 $GATK --java-options "-Xmx24G" CollectAlignmentSummaryMetrics \
	       --REFERENCE_SEQUENCE ${GENOMEIDX1} \
             --INPUT ${sample}_sort.bam --OUTPUT ${sample}_alignment_metrics.txt --VALIDATION_STRINGENCY LENIENT

  echo [$stime] "   * Insert size information"
 $GATK --java-options "-Xmx24G" CollectInsertSizeMetrics \
       --INPUT ${sample}_sort.bam 		   \
	       --OUTPUT ${sample}_insert_metrics.txt --Histogram_FILE ${sample}_histogram.pdf

   echo [$stime] "   * Insert size information"
    $GATK --java-options "-Xmx24G" AddOrReplaceReadGroups \
	       --INPUT ${sample}_sort.bam \
	       --OUTPUT ${sample}_readgroup.bam --RGLB lib1 --RGPL illumina --RGPU NONE --RGSM BG

   echo [$stime] "   * MarkDuplicates"
   $GATK --java-options "-Xmx24G" MarkDuplicates \
  	       --INPUT ${sample}_readgroup.bam \
           --OUTPUT ${sample}_dedup_reads.bam \
          --METRICS_FILE metrics.txt -AS true --VALIDATION_STRINGENCY LENIENT

#### BAM file Pre-processing and Variant calling ####
 echo [$stime] "   * BuildBamIndex"
$GATK --java-options "-Xmx24G" BuildBamIndex --INPUT ${sample}_dedup_reads.bam

 echo [$stime] "   * RealignerTargetCreator"
$GATK --java-options "-Xmx24G" BaseRecalibrator \
     -R ${GENOMEIDX1} -I ${sample}_dedup_reads.bam --known-sites ${DBSNP} -O ${sample}_recal_data.table

$GATK --java-options "-Xmx24G" ApplyBQSR -R ${GENOMEIDX1} \
  -I ${sample}_dedup_reads.bam \
  -bqsr ${sample}_recal_data.table \
         -O ${sample}_recal_reads.bam

  $GATK --java-options "-Xmx24G" BaseRecalibrator -R ${GENOMEIDX1} \
       -I ${sample}_recal_reads.bam \
      --known-sites ${DBSNP} -O ${sample}_post_recal_data.table

 $GATK --java-options "-Xmx24G" AnalyzeCovariates -R ${GENOMEIDX1} \
      -before ${sample}_recal_data.table \
     -after ${sample}_post_recal_data.table -plots ${sample}_recalibration_plots.pdf

 echo [$stime] "   *Germline variant calling - HaplotypeCaller"
$GATK --java-options "-Xmx24G" HaplotypeCaller --dbsnp ${DBSNP} -R ${GENOMEIDX1} \
     -I ${sample}_recal_reads.bam -O ${sample}.vcf

#### Variant Filter ####
#### GATK Hard Fillter ###

cat ${sample}.vcf | java -jar $PATH/SnpSift.jar filter "((QUAL >= 30 )&&(GEN[*].AD[1] >= 10)&&(GEN[*].AD[1] > GEN[*].AD[0]) &&(QD >= 2.0)&&((FS<=60)|(SOR <= 3.0))&&(MQ>=40)&&((! exists MQRankSum) |(MQRankSum >= -12.5))&&((! exists ReadPosRankSum) | (ReadPosRankSum >= -8.0)))" > ${sample}_filter.vcf

#### Variant annotation #######

java -Xmx64G -jar $PATH/snpEff.jar -v GRCh37.p13 ${sample}_filter.vcf  > ${sample}_ann.vcf

#### Filtering variants with minor allele frequency < 1% in normal population. 1000 Genome database is used for this purpose in this program
### Other population databases may also be used

java -Xmx64G -jar $PATH/SnpSift.jar annotate$PATH/hg38_1000GENOMES-phase_3_ncid.vcf.gz ${sample}_ann.vcf > ${sample}_1000G.vcf

cat ${sample}_1000G.vcf   |  java -Xmx64g -jar $PATH/SnpSift.jar filter "((! exists MAF) | (MAF <= 0.01))" > ${sample}_MAF01.vcf

### Mapping with a clinvar database for reporting pathogenic or likely pathogenic variants associated with a disease. 
## Any other in-house or proprietary databases may also be used.

java -Xmx64G -jar $PATH/SnpSift.jar annotate $PATH/GRCh38_latest_clinvar_ncid.vcf.gz ${sample}_MAF01.vcf > ${sample}_clinvar.vcf

# Predicting variant impact on gene and classify them based on ACMG guideline
$PATH/annovar/table_annovar.pl ${sample}_chno.vcf humandb/ -buildver hg38 -out ${sample}_myanno -remove -protocol refGene,cytoBand,dbnsfp41a,intervar_20180118 -operation g,r,f,f -nastring . -vcfinput -polish

java -Xmx64G -jar $PATH/SnpSift.jar extractFields ${sample}_myanno.vcf "CHROM" "POS" "REF" "ALT" "QUAL" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].FEATURE" "ANN[0].BIOTYPE" "ANN[0].AA" "CLNDN" "CLNSIG" "MAF" "CLIN_likely_pathogenic" "COSMIC_90" "CLIN_pathogenic" "DamagePredCount" "SIFT_pred" "SIFT4G_pred" "Polyphen2_HDIV_pred" "Polyphen2_HVAR_pred" "LRT_pred" "MutationTaster_pred" "MutationAssessor_pred" "FATHMM_pred" "PROVEAN_pred" "Interpro_domain" "InterVar_automated" "PVS1" "PS1" "PS2" "PS3" "PS4" "PM1" "PM2" "PM3" "PM4" "PM5" "PM6" "PP1" "PP2" "PP3" "PP4" "PP5" "BA1" "BS1" "BS2" "BS3" "BS4" "BP1" "BP2" "BP3" "BP4" "BP5" "BP6" "BP7" "GEN[*].GT" > all_filtered_variant_annot.txt
 
##### The report would be generated as text file and the same may be imported into SQL database generate a clinical report in PDF #####

done

