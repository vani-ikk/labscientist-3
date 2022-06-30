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
