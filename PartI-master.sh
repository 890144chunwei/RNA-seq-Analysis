#!/usr/bin/bash
#SBATCH --mem 32768
#SBATCH -p short
#SBATCH --output=/storage/goodell/projects/chunweic/slurm_out/220601_rnaseq_trx_%j.out
#SBATCH -e /storage/goodell/projects/chunweic/slurm_out/220601_rnaseq_trx_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=chunweic@bcm.edu # Email to which notifications will be sent

pwd; hostname; date

BCLDIR="/storage/goodell/bcl/220601_Chunwei_RNA_Trx/Files"
GENOMDIR="/storage/goodell/home/chunweic/mm10/STARgenome"
FASTQDIR="/storage/goodell/fastq/220601_Chunwei_RNA_Trx"
PROJECTDIR="/storage/goodell/projects/chunweic/220601_Chunwei_RNA_Trx"

bcl2fastq -R $BCLDIR -i "$BCLDIR/Data/Intensities/BaseCalls/" -o $FASTQDIR --sample-sheet "$PROJECTDIR/SampleSheet_220601.csv" --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --loading-threads 8 -p 16 --tiles s_[1]

cp -r "$FASTQDIR/220601" $PROJECTDIR
fastqc $PROJECTDIR/*gz
gunzip $PROJECTDIR/*gz
multiqc $PROJECTDIR/*zip -o $PROJECTDIR/

for FILE in $PROJECTDIR/*R1.fastq ;
do
  STAR --genomeDir $GENOMDIR/ --runMode alignReads --runThreadN 10 --readFilesIn $PROJECTDIR/${FILE%_R1.fastq}_R1.fastq $PROJECTDIR/${FILE%_R1.fastq}_R2.fastq \
  --outFileNamePrefix $FASTQDIR/$FILE_ --outSAMtype BAM Unsorted --outSAMunmapped Within \
  --limitGenomeGenerateRAM 1700000000000;
  samtools sort $PROJECTDIR/$FILE_Aligned.out.bam $PROJECTDIR/$FILE
  featureCounts -T 7 -p -s 2 -t exon -a $GENOMDIR/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf -o $PROJECTDIR/$FILE_fc.txt $PROJECTDIR/$FILE_sort.bam
  samtools index $PROJECTDIR/$FILE_sort.bam
  bamCoverage -b $PROJECTDIR/$FILE_sort.bam -o $PROJECTDIR/Coverage/$FILE_cov.bw --normalizeUsing RPKM ;
done

rm $PROJECTDIR/*fastq $$PROJECTDIR/*zip $$PROJECTDIR/*out.bam
