#!/usr/bin/bash
#SBATCH --mem 32768
#SBATCH -p short
#SBATCH --output=/storage/goodell/projects/chunweic/slurm_out/220601_rnaseq_trx_%j.out
#SBATCH -e /storage/goodell/projects/chunweic/slurm_out/220601_rnaseq_trx_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=chunweic@bcm.edu # Email to which notifications will be sent

pwd; hostname; date

BCLDIR=
GENOMDIR="/storage/goodell/home/chunweic/mm10/STARgenome"
FASTQDIR="/storage/goodell/projects/chunweic/220601_Chunwei_RNA_Trx"

for FILE in $FASTQDIR/*R1.fastq ;
do
  STAR --genomeDir $GENOMDIR/ --runMode alignReads --runThreadN 10 --readFilesIn $FASTQDIR/${FILE%_R1.fastq}_R1.fastq $FASTQDIR/${FILE%_R1.fastq}_R2.fastq \
  --outFileNamePrefix $FASTQDIR/$FILE_ --outSAMtype BAM Unsorted --outSAMunmapped Within \
  --limitGenomeGenerateRAM 1700000000000;
  samtools sort $FASTQDIR/$FILE_Aligned.out.bam $FASTQDIR/$FILE
  featureCounts -T 7 -p -s 2 -t exon -a $GENOMDIR/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf -o $FASTQDIR/$FILE_fc.txt $FASTQDIR/$FILE_sort.bam
  samtools index $FASTQDIR/$FILE_sort.bam
  bamCoverage -b $FASTQDIR/$FILE_sort.bam -bl $HOMEDIR/mm10-blacklist.v2.bed -o $FASTQDIR/Coverage/$FILE_cov.bw --normalizeUsing RPKM ;
done
