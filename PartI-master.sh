#!/usr/bin/bash
#SBATCH --mem 32768
#SBATCH -p short
#SBATCH --output=/storage/goodell/projects/chunweic/slurm_out/220601_rnaseq_trx_%j.out
#SBATCH -e /storage/goodell/projects/chunweic/slurm_out/220601_rnaseq_trx_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=chunweic@bcm.edu # Email to which notifications will be sent

pwd; hostname; date

GENOMDIR="/storage/goodell/home/chunweic/mm10/STARgenome"
FASTQDIR="/storage/goodell/projects/chunweic/220601_Chunwei_RNA_Trx2wk/"

STAR --genomeDir $GENOMDIR/ --runMode alignReads --runThreadN 10 --readFilesIn $FASTQDIR/Tcell_WT3_R1.fastq $FASTQDIR/Tcell_WT3_R2.fastq \
--outFileNamePrefix $FASTQDIR/Tcell_WT3 --outSAMtype BAM Unsorted --outSAMunmapped Within \
--limitGenomeGenerateRAM 1700000000000;

samtools sort $FASTQDIR/Tcell_WT3Aligned.out.bam $FASTQDIR/Tcell_WT3_sort

featureCounts -T 7 -p -s 2 -t exon -a $GENOMDIR/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf -o $FASTQDIR/Tcell_WT3_fc.txt \
$FASTQDIR/Tcell_WT3_sort.bam

STAR --genomeDir $GENOMDIR/ --runMode alignReads --runThreadN 10 --readFilesIn $FASTQDIR/Tcell_WT1_R1.fastq $FASTQDIR/Tcell_WT1_R2.fastq \
--outFileNamePrefix $FASTQDIR/Tcell_WT1 --outSAMtype BAM Unsorted --outSAMunmapped Within \
--limitGenomeGenerateRAM 1700000000000;

samtools sort $FASTQDIR/Tcell_WT1Aligned.out.bam $FASTQDIR/Tcell_WT1_sort

featureCounts -T 7 -p -s 2 -t exon -a $GENOMDIR/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf -o $FASTQDIR/Tcell_WT1_fc.txt \
$FASTQDIR/Tcell_WT1_sort.bam

STAR --genomeDir $GENOMDIR/ --runMode alignReads --runThreadN 10 --readFilesIn $FASTQDIR/Tcell_WT2_R1.fastq $FASTQDIR/Tcell_WT2_R2.fastq \
--outFileNamePrefix $FASTQDIR/Tcell_WT2 --outSAMtype BAM Unsorted --outSAMunmapped Within \
--limitGenomeGenerateRAM 1700000000000;

samtools sort $FASTQDIR/Tcell_WT2Aligned.out.bam $FASTQDIR/Tcell_WT2_sort

featureCounts -T 7 -p -s 2 -t exon -a $GENOMDIR/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf -o $FASTQDIR/Tcell_WT2_fc.txt \
$FASTQDIR/Tcell_WT2_sort.bam
