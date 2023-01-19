#!/usr/bin/bash
#SBATCH --mem 32768
#SBATCH -p short
#SBATCH --output=/storage/goodell/projects/chunweic/slurm_out/220601_rnaseq_trx_%j.out
#SBATCH -e /storage/goodell/projects/chunweic/slurm_out/220601_rnaseq_trx_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=chunweic@bcm.edu # Email to which notifications will be sent

pwd; hostname; date

GENOMDIR="/storage/goodell/home/chunweic/mm10/STARgenome"
FASTQDIR="/storage/goodell/projects/chunweic/220601_Chunwei_RNA_Trx/"

STAR --genomeDir $GENOMDIR/ --runMode alignReads --runThreadN 10 --readFilesIn $FASTQDIR/HSPC_germline_WT1_R1.fastq $FASTQDIR/HSPC_germline_WT1_R2.fastq \
--outFileNamePrefix $FASTQDIR/HSPC_germline_WT1 --outSAMtype BAM Unsorted --outSAMunmapped Within \
--limitGenomeGenerateRAM 1700000000000;

samtools sort $FASTQDIR/Tcell_WT3Aligned.out.bam $FASTQDIR/HSPC_germline_WT1

featureCounts -T 7 -p -s 2 -t exon -a $GENOMDIR/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf -o $FASTQDIR/HSPC_germline_WT1_fc.txt \
$FASTQDIR/HSPC_germline_WT1_sort.bam

samtools index $FASTQDIR/HSPC_germline_WT1_sort.bam
bamCoverage -b $FASTQDIR/HSPC_germline_WT1_sort.bam -bl $HOMEDIR/mm10-blacklist.v2.bed -o $FASTQDIR/Coverage/HSPC_germline_WT1_cov.bw --normalizeUsing RPKM
