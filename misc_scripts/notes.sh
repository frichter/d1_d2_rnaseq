

## download mm13

cd /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq

## run fastqc
## link: https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
module load fastqc/0.11.5
## generate an html report with symlink to user directory
time fastqc *fastq.gz


# M. musculus, GRCm38
# from https://ccb.jhu.edu/software/hisat2/index.shtml

## align with HISAT2 and STAR
## manual https://ccb.jhu.edu/software/hisat2/manual.shtml
module load hisat2/2.0.5

## what is the --dta tag? did I use this? Need this for stringtie
hisat2-build
hisat2 --time -x grcm38_snp_tran/genome_snp_tran -1 D1_D2_odd_repeat/D1_CTRL1_ATCACG_L001_R1.fastq.gz -2 D1_D2_odd_repeat/D1_CTRL1_ATCACG_L001_R2.fastq.gz -S D1_D2_odd_repeat/D1_CTRL1.sam
##--time prints walltime

## why isn't there a unction file?

## do I need to trim 5' or 3'?
## --phred33 or --phred64? do I have to specify?

## confirm that you are measuring retained introns: seems like it since gene feature
## in feature counts is more than 2x exon counts
--known-splicesite-infile, --novel-splicesite-outfile, --novel-splicesite-infile
## how long are the reads?

## is there a default for --un <path> --al <path> --un-conc <path> or --al-conc <path>?
--met-file <path> --met-stderr

## add readgroup IDs
--rg-id <text>

--p

HISAT2_PATH="/hpc/packages/minerva-common/hisat2/2.0.5/src/hisat2-2.0.5"
head $HISAT2_PATH/hisat2_extract_snps_haplotypes_UCSC.py

## according to HISAT2 website don't need to run make_grcm38_snp_tran.sh
## since " You can use the script to build the same HISAT2 index we provide."
## https://ccb.jhu.edu/software/hisat2/indexes.txt

## not sure which of these:
module load star/2.5.3a
module load rna-star/2.4.0d
## possibly use rnacocktail on local machine


### download mouse reference GTF
wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz

### de novo gtf with stringtie
module load samtools/1.4.1
module load stringtie/1.3.3b

## not ref guided
# stringtie D1_D2_odd_repeat/D1_CTRL1.sam -v -o D1_D2_odd_repeat/D1_CTRL1.gtf -p 12

## ref guided
time samtools view -u D1_D2_odd_repeat/D1_CTRL1.sam | samtools sort -T D1_D2_odd_repeat/D1_CTRL1.temp -o D1_D2_odd_repeat/D1_CTRL1.sorted.bam -@ 12
## 3 minutes
time stringtie D1_D2_odd_repeat/D1_CTRL1.sorted.bam -v -o D1_D2_odd_repeat/D1_CTRL1_refguided.gtf -p 12 \
-G Mus_musculus.GRCm38.90.gtf

# -p number of threads
# consider -A <gene_abund.tab> Gene abundances will be reported (tab
# delimited format) in the output file with the given name.
# -C just provides

## after generating gtf for single file, merge with "--merge" option
# https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual


## there are multiple D1_D2_whole_cell/D2CTRL3_
D2CTRL3_ACTTGA_L005, D2CTRL3_ACTTGA_L008


ln -s /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq .

##
