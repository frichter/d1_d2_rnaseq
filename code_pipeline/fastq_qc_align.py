# Felix Richter
# felix.richter@icahn.mssm.edu
# 10/31/2017
# Description: run alignment commands, FastQC scripts and anything else fastq
#              related (e.g., maybe trimming in the future)
###############################################################################


# module load python/3.5.0
# module load py_packages/3.5
# module load fastqc/0.11.5
# module load hisat2/2.0.5
# python


import re
import subprocess
import glob
import os


class fq_pair(object):
    def __init__(self, pair_file_loc, home_dir, hisat2_idx):
        """Create an object for the FASTQ file."""
        # assign the filenames to the instance
        self.pair_file_loc = pair_file_loc
        # nuclear and ribo filenames have an annoying 001 at the end
        if re.search("nuclear|ribo", pair_file_loc):
            self.r1 = pair_file_loc + "_R1_001.fastq.gz"
            self.r2 = pair_file_loc + "_R2_001.fastq.gz"
        else:
            self.r1 = pair_file_loc + "_R1.fastq.gz"
            self.r2 = pair_file_loc + "_R2.fastq.gz"
        # obtain the Lane from the filename
        self.lane = re.search("L[0-9]{3}", pair_file_loc)
        # obtain barcode
        self.barcode = re.search("[ATCG]{6}", pair_file_loc)
        # create the sam file
        self.prefix = re.sub("_[ATCG]{6}", "", pair_file_loc)
        # go to home directory
        self.home_dir = home_dir # "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/"
    def RunHISAT2(self):
        """ Run HISAT2. Manual:
            https://ccb.jhu.edu/software/hisat2/manual.shtml
        """
        # hisat2 index
        self.hisat2_idx = hisat2_idx # "grcm38_snp_tran/genome_snp_tran"
        ## confirm sam file isn't already made, then run hisat2
        if not os.path.exists(self.prefix + ".sam"):
            hisat2_cmd = ("time hisat2 --time -x %s -1 %s -2 %s -S %s.sam " + \
                "--un-conc %s_noPEalign -p 24") % \
                (self.hisat2_idx, self.r1, self.r2, self.prefix, self.prefix)
            print(hisat2_cmd)
            subprocess.call(hisat2_cmd, shell = True)
        else: print(self.prefix + ".sam already made")
        """
        other HISAT2 options to consider implementing
        trim 5' or 3' ends (--trim3 <int> and --trim5 <int>)
        --phred33 or --phred64
        --known-splicesite-infile, --novel-splicesite-outfile,
             --novel-splicesite-infile
        reads that fail to align: --un <path> --un-conc <path>
        alignment metrics: --met-file <path> --met-stderr
        readgroup IDs with --rg-id <text>
        """
    def FastQC(self):
        """ FastQC command. Manual:
            https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
        """
        fastqc_cmd = "time fastqc %s"
        print("running fastqc for %s" % self.r1)
        subprocess.call(fastqc_cmd % self.r1, shell=True)
        print("running fastqc for %s" % self.r2)
        subprocess.call(fastqc_cmd % self.r2, shell=True)


fq_file_loc = "D1_D2_whole_cell/D1CTRL2_CGATGT_L001"
fq_i = fq_pair(pair_file_loc=fq_file_loc,
               home_dir="/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/",
               hisat2_idx="grcm38_snp_tran/genome_snp_tran")
os.chdir(fq_i.home_dir)
fq_i.RunHISAT2()


# /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/

home_dir = "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/"
hisat2_idx = "grcm38_snp_tran/genome_snp_tran"

# if I loop will it wait for subprocess to finish before proceeding? Yes, looks like it
# Want to loop since multi-threading the hisat2 command
# D1_D2_nuclear D1_D2_whole_cell D1_D2_ribo
fq_file_iter = glob.iglob(home_dir + "D1_D2_ribo/*.fastq.gz")
fq_file_list = [re.sub("_R[12].*fastq.gz", "", i) for i in fq_file_iter]

# keep uniques
len(fq_file_list)
fq_file_list = list(set(fq_file_list))
len(fq_file_list)

# fq_file_list = ["D1_D2_whole_cell/D1CTRL2_CGATGT_L001", "D1_D2_whole_cell/D1CTRL3_CGATGT_L001"]
# fq_file_list = ["D1_D2_nuclear/D1-1_TGACCA_L004", "D1_D2_nuclear/D1-2_GCCAAT_L005"]

for fq_i_file_loc in fq_file_list:
    fq_i = fq_pair(fq_i_file_loc, home_dir, hisat2_idx)
    os.chdir(fq_i.home_dir)
    fq_i.RunHISAT2()
    fq_i.FastQC()


"""
#### possibly deprecated 11/18/2017:
#
# def SamtoolsSort(self):
#     sort_cmd = "time samtools view -u %s | samtools sort -T %s.temp -o %s.sorted.bam -@ 12"
# def StringTieGTF(self):
#     stringtie_cmd = "time stringtie %s.sorted.bam -v -o %s_refguided.gtf -p 12" + \
#         "-G Mus_musculus.GRCm38.90.gtf"
#     ## after generating gtf for single file, merge with "--merge" option
#     # https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
# def RunFeatureCounts(self):
#      Run Feature Counts. Might just do this with R
#
#     pass
"""
