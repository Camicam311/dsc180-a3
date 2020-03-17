import pandas as pd
import os
import shutil
import shlex
import sys
import subprocess as sp
import time
import matplotlib.pyplot as plt
import seaborn as sns
import re
import json
import time
import gzip
from ftplib import FTP

def download(datapath, fastq_files, bam_files, vcf_files, ref_file, ftp_path=None):
    # Create data directories
    if not os.path.exists(datapath):
        os.mkdir(datapath)
    if not os.path.exists(datapath + 'fastq_files/'):
        os.mkdir(datapath + 'fastq_files/')
    if not os.path.exists(datapath + 'bam_files/'):
        os.mkdir(datapath + 'bam_files/')
    if not os.path.exists(datapath + 'vcf_files/'):
        os.mkdir(datapath + 'vcf_files/')
    if not os.path.exists(datapath + 'ref_genome/'):
        os.mkdir(datapath + 'ref_genome/')
    
    f = None
    if ftp_path:
        f = FTP(ftp_path)
        f.login()

    # Copy or download data
    print("--------------------Downloading--------------------")
    print("Downloading FASTQ files...")
    t1 = time.time()
    download_fastqs(fastq_files, ref_file, datapath, f)
    t2 = time.time()
    print("FASTQ files done. Time: {0:.0f} s".format(t2 - t1))
    print("Downloading BAM files...")
    download_bams(bam_files, datapath + 'bam_files/', f)
    t1 = time.time()
    print("BAM files done. Time: {0:.0f} s".format(t1 - t2))
    print("Downloading VCF files...")
    download_vcfs(vcf_files, datapath + 'vcf_files/', f)
    t2 = time.time()
    print("VCF files done. Time: {0:.0f} s".format(t2 - t1))
        
def process(datapath, ref_file, options):
    
    print("--------------------Processing---------------------")
    # Still having issues with FASTQ and BAM
#     print("Processing FASTQ files...")
#     t1 = time.time()
#     process_fastq(datapath, ref_file)
#     t2 = time.time()
#     print("FASTQ files done. Time: {0:.0f} s".format(t2 - t1))
#     print("Processing BAM files...")
#     process_bam(datapath, ref_file)
#     t1 = time.time()
#     print("BAM files done. Time: {0:.0f} s".format(t1 - t2))
    print("Processing VCF files...")
    process_vcf(datapath, **options)
    t2 = time.time()
    print("VCF files done. Time: {0:.0f} s".format(t2 - t1))
    
def download_fastqs(fastq_files, ref_file, datapath, ftp=None):
    """
    Copies/Downloads FASTQ files from defined paths to the local data directory

    Keyword arguments:
    fastq_files -- array of strings representing paths to fastq files
    ref_file -- string representing the path to the reference genome
    datapath -- string representing the local path to copy/download data into
    """
    if ftp:
        # Download from ftp server
        for fastq_file in fastq_files:
            # Split file up into directory and file name
            lastsep = fastq_file.rindex("/") + 1
            fqpath = fastq_file[:lastsep]
            fqname = fastq_file[lastsep:]
            with open(datapath + 'fastq_files/' + fqname, 'wb') as fp:
                # Change directory then download file
                ftp.cwd(fqpath)
                ftp.retrbinary('RETR ' + fqname, fp.write)
                fp.close()
        # Split file up into directory and file name for reference genome
        lastsep = ref_file.rindex("/") + 1
        refpath = ref_file[:lastsep]
        refname = ref_file[lastsep:]
        with open(datapath + 'ref_genome/' + refname, 'wb') as fp:
            # Change directory then download file
            ftp.cwd(refpath)
            ftp.retrbinary('RETR ' + refname, fp.write)
            fp.close()
    else:
        # Copy files into local directory
        for fastq_file in fastq_files:
            cmd = shlex.split('cp -r ' + fastq_file + ' ' + datapath + 'fastq_files/')
            sp.run(cmd)
        for ext in ['', '.gz', '.amb', '.ann', '.bwt', '.pac', '.sa', '.fai']:
            cmd = shlex.split('cp -r ' + ref_file + ext + ' ' + datapath + 'ref_genome/')
            sp.run(cmd)

def download_bams(bam_files, datapath, ftp=None):
    """
    Copies/Downloads BAM files from defined paths to the local data directory

    Keyword arguments:
    bam_files -- array of strings representing paths to bam files
    datapath -- string representing the local path to copy/download data into
    """
    if ftp:
        for bam_file in bam_files:
            lastsep = bam_file.rindex("/") + 1
            bampath = bam_file[:lastsep]
            bamname = bam_file[lastsep:]
            with open(datapath + bamname, 'wb') as fp:
                ftp.cwd(bampath)
                ftp.retrbinary('RETR ' + bamname, fp.write)
                fp.close()
    else:
        for bam_file in bam_files:
            cmd = shlex.split('cp -r ' + bam_file + ' ' + datapath)
            sp.run(cmd)

def download_vcfs(vcf_files, datapath, ftp=None):
    """
    Copies/Downloads VCF files from defined paths to the local data directory,
    and unpacks zipped files and edits the first line to make sure the VCF
    version is compatible with vcftools

    Keyword arguments:
    vcf_files -- array of strings representing paths to vcf files
    datapath -- string representing the local path to copy/download data into
    """
    if ftp:
        for vcf_file in vcf_files:
            lastsep = vcf_file.rindex("/") + 1
            vcfpath = vcf_file[:lastsep]
            vcfname = vcf_file[lastsep:]
            with open(datapath + vcfname, 'wb') as fp:
                ftp.cwd(vcfpath)
                ftp.retrbinary('RETR ' + vcfname, fp.write)
                fp.close()
            with open(datapath + vcfname + '.tbi', 'wb') as fp:
                ftp.cwd(vcfpath)
                ftp.retrbinary('RETR ' + vcfname + '.tbi', fp.write)
                fp.close()
    else:
        for vcf_file in vcf_files:
            cmd = shlex.split('cp -r ' + vcf_file + ' ' + datapath)
            sp.run(cmd)
            cmd = shlex.split('cp -r ' + vcf_file + '.tbi ' + datapath)
            sp.run(cmd)

def process_fastq(datapath, ref_file):
    ref_file = re.findall("([^/]+$)", ref_file)[-1]
    path_to_ref = datapath + 'ref_genome/' + ref_file
    for fastq_file in os.listdir(datapath + 'fastq_files/'):
        fname = re.findall(".+?((?=\.fastq)|(?=\.fq))", fastq_file)[-1]
        cmd = shlex.split('bwa mem ' + path_to_ref + ' ' + datapath + 'fastq_files/'\
                          + fastq_file + ' | samtools sort -o ' + fname + '.bam -')
        sp.run(cmd)

def process_bam(datapath, ref_file):
    ref_file = re.findall("([^/]+$)", ref_file)[-1]
    path_to_ref = datapath + 'ref_genome/' + ref_file
    cmd = 'samtools mpileup ' + path_to_ref + ' '
    for bam_file in os.listdir(datapath + 'bam_files/'):
         cmd += bam_file + ' '
    cmd += datapath + 'vcf_files/bam_output.vcf'
    cmd = shlex.split(cmd)
    sp.run(cmd)

def process_vcf(datapath, maf=0.05, geno=0.1, mind=0.05):
    """
    Processes VCF files and outputs a PCA plot

    Keyword arguments:
    datapath -- String representing path to directory containing VCF files
    maf -- Minimum frequency of an allele for it to be included in PCA
    geno -- Maximum missingness of variant before variant is excluded
    mind -- Maximum missingness of sample before sample is excluded
    """
    filter_and_combine_vcfs(datapath + 'vcf_files/')
    if not os.path.exists(datapath + 'output/'):
        os.mkdir(datapath + 'output/')
    if not os.path.exists(datapath + 'plots/'):
        os.mkdir(datapath + 'plots/')
    if not os.path.exists(datapath + 'final_vcfs/'):
        os.mkdir(datapath + 'final_vcfs/')
    filter_and_combine_vcfs()
    # merge VCFs into one VCF
    cmd = shlex.split('bcftools merge ' + ' '.join(files) + ' | bgzip -c > ' + datapath + 'merged_vcf.vcf ')
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    while proc.poll() == None:
        time.sleep(0.1)
    cmd = shlex.split('bcftools index -t ' + datapath + 'merged_vcf.vcf')
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    while proc.poll() == None:
        time.sleep(0.1)
    # run PCA
    cmd = 'plink2 \
        --vcf ' + datapath + 'final_vcfs/merged_vcf.vcf \
        --make-bed \
        --pca \
        --snps-only \
        --maf ' + str(maf) + ' \
        --geno ' + str(geno) + ' \
        --mind ' + str(mind) + ' \
        --recode \
        --out ' + datapath + 'output/final'
    cmd = shlex.split(cmd)
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    t1 = time.time()
    while proc.poll() == None:
        print('Finding PCA... Time: {0:1f} s'.format(time.time() - t1), end='\r')
        time.sleep(0.1)
    print()
    vec = pd.read_csv(datapath + 'output/final.eigenvec', sep='\t')
    val = pd.read_csv(datapath + 'output/final.eigenval')
    pca_plot = sns.scatterplot(vec['PC1'], vec['PC2'])
    pca_plot.figure.savefig( + 'plots/final.png')
    print('PCA Plot saved to ' + datapath + 'plots/final.png')

def fix_vcf_version(datapath):
    """
    Changes all the VCF files in the directory to be compatible with vcftools

    Keyword arguments:
    datapath -- string representing the absolute path to the directory containing vcf files
    """
    to_fix = []
    t1 = time.time()
    # Prepare vcf files for analysis by making sure they are compatible with vcftools
    for vcf_file in os.listdir(datapath):
        # if zipped, unzip vcf file and assert that version is < 4.2
        if vcf_file[-7:] == '.vcf.gz':
            with gzip.open(datapath + vcf_file, 'rt') as f:
                head = f.readline()
                # if version not compatible with vcftools, have to fix it
                if 'VCFv4.3' in head:
                    fname = re.findall("([^\/]+$)", vcf_file)[-1]
                    # unpack file to replace version number
                    cmd = shlex.split('gzip -d ' + datapath + fname)
                    sp.run(cmd)
                    fname = fname[:-3]
                    # replace version number
                    cmd = shlex.split("sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' " + \
                                     datapath + fname + ' > ' + datapath + fname[:-3] + '.tmp.vcf')
                    sp.run(cmd)
                    cmd = shlex.split("rm " + datapath + fname)
                    sp.run(cmd)
                    cmd = shlex.split("mv " + datapath + fname[:-3] + '.tmp.vcf ' + datapath + fname)
                    sp.run(cmd)
        # assert that version is < 4.2
        elif vcf_file[-4:] == '.vcf':
            with open(datapath + vcf_file, 'rt') as f:
                head = f.readline()
                if 'VCFv4.3' in head:
                    fname = re.findall("([^\/]+$)", vcf_file)[-1][:-3]
                    # replace version number
                    cmd = shlex.split("sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' " + \
                                     datapath + fname + ' > ' + datapath + fname[:-3] + '.new.vcf')
                    sp.run(cmd)
                    cmd = shlex.split("rm " + datapath + fname)
                    sp.run(cmd)
                    cmd = shlex.split("mv " + datapath + fname[:-3] + '.new.vcf ' + datapath + fname)
                    sp.run(cmd)

def filter_and_combine_vcfs(datapath, snp_treshold=0.05):
    """
    Filters out variants below a specific treshold, combines vcfs into a single vcfs

    Keyword arguments:
    datapath -- string representing the absolute path to the directory containing vcf files
    snp_treshold -- how common SNP must be to be included
    """
    cmd = shlex.split('./src/filter_vcfs.sh')
    proc = sp.Popen(cmd)
    t1 = time.time()
    while proc.poll() == None:
        print('Filtering VCFs... Time: {0:1f} s'.format(time.time() - t1), end='\r')
        time.sleep(0.1)
    print()
    