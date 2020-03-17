import pandas as pd
import os
import shutil
import shlex
import sys
import subprocess as sp
import time
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import re
import json
import time
import gzip
from ftplib import FTP

def download(datapath, fastq_files, bam_files, vcf_files, ref_file, ftp_path=None):
    """
    Downloads all the files

    Keyword arguments:
    datapath -- string representing the absolute path to the directory to download data into
    fastq_files -- list of fastq files to download
    bam_files -- list of bam files to download
    vcf_files -- list of vcf files to download
    ref_file -- path to reference genome file to download
    ftp_path -- path to ftp server, if this is defined then download will use ftp
    """
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
    """
    Processes all the files

    Keyword arguments:
    datapath -- string representing the absolute path to the directory to download data into
    ref_file -- path to reference genome file, which has been downloaded locally
    options -- options for the PCA analysis
    """
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
    t1 = time.time()
    print("Processing VCF files...")
    process_vcf(datapath, **options)
    t2 = time.time()
    print("VCF files done. Time: {0:.0f} s".format(t2 - t1))
    print("Creating plot...")
    create_plot()
    print("Plot done. Saved to final_plot.png")
    
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
    """
    Processes FASTQ files and converts them into BAM files

    Keyword arguments:
    datapath -- path pointing to where the FASTQ files have been downloaded
    ref_file -- path to reference genome file
    """
    ref_file = re.findall("([^/]+$)", ref_file)[-1]
    path_to_ref = datapath + 'ref_genome/' + ref_file
    for fastq_file in os.listdir(datapath + 'fastq_files/'):
        fname = re.findall(".+?((?=\.fastq)|(?=\.fq))", fastq_file)[-1]
        cmd = shlex.split('bwa mem ' + path_to_ref + ' ' + datapath + 'fastq_files/'\
                          + fastq_file + ' | samtools sort -o ' + fname + '.bam -')
        sp.run(cmd)

def process_bam(datapath, ref_file):
    """
    Processes BAM files and converts them into VCF files

    Keyword arguments:
    datapath -- path pointing to where the VCF files have been downloaded
    ref_file -- path to reference genome file
    """
    ref_file = re.findall("([^/]+$)", ref_file)[-1]
    path_to_ref = datapath + 'ref_genome/' + ref_file
    cmd = 'samtools mpileup ' + path_to_ref + ' '
    for bam_file in os.listdir(datapath + 'bam_files/'):
         cmd += bam_file + ' '
    cmd += datapath + 'vcf_files/bam_output.vcf'
    cmd = shlex.split(cmd)
    sp.run(cmd)

def process_vcf(datapath, maf=0.05, geno=0.1, mind=0.05, sample=0.01):
    """
    Processes VCF files into PCA eigenvector and eigenvalue

    Keyword arguments:
    datapath -- String representing path to directory containing VCF files
    maf -- Minimum frequency of an allele for it to be included in PCA
    geno -- Maximum missingness of variant before variant is excluded
    mind -- Maximum missingness of sample before sample is excluded
    """
    if not os.path.exists(datapath + 'output/'):
        os.mkdir(datapath + 'output/')
    if not os.path.exists(datapath + 'final_vcfs/'):
        os.mkdir(datapath + 'final_vcfs/')
    # Filter and merge VCFs
    filter_and_combine_vcfs()
    # Take a small sample from the merged VCF
    cmd = 'gatk SelectVariants \
        -V ' + datapath + '/merged_vcf.vcf \
        -fraction ' + str(sample) + ' \
        -O ' + datapath + '/final.vcf'
    cmd = shlex.split(cmd)
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    t1 = time.time()
    while proc.poll() == None:
        print('Sampling from VCF... Time: {0:.0f} s'.format(time.time() - t1), end='\r')
        time.sleep(0.1)
    print()
    # run PCA
    cmd = 'plink2 \
        --vcf ' + datapath + '/final.vcf \
        --make-bed \
        --pca approx biallelc\
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
        print('Finding PCA... Time: {0:.0f} s'.format(time.time() - t1), end='\r')
        time.sleep(0.1)
    print()

def filter_and_combine_vcfs():
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
        print('Filtering VCFs... Time: {0:.0f} s'.format(time.time() - t1), end='\r')
        time.sleep(1)
    print()
    
    
def create_plot():
    """
    Creates the final plots of PCA colored by population
    """
    if not os.path.exists('plots/'):
        os.mkdir('plots/')
    # Read PCA values and populations
    vec = pd.read_csv('./data/output/final.eigenvec', sep=' ', header=None)
    vec = vec.drop(0, axis=1)
    val = pd.read_csv('./data/output/final.eigenval', header=None)
    sampop = pd.read_excel('pop_data/sample_populations.xlsx')
    popnames = pd.read_csv('pop_data/population_names.tsv', sep='\t')

    # Clean column names
    cols = ['Sample']
    cols.extend(['PCA' + str(i) for i in range(1, 21)])
    vec.columns = cols

    # Join PCA vector with population values
    vec = pd.merge(vec, sampop[['Sample', 'Population', 'Population Description']], left_on='Sample', right_on='Sample')
    vec = pd.merge(vec, popnames[['Population Code', 'Super Population']], left_on='Population', right_on='Population Code')

    # Define colors for populations
    colors = ['#FF6666', '#FF8000', '#B2FF66', '#6666FF', '#FF66B2']

    # Create plot
    sns.set_style("whitegrid", {'axes.grid' : False})

    fig = plt.figure(figsize=(6,6))
    ax = Axes3D(fig)
    for i, pop in enumerate(vec['Super Population'].unique()):
        subset = vec[vec['Super Population'] == pop]
        ax.scatter(subset['PCA1'], subset['PCA2'], subset['PCA3'], c=colors[i], marker='o', label=pop)
        ax.set_xlabel('PCA1')
        ax.set_ylabel('PCA2')
        ax.set_zlabel('PCA3')

    ax.legend()
    # Save figure
    fig.savefig('plots/final_plot.png')
    
    
