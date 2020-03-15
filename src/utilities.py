import os
import shutil
import shlex

def download_fastqs(fastq_files, ref_file, datapath):
    """
    Copies/Downloads FASTQ files from defined paths to the local data directory

    Keyword arguments:
    fastq_files -- array of strings representing paths to fastq files
    ref_file -- string representing the path to the reference genome
    datapath -- string representing the local path to copy/download data into
    """
    for fastq_file in fastq_files:
        cmd = shlex.split('cp -r ' + fastq_file + ' ' + datapath)
        sp.run(cmd)

def download_bams(bam_files, datapath):
    """
    Copies/Downloads BAM files from defined paths to the local data directory

    Keyword arguments:
    bam_files -- array of strings representing paths to bam files
    datapath -- string representing the local path to copy/download data into
    """
    for bam_file in bam_files:
        cmd = shlex.split('cp -r ' + bam_file + ' ' + datapath)
        sp.run(cmd)


def download_vcfs(vcf_files, datapath):
    """
    Copies/Downloads VCF files from defined paths to the local data directory,
    and unpacks zipped files and edits the first line to make sure the VCF
    version is compatible with vcftools

    Keyword arguments:
    vcf_files -- array of strings representing paths to vcf files
    datapath -- string representing the local path to copy/download data into
    """
    for vcf_file in vcf_files:
        cmd = shlex.split('cp -r ' + vcf_file + ' ' + datapath)
        sp.run(cmd)
        cmd = shlex.split('cp -r ' + vcf_file + '.tbi ' + datapath)
        sp.run(cmd)


def process_fastq(datapath):
    for fastq_file in os.listdir(datapath + 'fastq_files/'):
        cmd = shlex.split('bwa mem ' + path_to_ref + ' ' + fastq_file + ' | samtools sort -o output.bam -')
        proc = sp.run(cmd)

def process_bam(datapath):
    for bam_file in os.listdir(datapath + 'bam_files/'):
        pass

def process_vcfs(datapath, maf=0.05, geno=0.1, mind=0.05):
    """
    Processes VCF files and outputs a PCA plot

    Keyword arguments:
    datapath -- String representing path to directory containing VCF files
    maf -- Minimum frequency of an allele for it to be included in PCA
    geno -- Maximum missingness of variant before variant is excluded
    mind -- Maximum missingness of sample before sample is excluded
    """
    fix_vcf_version(datapath)
    filter_and_combine_vcfs(datapath)

    for vcf_file in os.listdir(datapath + 'vcf_files/'):
        outname = re.findall('chr[0-9]{1,2}|$', vcf_file)[0]
        outname = outname if outname != '' else 'unknown'
        os.makedirs(outpath + outname, exist_ok=True)
        # shlex.split parses a terminal command into a list for subprocess
        cmd = shlex.split( path_to_plink + '\
            --vcf ' + datapath + vcf_file + ' \
            --make-bed \
            --pca approx biallelic-var-wts\
            --snps-only \
            --maf ' + maf + ' \
            --geno ' + geno + ' \
            --mind ' + mind + ' \
            --recode \
            --out ' + outpath + outname + '/' + outname)
        sp.run(cmd)
        vec = pd.read_csv(outpath + outname + '/' + outname + '.eigenvec', sep='\t')
        val = pd.read_csv(outpath + outname + '/' + outname + '.eigenval')
        pca_plot = sns.scatterplot(vec['PC1'], vec['PC2'])
        pca_plot.figure.savefig(outpath + 'plots/' + outname + '.png')
        print('PCA Plot saved to ' + outpath + 'plots/')

def fix_vcf_version(datapath):
    """
    Changes all the VCF files in the directory to be compatible with vcftools

    Keyword arguments:
    datapath -- string representing the absolute path to the directory containing vcf files
    """
    # Prepare vcf files for analysis by making sure they are compatible with vcftools
    for vcf_file in os.listdir(datapath):
        # if zipped, unzip vcf file and assert that version is < 4.2
        if vcf_file[-7:] == '.vcf.gz':
            with gzip.open(vcf_file, 'rt') as f:
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
                                     datapath + fname)
                    sp.run(cmd)
        # assert that version is < 4.2
        elif vcf_file[-4:] == '.vcf':
            with open(vcf_file, 'rt') as f:
                head = f.readline()
                if 'VCFv4.3' in head:
                    fname = re.findall("([^\/]+$)", vcf_file)[-1][:-3]
                    # replace version number
                    cmd = shlex.split("sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' " + \
                                     datapath + fname)
                    sp.run(cmd)

def filter_and_combine_vcfs(datapath, snp_treshold=0.05):
    """
    Filters out variants below a specific treshold, combines vcfs into a single vcfs

    Keyword arguments:
    datapath -- string representing the absolute path to the directory containing vcf files
    snp_treshold -- how common SNP must be to be included
    """
    



