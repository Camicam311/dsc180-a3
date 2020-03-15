import numpy as np
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


def main():
    for arg in sys.argv:
        if arg == 'data':
            data()
        elif arg == 'data-test':
            data_test()
        elif arg == 'process':
            process()
    
def data():
    
    # This code could be changed to download from 1000 Genomes FTP server
    cfg = json.load(open('config/data.json'))
    
    datapath = cfg['datapath']
    if not os.path.exists(datapath):
        os.mkdir(datapath)
    
    procs = []
    for vcf_file in cfg['vcf_files']:
        cmd = shlex.split('cp -r ' + vcf_file + ' ' + datapath)
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        procs.append(proc)
        
    
    print("Copying Files...")
    copytime = 0
    while procs[-1].poll() == None:
        time.sleep(1)
        copytime += 1
        loadingstr = 'Copying'
        for i in range(copytime % 3):
            loadingstr += '.'
        print(loadingstr, end = '\r')
    print('                ', end = '\r')
    print('Done', end = '\r')
    
        
def data_test():
    
    cfg = json.load(open('config/data-test.json'))
    
    datapath = cfg['datapath']
    if not os.path.exists(datapath):
        os.mkdir(datapath)
    
    procs = []
    for vcf_file in cfg['vcf_files']:
        cmd = shlex.split('cp -r ' + vcf_file + ' ' + datapath + re.findall("/.*|$", vcf_file)[-1][1:])
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        procs.append(proc)
    
    
    print("Copying Files...")
    copytime = 0
    while procs[-1].poll() == None:
        time.sleep(1)
        copytime += 1
        loadingstr = 'Copying'
        for i in range(copytime % 3):
            loadingstr += '.'
        print(loadingstr, end = '\r')
    print('                ', end = '\r')
    print('Done')
        
def process():
    
    cfg = json.load(open('config/process.json'))
    datapath = cfg['datapath']
    outpath = cfg['outpath']
    path_to_plink = cfg['path_to_plink']
    os.makedirs(outpath + 'plots/', exist_ok=True)
    for vcf_file in os.listdir(datapath):
        outname = re.findall('chr[0-9]{1,2}|$', vcf_file)[0]
        outname = outname if outname != '' else 'unknown'
        os.makedirs(outpath + outname, exist_ok=True)
        # shlex.split parses a terminal command into a list for subprocess
        cmd = shlex.split( path_to_plink + '\
            --vcf ' + datapath + vcf_file + ' \
            --make-bed \
            --pca approx biallelic-var-wts\
            --snps-only \
            --maf 0.05 \
            --geno 0.1 \
            --mind 0.05 \
            --recode \
            --out ' + outpath + outname + '/' + outname)
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        
        runtime = 0
        while proc.poll() == None:
            time.sleep(1)
            runtime += 1
            print('Finding PCA for ' + outname + '. Runtime: %d secs' % (runtime), end = '\r')
        print('                ', end = '\r')
        print('Done')
        vec = pd.read_csv(outpath + outname + '/' + outname + '.eigenvec', sep='\t')
        val = pd.read_csv(outpath + outname + '/' + outname + '.eigenval')
        pca_plot = sns.scatterplot(vec['PC1'], vec['PC2'])
        pca_plot.figure.savefig(outpath + 'plots/' + outname + '.png')
        print('PCA Plot saved to ' + outpath + 'plots/')
        
main()