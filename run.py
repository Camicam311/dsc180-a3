import os
import sys
sys.path.append('./src/')
from utilities import *


def main():
    # Check argument length
    if len(sys.argv) != 2: 
        print('Usage:\n\
               python ./src/run.py {run-project | run-ftp | test-project} ')
        return
    cfg_opt = sys.argv[1]
    cfg_file = None
    # Make sure config file exists
    if cfg_opt == 'run-project':
        cfg_file = './config/data.json'
    elif cfg_opt == 'test-project':
        cfg_file = './config/data-test.json'
    elif cfg_opt == 'run-ftp':
        cfg_file = './config/data-ftp.json'
    else:
        print('Usage:\n\
               python ./src/run.py {run-project | run-ftp | test-project} ')
        return
    
    # Load config
    cfg = json.load(open(cfg_file))
    
    # Download data
    datapath = cfg['datapath']
    download(datapath, **cfg['data'])
    
    # Process data
    process(datapath, cfg['data']['ref_file'], cfg['options'])
        
main()