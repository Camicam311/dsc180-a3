import os
import sys
sys.path.append('./src/')
from utilities import *


def main():
    # Check argument length
    if len(sys.argv) != 2: 
        print('Usage:\n\
               python ./src/run.py {data | data-test | data-ftp} ')
        return
    cfg_opt = sys.argv[1]
    cfg_file = None
    # Make sure config file exists
    if cfg_opt == 'data':
        cfg_file = './config/data.json'
    elif cfg_opt == 'data-test':
        cfg_file = './config/data-test.json'
    elif cfg_opt == 'data-test':
        cfg_file = './config/data-ftp.json'
    else:
        print('Usage:\n\
               python ./src/run.py {data | data-test | data-ftp} ')
        return
    
    # Load config
    cfg = json.load(open(cfg_file))
    
    # Download data
    datapath = cfg['datapath']
    download(datapath, **cfg['data'])
    
    # Process data
    if sys.argv[2] == 'process':
        process(datapath, cfg['data']['ref_file'], cfg['options'])
        
main()