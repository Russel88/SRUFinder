import os
import subprocess
import logging
import sys
import shutil

class Bprom(object):
    
    def __init__(self, obj):
        self.master = obj

    def run(self):
        '''
        Running Bprom to predcit promoters
        '''

        logging.info('Predicting promoters in flanking sequences')

        # Bprom
        if shutil.which('bprom') is not None:
            for flfile in os.listdir(self.master.out+'flanking'):
                subprocess.run(['bprom', 
                                os.path.join(self.master.out, 'flanking', flfile),
                                os.path.join(self.master.out, 'flanking', flfile)+'_bprom'])
        else:
            logging.error('Bprom not found. Skipping promoter prediction')
