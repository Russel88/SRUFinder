import os
import subprocess
import logging
import sys

class Bprom(object):
    
    def __init__(self, obj):
        self.master = obj

    def run(self):
        '''
        Running Bprom to predcit promoters
        '''

        logging.info('Predicting promoters in flanking sequences')

        # Bprom
        subprocess.run(['bprom', 
                        self.master.out+'flanking.fna',
                        self.master.out+'promoters'])
    
