import os
import subprocess
import logging
import sys

class Fold(object):
    
    def __init__(self, obj):
        self.master = obj

    def run(self):
        '''
        Running RNAfold to predict secondary structure of repeats
        '''

        logging.info('Predicting secondary structure of repeat sequences')

        # RNAfold
        os.mkdir(self.master.out+'RNAfold')

        subprocess.run(['RNAfold', 
                        '-i', self.master.out+'repeats.fna',
                        '-o', os.path.join(self.master.out, 'RNAfold', 'repeat')])
    
