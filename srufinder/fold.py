import os
import subprocess
import logging
import sys
import shutil

class Fold(object):
    
    def __init__(self, obj):
        self.master = obj

    def run(self):
        '''
        Running RNAfold to predict secondary structure of repeats
        '''

        logging.info('Predicting secondary structure of repeat sequences')

        # RNAfold
        if shutil.which('RNAfold') is not None:
            os.mkdir(self.master.out+'RNAfold')

            subprocess.run(['RNAfold', 
                            '-i', self.master.out+'repeats.fna',
                            '-o', os.path.join(self.master.out, 'RNAfold', 'repeat')])
        else:
            logging.error('RNAfold not found. Skipping secondary structure prediction')
