import os
import subprocess
import logging
import sys

class Blast(object):
    
    def __init__(self, obj):
        self.master = obj

    def make_db(self):
        '''
        Make a BLAST database from the masked input sequence
        '''

        logging.debug('Making BLAST database')

        subprocess.run(['makeblastdb', 
                        '-dbtype', 'nucl', 
                        '-in', self.master.out+'masked.fna',
                        '-out', self.master.out+'masked'], 
                        stdout=subprocess.DEVNULL)
    
    def run(self):
        '''
        BLASTing repeat database against the masked input sequence
        '''

        # Make the database
        self.make_db()

        logging.info('BLASTing repeats')

        # BLASTn
        subprocess.run(['blastn', 
                        '-task', 'blastn-short', 
                        '-word_size', str(self.master.word_size), 
                        '-query', self.master.repeatdb,
                        '-db', self.master.out+'masked',
                        '-outfmt', '6',
                        '-out', self.master.out+'blast.tab',
                        '-num_threads', str(self.master.threads)])
