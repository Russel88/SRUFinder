import os
import subprocess
import logging
import sys

class Blast(object):
    
    def __init__(self, obj):
        self.master = obj

    def make_db(self, indb, outdb):
        '''
        Make a BLAST database from the masked input sequence
        '''

        logging.debug('Making BLAST database')

        subprocess.run(['makeblastdb', 
                        '-dbtype', 'nucl', 
                        '-in', indb,
                        '-out', outdb], 
                        stdout=subprocess.DEVNULL)
    
    def run(self):
        '''
        BLASTing repeat database against the masked input sequence
        '''

        # Make the database
        self.make_db(self.master.out+'masked.fna', self.master.out+'masked')

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
    
    def run_spacer(self):
        '''
        BLASTing repeat database against the masked input sequence
        '''

        # Make the database
        self.make_db(self.master.fasta, self.master.out+'genome')

        logging.info('BLASTing spacers against self')

        # BLASTn
        subprocess.run(['blastn', 
                        '-task', 'blastn-short', 
                        '-word_size', str(self.master.word_size), 
                        '-query', self.master.out+'spacers.fa',
                        '-db', self.master.out+'genome',
                        '-outfmt', '6',
                        '-perc_identity', str(self.master.spacer_identity),
                        '-qcov_hsp_perc', str(self.master.spacer_coverage),
                        '-out', self.master.out+'blast_spacers.tab',
                        '-num_threads', str(self.master.threads)])
