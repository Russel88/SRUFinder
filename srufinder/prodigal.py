import os
import subprocess
import logging
import sys
import re

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq

class Prodigal(object):
    
    def __init__(self, obj):
        self.master = obj

    def run(self):
        '''
        Running prodigal to predict ORFs in the input sequence.
        Then check if the run was succesful, load the gff file, 
        and create a sequence masked with the ORFs
        '''

        logging.info('Predicting ORFs with prodigal')

        # Run prodigal
        with open(self.master.out+'prodigal.gff', 'w') as prodigal_out:
            subprocess.run(['prodigal', 
                            '-i', self.master.fasta, 
                            '-p', self.master.prod,
                            '-f', 'gff'], 
                            stdout=prodigal_out, 
                            stderr=subprocess.DEVNULL)

        # Check if succesful
        self.check()
        
        # Load genes and filter
        self.get_genes()

        # Mask fasta
        self.mask()

    def check(self):
        '''
        Check if the prodigal output has a size larger than 0 
        else terminate
        '''

        logging.debug('Checking prodigal output')

        # Check prodigal output
        if os.stat(self.master.out+'prodigal.gff').st_size == 0:
            logging.critical('Prodigal failed!')
            sys.exit()

    def get_genes(self):
        '''
        Load the prodigal gff file and save a dataframe,
        where low confident ORFs have been removed
        '''

        logging.debug('Loading prodigal GFF')

        genes = pd.read_csv(self.master.out+'prodigal.gff', sep='\t|;[a-z,A-Z,_]*=', comment="#", engine='python')
       
        # Remove low confidence
        self.genes = genes[genes.iloc[:,14] >= self.master.orf]

    def mask(self):
        '''
        Masking input by replacing all ORF sequences with N's
        '''

        logging.info('Masking input sequence')
        
        with open(self.master.out+'masked.fna', 'w') as out_file:
            falist = SeqIO.parse(open(self.master.fasta, 'r'), 'fasta')
            # For each sequence
            for fas in falist:
                name = str(fas.id)
                Xsub = self.genes[[x == name for x in self.genes.iloc[:,0]]]
                # Only fastas found in Xtable
                if not Xsub.empty:
                    seq = str(fas.seq)
                    # Each row of Xtable
                    for row in Xsub.itertuples():
                        # From where to where
                        Xfrom = row[4]
                        Xto = row[5]
                        # New sequence
                        seq1 = seq[:int(Xfrom) - 1]
                        seqX = 'N'*(Xto-Xfrom+1)
                        seq2 = seq[int(Xto):]
                        seq = seq1+seqX+seq2
                    fas.seq = Seq(seq)
                # Write sequence
                SeqIO.write(fas, out_file, "fasta")

