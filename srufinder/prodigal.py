import os
import subprocess
import logging
import sys
import re

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from shutil import copyfile

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

        try:
            genes = pd.read_csv(self.master.out+'prodigal.gff', sep='\t|;[a-z,A-Z,_]*=', comment="#", engine='python', header=None)
       
            # Remove low confidence
            self.genes = genes[genes.iloc[:,14] >= self.master.orf]
        
            self.noorf = False
        
        except:
            logging.warning('No ORFs found. Skipping masking')
            self.noorf = True

    def mask(self):
        '''
        Masking input by replacing all ORF or non-ORF sequences with N's
        '''

        if self.noorf:
            copyfile(self.master.fasta, self.master.out+'masked.fna')
        else:

            logging.info('Masking input sequence')
            
            with open(self.master.out+'masked.fna', 'w') as out_file:
                falist = SeqIO.parse(open(self.master.fasta, 'r'), 'fasta')
                # For each sequence
                for fas in falist:
                    
                    name = str(fas.id)
                    seq = str(fas.seq)
                    Xsub = self.genes[[x == name for x in self.genes.iloc[:,0]]]
                    
                    # If non-ORFs should be masked
                    if self.master.in_orf:
                        # If no ORFs, all is intergenic
                        if Xsub.empty:
                            seq = 'N'*len(seq)
                        else:
                            # Ensure the order is correct
                            Xsub = Xsub.sort_values(by = 3)
                            # Get intergenic regions
                            starts = [1] + list(Xsub.iloc[:,4])
                            ends = list(Xsub.iloc[:,3]) + [len(fas.seq)]
                            # Combine
                            overlap = [x[1]-x[0] for x in zip(starts, ends)]
                            intergenic = zip(starts, ends, overlap)
                            # Remove empty intergenic regions
                            intergenic = [x for x in intergenic if x[2] > 1]
                            # For each intergenic region
                            for reg in intergenic:
                                Xfrom = int(reg[0])
                                Xto = int(reg[1])
                                # New sequence
                                seq1 = seq[:Xfrom]
                                seqX = 'N'*(Xto-Xfrom-1)
                                seq2 = seq[Xto-1:]
                                seq = seq1+seqX+seq2
                    
                    # If ORFs should be masked    
                    else:
                        # Only fastas found in Xtable
                        if not Xsub.empty:
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

