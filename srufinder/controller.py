import os
import logging
import sys

import pandas as pd

from Bio import SeqIO

class Controller(object):

    def __init__(self, args):
       
        self.fasta = args.input
        self.out = args.output
        self.prod = args.prodigal
        self.db = args.db
        self.orf = args.orf

        # Logger
        logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        logging.info('Running SRUFinder version 0.0.1')

        # Force consistency
        self.out = os.path.join(self.out, '')

        # Check databases
        self.check_db()
        
        # Check input and output
        self.check_input()
        self.check_out()

        # Write arguments
        da = vars(args)
        f = open(self.out+'arguments.tab', 'w')
        for k, v in da.items():
            f.write('{}:\t{}\n'.format(k, v))
        f.close()

    def check_out(self):

        if not self.redo:
            try:
                os.mkdir(self.out)
            except FileExistsError:
                logging.error('Directory '+self.out+' already exists')
                sys.exit()

    def check_input(self):

        if not self.check_inp:
            if os.path.isfile(self.fasta):
                if not self.is_fasta():
                    logging.error('Input file is not in fasta format')
                    sys.exit()
            else:
                logging.error('Could not find input file')
                sys.exit()

    def is_fasta(self):
        
        try:
            with open(self.fasta, 'r') as handle:
                fa = SeqIO.parse(handle, 'fasta')
                [float(x.id) for x in fa]
                logging.error('Numeric fasta headers not supported')
                return False
        except:
            with open(self.fasta, 'r') as handle:
                fa = SeqIO.parse(handle, 'fasta')
                return any(fa)

    def check_db(self):
        
        if self.db == '':
            try:
                self.db = os.environ['SRUFINDER_DB']
            except:
                logging.error('Could not find database directory')
                sys.exit()

        self.repeatdb = os.path.join(self.db, "repeats.fa")

