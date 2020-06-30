import os
import subprocess
import logging
import sys
import re

import pandas as pd

from Bio import pairwise2

class Cluster(object):
    
    def __init__(self, obj):
        self.master = obj

    def run(self):
        '''
        Load the BLAST table and run the different clustering steps
        '''

        # Load blast table and add lengths
        self.df = pd.read_csv(self.master.out+'blast.tab', sep='\t', header=None,
            names=('Repeat', 'Acc', 'Identity', 'Alignment', 'Mismatches', 'Gaps',
                   'Repeat_start', 'Repeat_end', 'Acc_start', 'Acc_end', 'Evalue', 'Score'))
        
        # Calculate coverage
        self.df = self.df.merge(self.master.len_df, left_on='Repeat', right_index=True)
        self.df['Coverage'] = (self.df['Alignment']-self.df['Gaps'])/self.df['Repeat_len']*100

        # Filter by identity and coverage
        self.df = self.df[self.df['Identity'] >= self.master.identity]
        self.df = self.df[self.df['Coverage'] >= self.master.coverage_part]

        # Check if any matches
        if len(self.df) == 0:
            logging.info('No matches with identity >= {}% found'.format(self.master.identity))
            sys.exit()

        # Create new columns
        self.df['Min'] = [min(x,y) for x,y in zip(self.df['Acc_start'],self.df['Acc_end'])]
        self.df['Max'] = [max(x,y) for x,y in zip(self.df['Acc_start'],self.df['Acc_end'])]

        # Keep only best matches if overlapping
        self.remove_overlap()

        # Add repeats
        self.add_repeats()

        # Cluster matches in arrays
        self.cluster_adj()

        # Append partial matches
        self.append_partial()

        # Round
        self.df_appended = self.df_appended.round({'Identity': 1, 'Coverage': 1})
        self.df_appended['Evalue'] = ['{:0.1e}'.format(x) for x in self.df_appended['Evalue']]

        # Split in SRU and arrays
        count_dict = self.df_appended.groupby('Cluster')['Cluster'].count().to_dict()
        cluster_sru = [x for x in count_dict if count_dict[x] == 1]
        cluster_array = [x for x in count_dict if count_dict[x] > 1]
       
        # If any SRUs
        if len(cluster_sru) > 0:
        
            self.df_sru = self.df_appended[self.df_appended['Cluster'].isin(cluster_sru)]
            self.add_flank()
            self.df_sru.to_csv(self.master.out+'SRUs.tab', index=False, sep='\t')
       
        # If any arrays
        if len(cluster_array) > 0:

            self.df_array = self.df_appended[self.df_appended['Cluster'].isin(cluster_array)]
            self.convert_array()
            self.df_arrays.to_csv(self.master.out+'arrays.tab', index=False, sep='\t')

        logging.info('Found {} SRU(s) and {} CRISPR array(s)'.format(len(cluster_sru), len(cluster_array)))

    def overlap(self,x,y):
        '''
        Evaluate whether two (start,end) tuples overlap
        '''
        return x[0] <= y[1] and y[0] <= x[1]

    def overlap_any(self,x,ll):
        '''
        Evaluate whether a (start,end) tuple has overlap with any in a list of tuples
        '''
        return any([self.overlap(x,y) for y in ll])

    def dist(self,x,y):
        '''
        Calculate distance between two (start,end) tuples
        '''
        return y[0]-x[1] if y[0]>x[1] else x[0]-y[1]

    def dist_all(self,x,ll):
        '''
        Get all distances between a (start,end) tuple and a list of tuples
        '''
        return [self.dist(x,y) for y in ll]
    
    def identity(self,x,y):
        '''
        Calculate identity between two sequences
        '''
        align = pairwise2.align.globalxs(x, y, -1, -1, penalize_end_gaps=False)
        return(align[0][2]/min(len(x), len(y))*100)

    def identity_all(self,x,ll):
        '''
        Calculate identity between one sequence and a list of sequences
        '''
        return([self.identity(x, y) for y in ll])

    def remove_overlap(self):
        '''
        If matches overlap keep only the best
        '''

        logging.info('Removing overlapping matches')

        # Sort by alignment quality
        self.df = self.df.sort_values(['Acc', 'Score', 'Coverage'], ascending=False) 

        overlap_lst = []
        for i in set(self.df['Acc']):
            tmp = self.df[self.df['Acc'] == i]
            
            # First remove those with similar start or end
            tmp = tmp.drop_duplicates('Min')
            tmp = tmp.drop_duplicates('Max')
            tmp = tmp.drop_duplicates('Acc_start')
            tmp = tmp.drop_duplicates('Acc_end')

            # Then traverse through matches comparing only with previous
            pos = tmp[['Min','Max']].values
            keep = []
            matches_all = []
            # For each match
            for ind, k in enumerate(pos):
                # If no overlaps with any previous, keep
                if not self.overlap_any(k, matches_all):
                    keep.append(ind)
                    matches_all.append(k)

            overlap_lst.append(tmp.iloc[keep,:])
        
        # If several contigs, concatenate
        self.df_overlap = pd.concat(overlap_lst)

        # Write
        self.df_overlap.to_csv(self.master.out+'blast_best.tab', index=False, sep='\t')
    
    def cluster_adj(self):
        '''
        Cluster adjacent matches into arrays
        '''

        logging.info('Clustering matches')

        # Sort by position
        self.df_overlap = self.df_overlap.sort_values('Min')

        # Split in high and low coverage
        self.df_overlap_compl = self.df_overlap[self.df_overlap['Coverage'] >= self.master.coverage]
        self.df_overlap_part = self.df_overlap[self.df_overlap['Coverage'] < self.master.coverage]

        if len(self.df_overlap_compl) == 0:
            logging.info('No matches with coverage >= {}% found'.format(self.master.coverage))
            sys.exit()

        cluster_df_lst = []
        cluster = 0
        # For each contig
        for i in set(self.df_overlap_compl['Acc']):
            tmp = self.df_overlap_compl[self.df_overlap_compl['Acc'] == i]

            pos = tmp[['Min','Max']].values
            cluster_list = []
            # Loop over complete matches
            for ind, k in enumerate(pos):
                # Keep first match
                if ind == 0:
                    cluster_list.append(cluster)
                    arrays_cluster = [k]
                else:
                    # If match within Xbp of any previous, add match to current cluster
                    if min(self.dist_all(k, arrays_cluster)) <= self.master.max_dist:
                        cluster_list.append(cluster)
                        arrays_cluster.append(k)
                    # If match > Xbp from previous, initiate new cluster
                    else:
                        cluster += 1
                        arrays_cluster = [k]
                        cluster_list.append(cluster)
            
            tmp.insert(len(tmp.columns), 'Cluster', cluster_list)
            cluster_df_lst.append(tmp)
            
            # Increment cluster ID for next acc
            cluster += 1


        # If several contigs, concatenate
        self.df_cluster = pd.concat(cluster_df_lst)

    def append_partial(self):
        '''
        Check if there are any partial matches near clusters
        '''

        append_lst = []
        # For each cluster
        for cl in set(self.df_cluster['Cluster']):
            tmp = self.df_cluster[self.df_cluster['Cluster'] == cl]
            tmp_part = self.df_overlap_part[self.df_overlap_part['Acc'] == list(tmp['Acc'])[0]]

            # Distances between cluster position and partial matches
            cluster_start = min(tmp['Min'])
            cluster_end = max(tmp['Max'])

            dists = self.dist_all((cluster_start, cluster_end), zip(list(tmp_part['Min']),list(tmp_part['Max'])))
            
            # Add partial matches
            part_adj = tmp_part[[x < self.master.max_dist and x > 0 for x in dists]]
            
            # Only those with similar sequences
            idents = part_adj.apply(lambda row: any([k >= self.master.identity for k in self.identity_all(str(row['Sequence']), [str(x) for x in tmp['Sequence'].values])]), axis=1)
            part_adj = part_adj[idents.values] 
            part_adj.insert(len(part_adj.columns), 'Cluster', cl)
            tmp = pd.concat([tmp, part_adj])

            append_lst.append(tmp)

        self.df_appended = pd.concat(append_lst)
        self.df_appended = self.df_appended.sort_values(['Acc', 'Min']) 
        self.df_appended = self.df_appended.drop(columns=['Acc_start','Acc_end'])
        self.df_appended = self.df_appended.rename(columns={'Min':'Start', 'Max':'End'})

    def get_sequence(self, acc, start, end):
        '''
        Return sequence from position information
        '''

        if start < 1:
            start = 1
        
        return(self.master.sequences[str(acc)][(start-1):end])

    def add_repeats(self):
        '''
        Add repeats to the no-overlap dataframe
        '''
        
        self.df_overlap.insert(len(self.df_overlap.columns), 'Sequence', self.df_overlap.apply(lambda row: self.get_sequence(row['Acc'], row['Min'], row['Max']), axis=1))
    
    def add_flank(self):
        '''
        Add flanking sequences to the SRU dataframe
        '''
        
        self.df_sru.insert(len(self.df_sru.columns), 'Left_flank', self.df_sru.apply(lambda row: self.get_sequence(row['Acc'], row['Start']-1-self.master.flank, row['Start']-1), axis=1))
        self.df_sru.insert(len(self.df_sru.columns), 'Right_flank', self.df_sru.apply(lambda row: self.get_sequence(row['Acc'], row['End']+1, row['End']+1+self.master.flank), axis=1))

    def convert_array(self):
        '''
        Add spacers to arrays
        '''

        f = open(self.master.out+'spacers.fa', 'w')

        # For each array
        cls = set(self.df_array['Cluster'])
        dict_lst = []
        for cl in cls:
            tmp = self.df_array[self.df_array['Cluster'] == cl]
            acc = list(tmp['Acc'])[0]
            n = 0

            # Get spacers
            spacers = [str(self.get_sequence(acc, x[0]+1, x[1]-1)) for x in zip(tmp['End'][:(len(tmp)-1)], tmp['Start'][1:])]

            for sp in spacers:
                n += 1
                f.write('>{}_{}:{}\n'.format(acc, cl, n))
                f.write('{}\n'.format(sp))

            # Compile
            dict_lst.append({'Acc': acc,
                            'Start': min(tmp['Start']),
                            'End': max(tmp['End']),
                            'Cluster': cl,
                            'Repeats': [str(x) for x in list(tmp['Sequence'])],
                            'Repeat_types': [re.sub(':.*','',x) for x in list(tmp['Repeat'])],
                            'Spacers': spacers})
        
        self.df_arrays = pd.DataFrame(dict_lst)

        f.close()

        # Add flanks
        self.df_arrays.insert(len(self.df_arrays.columns), 'Left_flank', self.df_arrays.apply(lambda row: self.get_sequence(row['Acc'], row['Start']-1-self.master.flank, row['Start']-1), axis=1))
        self.df_arrays.insert(len(self.df_arrays.columns), 'Right_flank', self.df_arrays.apply(lambda row: self.get_sequence(row['Acc'], row['End']+1, row['End']+1+self.master.flank), axis=1))

