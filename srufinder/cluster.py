import os
import subprocess
import logging
import sys

import pandas as pd

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

        # Cluster matches in arrays
        self.cluster_adj()

        # Append partial matches
        self.append_partial()

        # Add sequences
        self.add_sequences()

        # Round
        self.df_appended = self.df_appended.round({'Identity': 1, 'Coverage': 1})
        self.df_appended['Evalue'] = ['{:0.1e}'.format(x) for x in self.df_appended['Evalue']]

        # Write result
        self.df_appended.to_csv(self.master.out+'clusters.tab', index=False, sep='\t')

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
    
    def remove_overlap(self):
        '''
        If matches overlap keep only the best
        '''

        logging.info('Removing overlapping matches')

        # Sort by alignment quality
        self.df = self.df.sort_values(['Acc', 'Score', 'Alignment'], ascending=False) 

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
        cluster_lst = []
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

            # Distances between cluster position and partial matches
            cluster_start = min(tmp['Min'])
            cluster_end = max(tmp['Max'])

            dists = self.dist_all((cluster_start, cluster_end), zip(list(self.df_overlap_part['Min']),list(self.df_overlap_part['Max'])))
            
            # Add partial matches
            part_adj = self.df_overlap_part[[x < self.master.max_dist and x > 0 for x in dists]]
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

        return(self.master.sequences[acc][(start-1):end])

    def add_sequences(self):
        '''
        Add sequences to the clustering dataframe
        '''
        
        self.df_appended['Sequence'] = self.df_appended.apply(lambda row: self.get_sequence(row['Acc'], row['Start'], row['End']), axis=1)
