# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import tempfile
import os
import subprocess

import biom
import skbio
import numpy as np

def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)

# We should consider moving these functions to scikit-bio. They're part of
# the private API here for now.
def phylogenetic_metrics():
    return {'unweighted_unifrac', 'weighted_unifrac', 'normalized_weighted_unifrac',
            'braycurtis', 'sorensen', 'canberra', 'chi_squared', 'coeff_similarity',
            'complete_tree', 'euclidean', 'p_st', 'gower', 'hellinger', 'kulczynski',
            'lennon', 'manhattan', 'mnnd', 'mpd', 'morisita_horn', 'pearson', 'raohp',
            'soergel', 'jaccard', 'ruzicka', 'tamas_coeff', 'weighted_corr', 
            'whittaker', 'yue_clayton'}

def non_phylogenetic_metrics():
    return {'braycurtis', 'sorensen', 'canberra', 'chi_squared', 'coeff_similarity',
            'euclidean', 'f_st', 'gower', 'hellinger', 'kulczynski',
            'lennon', 'manhattan', 'morisita_horn', 'pearson', 'raohp',
            'soergel', 'jaccard', 'ruzicka', 'tamas_coeff', 'weighted_corr', 
            'whittaker', 'yue_clayton'}

def all_metrics():
    return phylogenetic_metrics() | non_phylogenetic_metrics()


def beta_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                      metric: str, weighted: bool)-> skbio.DistanceMatrix:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)
    if table.is_empty():
        raise ValueError("The provided table object is empty")

    # Write table to temp file
    with tempfile.TemporaryDirectory() as temp_dir_name:
        table_fp = os.path.join(temp_dir_name, 'otu_table.tsv')
        newick_fp = os.path.join(temp_dir_name, 'tree.newick')
        with open(table_fp, 'w') as out_table, open(newick_fp, 'w') as newick:
            # This is easy, just write to newick
            phylogeny.write(newick)
            # We have to iterate through each sample
            out_table.write("\t" + "\t".join(table.ids(axis='observation')))
            for sample_id in table.ids(axis='sample'):
                row = table.data(sample_id)
                out_table.write("\n" + str(sample_id) + "\t" + \
                        "\t".join([str(x) for x in row]))
    # Run ExpressBetaDiversity on them
        name_map = {'braycurtis': 'Bray-Curtis',
                    'sorensen': 'Bray-Curtis',
                    'canberra': 'Canberra',
                    'chi_squared': 'Chi-squared',
                    'coeff_similarity': 'CS',
                    'complete_tree': 'CT',
                    'euclidean': 'Euclidean',
                    'f_st': 'Fst',
                    'p_st': 'Fst',
                    'gower': 'Gower',
                    'hellinger': 'Hellinger',
                    'kulczynski': 'Kulczynski',
                    'lennon': 'Lennon',
                    'manhattan': 'Manhattan',
                    'weighted_unifrac': 'Manhattan',
                    'mnnd': 'MNND',
                    'mpd': 'MPD',
                    'morisita_horn': 'Morisita-Horn',
                    'normalized_weighted_unifrac': 'NWU',
                    'pearson': 'Pearson',
                    'raohp': 'RaoHp',
                    'soergel': 'Soergel',
                    'jaccard': 'Soergel',
                    'unweighted_unifrac': 'Soergel',
                    'ruzicka': 'Soergel',
                    'tamas_coeff': 'TC',
                    'weighted_corr': 'WC',
                    'whittaker': 'Whittaker',
                    'yue_clayton': 'Yue-Clayton'
                   }
        if weighted:
            weighted = "-w"
        else:
            weighted = ""
        cmd = 'ExpressBetaDiversity -t tree.newick -s otu_table.tsv %s -c %s' \
                                                  % (weighted, name_map[metric])
        subprocess.run(cmd, cwd=temp_dir_name, shell=True)
        with open(os.path.join(temp_dir_name, 'output.diss'), 'r') as dist_file:
            nsamples = int(dist_file.readline())
            dist_mat = np.zeros((nsamples, nsamples))
            ids = []
            for i, line in enumerate(dist_file):
                ids.append(line.split("\t")[0].strip())
                for j, dist in enumerate(line.split("\t")[1:]):
                    dist_mat[i,j] = float(dist)
                    dist_mat[j,i] = float(dist)

    # Suck the data matrix back in
    # Return a DistanceMatrix object
    results = skbio.DistanceMatrix(dist_mat, ids)
    return results



def beta(table: biom.Table, metric: str, weighted: bool)-> skbio.DistanceMatrix:
    if metric not in non_phylogenetic_metrics():
        raise ValueError("Unknown metric: %s" % metric)
    if table.is_empty():
        raise ValueError("The provided table object is empty")

    # Write table to temp file
    with tempfile.TemporaryDirectory() as temp_dir_name:
        table_fp = os.path.join(temp_dir_name, 'otu_table.tsv')
        with open(table_fp, 'w') as out_table:
            # We have to iterate through each sample
            out_table.write("\t" + "\t".join(table.ids(axis='observation')))
            for sample_id in table.ids(axis='sample'):
                row = table.data(sample_id)
                out_table.write("\n" + str(sample_id) + "\t" + \
                        "\t".join([str(x) for x in row]))
    # Run ExpressBetaDiversity on them
        name_map = {'braycurtis': 'Bray-Curtis',
                    'sorensen': 'Bray-Curtis',
                    'canberra': 'Canberra',
                    'chi_squared': 'Chi-squared',
                    'coeff_similarity': 'CS',
                    'euclidean': 'Euclidean',
                    'f_st': 'Fst',
                    'gower': 'Gower',
                    'hellinger': 'Hellinger',
                    'kulczynski': 'Kulczynski',
                    'lennon': 'Lennon',
                    'manhattan': 'Manhattan',
                    'morisita_horn': 'Morisita-Horn',
                    'pearson': 'Pearson',
                    'raohp': 'RaoHp',
                    'soergel': 'Soergel',
                    'jaccard': 'Soergel',
                    'ruzicka': 'Soergel',
                    'tamas_coeff': 'TC',
                    'weighted_corr': 'WC',
                    'whittaker': 'Whittaker',
                    'yue_clayton': 'Yue-Clayton'
                   }
        if weighted:
            weighted = "-w"
        else:
            weighted = ""
        cmd = 'ExpressBetaDiversity -s otu_table.tsv %s -c %s' \
                                          % (weighted, name_map[metric])
        subprocess.run(cmd, cwd=temp_dir_name, shell=True)
        with open(os.path.join(temp_dir_name, 'output.diss'), 'r') as dist_file:
            nsamples = int(dist_file.readline())
            dist_mat = np.zeros((nsamples, nsamples))
            ids = []
            for i, line in enumerate(dist_file):
                ids.append(line.split("\t")[0].strip())
                for j, dist in enumerate(line.split("\t")[1:]):
                    dist_mat[i,j] = float(dist)
                    dist_mat[j,i] = float(dist)

    # Suck the data matrix back in
    # Return a DistanceMatrix object
    results = skbio.DistanceMatrix(dist_mat, ids)

    return results
