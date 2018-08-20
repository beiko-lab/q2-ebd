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
import sys
import pkg_resources

import holoviews as hv
hv.extension("bokeh")
from bokeh.io import save

import biom
import skbio
import numpy as np
import q2templates

from q2_types.distance_matrix import DistanceMatrixDirectoryFormat

TEMPLATES = pkg_resources.resource_filename('q2_ebd', 'assets')

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

def plot(output_dir: str, distance_matrix: DistanceMatrixDirectoryFormat)-> None:
    measures = []
    n_measures = len(distance_matrix)
    n_samples = None
    i = 0
    scatter_dict = {}
    scatter_opts = dict(height=500, width=500, tools=['hover', 'box_select'], 
                        framewise=True, axiswise=True, size=5)
    for matrix in distance_matrix:
        action_f = open(str(matrix)+"/../provenance/action/action.yaml",'r')
        while (action_f.readline() != "    parameters:\n"):
            pass
        param = action_f.readline()
        md = {}
        while "    -" in param:
            key, val = param.strip()[1:].strip().split(":")
            val = val.strip()
            md[key] = val
            param = action_f.readline()
        measure = md['metric'] + ", " + md['weighted']
        measures.append(measure)
        dm = matrix.file.view(skbio.DistanceMatrix)
        coords = skbio.stats.ordination.pcoa(dm)
        if n_samples is None:
            n_samples = dm.shape[1]
            samples = coords.samples.index
        assert (coords.samples.index == samples).all(), "sample order mismatch, are these all from the same analysis?"
        scatter_dict[measure] = hv.Points(coords.samples[['PC1','PC2']]).options(**scatter_opts) # This is where we'd add in colour
        i += 1
    ds = hv.HoloMap(scatter_dict, kdims=['measure'])
    renderer = hv.renderer('bokeh')
    renderer.save(ds, output_dir + '/ebd')
#    plot = renderer.get_plot(ds).state
#    save(plot, output_dir + '/ebd.html')
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={})

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

def cluster_distance_matrices(dist_mats: DistanceMatrixDirectoryFormat)-> None:
    print(dist_mats)
    pass

