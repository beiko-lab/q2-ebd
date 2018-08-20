#!/usr/bin/bash


table_file=filtered_table.qza
metadata_file=METADATA.txt
tree_file=rooted_tree.qza
mkdir distance_matrices
#This will throw a bunch of errors for distance measures that can't/shouldn't exist (like a nonphylogenetic UniFrac, for example), but just ignore them
for metric in 'braycurtis' 'canberra' 'chi_squared' 'coeff_similarity' 'complete_tree' 'euclidean' 'f_st' 'p_st' 'gower' 'hellinger' 'kulczynski' 'lennon' 'manhattan' 'mnnd' 'mpd' 'morisita_horn' 'normalized_weighted_unifrac' 'pearson' 'raohp' 'soergel' 'tamas_coeff' 'weighted_corr' 'whittaker' 'yue_clayton'; do

	qiime ebd beta --i-table ${table_file} --p-weighted --p-metric ${metric} --o-distance-matrix distance_matrices/${metric}_nonphylogenetic_weighted.qza
	qiime ebd beta --i-table ${table_file} --p-no-weighted --p-metric ${metric} --o-distance-matrix distance_matrices/${metric}_nonphylogenetic_unweighted.qza
	qiime ebd beta-phylogenetic --i-table ${table_file} --i-phylogeny ${tree_file} --p-weighted --p-metric ${metric} --o-distance-matrix distance_matrices/${metric}_phylogenetic_weighted.qza
	qiime ebd beta-phylogenetic --i-table ${table_file} --i-phylogeny ${tree_file} --p-no-weighted --p-metric ${metric} --o-distance-matrix distance_matrices/${metric}_phylogenetic_unweighted.qza

        qiime diversity pcoa --i-distance-matrix distance_matrices/${metric}_phylogenetic_unweighted.qza --o-pcoa distance_matrices/${metric}_phylogenetic_unweighted_pcoa.qza
        qiime diversity pcoa --i-distance-matrix distance_matrices/${metric}_phylogenetic_weighted.qza --o-pcoa distance_matrices/${metric}_phylogenetic_weighted_pcoa.qza
        qiime diversity pcoa --i-distance-matrix distance_matrices/${metric}_nonphylogenetic_unweighted.qza --o-pcoa distance_matrices/${metric}_nonphylogenetic_unweighted_pcoa.qza
        qiime diversity pcoa --i-distance-matrix distance_matrices/${metric}_nonphylogenetic_weighted.qza --o-pcoa distance_matrices/${metric}_nonphylogenetic_weighted_pcoa.qza

        qiime emperor plot --i-pcoa distance_matrices/${metric}_phylogenetic_unweighted_pcoa.qza --o-visualization distance_matrices/${metric}_phylogenetic_unweighted.qzv --m-metadata-file ${metadata_file}
        qiime emperor plot --i-pcoa distance_matrices/${metric}_phylogenetic_weighted_pcoa.qza --o-visualization distance_matrices/${metric}_phylogenetic_weighted.qzv --m-metadata-file ${metadata_file}
        qiime emperor plot --i-pcoa distance_matrices/${metric}_nonphylogenetic_unweighted_pcoa.qza --o-visualization distance_matrices/${metric}_nonphylogenetic_unweighted.qzv --m-metadata-file ${metadata_file}
        qiime emperor plot --i-pcoa distance_matrices/${metric}_nonphylogenetic_weighted_pcoa.qza --o-visualization distance_matrices/${metric}_nonphylogenetic_weighted.qzv --m-metadata-file ${metadata_file}

done
