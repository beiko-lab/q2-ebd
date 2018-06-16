# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Str, Properties, Choices, Int, Bool, Range,
                           Float, Set, Visualization, Metadata, MetadataColumn,
                           Categorical, Citations)

import q2_ebd
from q2_ebd._method import phylogenetic_metrics, non_phylogenetic_metrics
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.tree import Phylogeny, Rooted
from q2_types.ordination import PCoAResults


citations = Citations.load('citations.bib', package='q2_ebd')

plugin = Plugin(
    name='ebd',
    version=q2_ebd.__version__,
    website='https://github.com/beiko-lab/q2-ebd',
    package='q2_ebd',
    description=('This QIIME 2 plugin supports metrics for calculating '
                 'and exploring community alpha and beta diversity through '
                 'statistics and visualizations in the context of sample '
                 'metadata.'),
    short_description='Plugin for exploring community diversity.',
)

plugin.methods.register_function(
    function=q2_ebd.beta_phylogenetic,
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(phylogenetic_metrics()),
                'weighted': Bool},
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.'),
        'phylogeny': ('Phylogenetic tree containing tip identifiers that '
                      'correspond to the feature identifiers in the table. '
                      'This tree can contain tip ids that are not present in '
                      'the table, but all feature ids in the table must be '
                      'present in this tree.')
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.',
        'weighted': 'True if you wish to use the weighted version of the specific measure.'
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity (phylogenetic)',
    description=("Computes a user-specified phylogenetic beta diversity metric"
                 " for all pairs of samples in a feature table."),
    citations=[citations['parks2013measures']]
)


plugin.methods.register_function(
    function=q2_ebd.beta,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'metric': Str % Choices(non_phylogenetic_metrics()),
                'weighted': Bool},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.')
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.',
        'weighted': 'True if you wish to use the weighted version of the specific measure.'
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity',
    description=("Computes a user-specified beta diversity metric for all "
                 "pairs of samples in a feature table."),
    citations=[citations['parks2013measures']]
)

