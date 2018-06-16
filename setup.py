# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer

setup(
    name="q2-ebd",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={'q2_ebd': ['citations.bib'],
                  'q2_ebd.tests': [
                      'data/*'
                  ]},
    author="Michael Hall",
    author_email="hallm2533@gmail.com",
    description="ExpressBetaDiversity wrapper.",
    license='BSD-3-Clause',
    url="https://qiime2.org",
    entry_points={
        'qiime2.plugins': ['q2-ebd=q2_ebd.plugin_setup:plugin']
    },
    zip_safe=False,
)
