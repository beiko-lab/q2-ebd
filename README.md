# QIIME 2 ExpressBetaDiversity wrapper plugin

This is a QIIME 2 plugin. For details on QIIME 2, see https://qiime2.org.

This QIIME 2 plugin is a wrapper for Donovan Parks' ExpressBetaDiversity program that computes the non-phylogenetic (count-based) beta diversity measures, as well as their phylogenetic equivalents (if a tree is supplied). Its repo is here: https://github.com/dparks1134/ExpressBetaDiversity

This wrapper handles the shuttling between QIIME's internal bits and ExpressBetaDiversity. You just need to have ExpressBetaDiversity built and the binary put somewhere in your PATH.

Activate your QIIME conda environment and `make install` in this directory and this plugin should then be available to your QIIME installation. Commands are available with `qiime ebd`.
