# Simulation Branch of MosaiCatcher

This branch is only meant to test the performance of MC on simulated SV calls.

## Usage:

This pipeline automatically checks out the [Mosaicatcher pipeline](https://github.com/friendsofstrandseq/pipeline)
hence uses the latest updates during calling.

However, the Snakeamke rules for **segmentation** and **SV classification** might need to be
updated in this `Snakefile`.

## TODO

 * currently the evaluation is done only on the position of loci, irrespective of
   the calls being made in single cells. This should be added later

