# HyperloopFinder
An efficient and effective python tool for detecting regulatory multi-way chromatin contacts

## Overview
Recent advances in chromatin conformation capture technologies, such as SPRITE and Pore-C, have enabled the detection of simultaneous contacts among multiple chromatin loci. This has made it possible to investigate the cooperative transcriptional regulation involving multiple genes and regulatory elements at the resolution of a single molecule. However, these technologies are unavoidably subject to the random polymer looping effect and technical biases, making it challenging to distinguish genuine regulatory relationships directly from random polymer interactions.

Therefore, we developed HyperloopFinder, a efficient and effective tool for identifying regulatory multi-way chromatin contacts (hyperloops) by jointly modeling the random polymer looping effect and technical biases to estimate the statistical significance of multi-way contacts.

## Required Package

HyperloopFinder could be installed in a linux-like system. The HyperloopFinder requires the following dependencies. We recommend to use [Anaconda python distribution](https://www.anaconda.com/what-is-anaconda/) for installation of the below python packages.

1. Python (tested on 3.6.13)
2. numpy (tested on 1.18.1)
3. pandas (tested on 1.1.5)
4. matplotlib (tested on 3.1.1)
5. networkx (tested on 2.4)
6. scipy (tested on 1.5.2)
7. fithic (tested on 2.0.7)

HyperloopFinder also require fpgrowth, please download it from https://borgelt.net/fpgrowth.html and put it in HyperloopFinder/utils directory.

## Installation

Download HyperloopFinder by

```shell
git clone https://github.com/gaolabXDU/HyperloopFinder
```

## Usage

In order to facilitate the use of the software, the user needs to configure the hyperloopfinder.sh script in the code folder with the following main configured parameters:


### task_suffix
Task suffix for naming the folder of output.

### start_step
Task can be start at a user-specified step.
Accepts:
0: detect loop
1: bin and split cluster file
2: mine frequent pattern
3: test connectivity by significant loops
4: test the significance of hyperloops
5: joint all results from different chromosomes

### filter_by_loop = 1
Determine whether to use pairwise loops for filtering multi-way chromatin contacts.

Accepts: 0 for not, 1 for yes.

### correct_bias = 1
Determine whether to correct technical biases.

Accepts: 0 for no, 1 for yes

### fithic_pass = 2
Determine the number of spline passes of fithic software.

Accepts: integer values greater than or equal to 1.


### chroms=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')

Specify the list of chromosomes to be tested for hyperloop.

Accepts: the names of the chromosomes must have appeared in the input data.

### assembly = hg38
Specify the reference genome version of input data.

### downweighting = none
Accepts: none only, this parameter is only passed to get_sprite_contacts.py script.

### num_thread = 1
Specify the number of threads that used by HyperloopFinder. Using one thread for a chromosome usually increases the speed of the software operation dramatically.

Accepts: integer values greater than or equal to 1.

### min_cluster_size = 3

Specify the minimum cluster size to be used in hyperloop mining.

### max_cluster_size = 100
Specify the maximum cluster size to be used in hyperloop mining.

### resolution = 25000
Specify the resolution of bin size for binning interaction loci of multi-way contacts.

### min_support_count=5
Specify the minimum support count parameter for the fpgrowth algorithm. 

### max_hyperloop_size=3
Specify the maximum number of interaction loci included in the hyperloop.

### loop_lowerbound=25000

Integer representing lower bound for the intrachromosomal interactions to be considered in base pairs.

### loop_upperbound=10000000
Integer representing upper bound for the intrachromosomal interactions to be considered in base pairs.

### loop_fdr_cut=0.05
False discovery rate threshold for determining whether a pairwise loop detected by fithic is significant.

### cluster_file=../data/GM12878/human.combined.mapq-ge10_splited_up10mb.clusters
Specify the location of the input file for the multi-way clusters.

### result_dir=../result/GM12878/hyperloops_results_${resolution}_${downweighting}_$task_suffix
Specify the location and format of the output folder.

### loop_dir=${result_dir}/loops
Specify the output folder for fithic.

### hyperloop_dir=${result_dir}/hyperloops
Specify the output folder for HyperloopFinder.

### chromosome_size_file=../data/hg38.chrom.sizes.txt
Specify the size of each chromosome.
