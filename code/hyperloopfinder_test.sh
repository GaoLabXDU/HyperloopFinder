#!/bin/bash

task_suffix=porec_gm12878_25kb_grch38
# Task suffix for naming the folder of output. Users can set a friendly name prefix to distinguish other tasks. This prefix will be combined with other important parameters to set the name of the output directory.

start_step=0
# The task can be started at a user-specified step. This parameter allows the user to start executing the pipeline in the middle step. When the program runs incorrectly, on the premise of ensuring that the previous steps are executed correctly, the user can directly start to execute the error step after modifying the parameters, instead of starting the execution again.
# 0: detect loop
# 1: bin and split cluster file
# 2: mine frequent pattern
# 3: test connectivity by significant loops
# 4: test the significance of hyperloops
# 5: joint all results from different chromosomes


filter_by_loop=1
# Determine whether to use pairwise loops for filtering multi-way chromatin contacts. Our method supports the detection of significant multi-way interactions in all candidate sets, but this is usually not a good choice, because it can cause huge time and resource consumption and even cause computer crashes.
# Accepts: 0 for not, 1 for yes.

correct_bias=1
# Determine whether to correct technical biases. It is recommended to correct the technical biases because the technical biases usually have a great impact on the significance estimation of hyperloops.
# Accepts: 0 for no, 1 for yes

fithic_pass=2
# Required by Fit-Hi-C2, you can use the default values, and users generally do not need to modify it.

# chroms=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
chroms=('chr22')
# Specify the list of chromosomes to be tested for hyperloopFinder.
# Accepts: the names of the chromosomes must have appeared in the input data.

assembly=hg38
# Specify the reference genome version of cluster data.

num_thread=1
# Specify the number of threads that used by HyperloopFinder. Using one thread for one chromosome usually increases the speed of the software operation dramatically. But be careful, it may exceed the limit of the user's machine memory.
# Accepts: integer values greater than or equal to 1.

min_cluster_size=3
# Specify the minimum cluster size to be used in hyperloop mining. The default value is 3.

max_cluster_size=100
# Specify the maximum cluster size to be used in hyperloop mining. Very large clusters may be caused by technical bias. Sometimes it is necessary to filter out large clusters, especially for SPRTIE data.

resolution=25000
# Specify the resolution of bin size for binning interaction loci of multi-way contacts.

min_support_count=5
# Specify the minimum support count parameter for the fpgrowth algorithm. The FP-growth algorithm finds all frequent itemsets whose support is greater than the minimum support count. Here the minimum support count refers to the minimum frequency of multi-way interactions. The minimum support count is greater than or equal to 1. Smaller support thresholds usually produce a large number of results, so be careful with setting small support counts.

max_hyperloop_size=15
# Specify the maximum number of interaction loci included in the hyperloop. This limits the number of sites contained in the maximum hyperloop and can be set as needed.

loop_lowerbound=$((2*$resolution))
# Integer representing lower bound of distance for the intrachromosomal interactions to be considered in base pairs.

loop_upperbound=2000000
# Integer representing upper bound of distance for the intrachromosomal interactions to be considered in base pairs.

loop_fdr_cut=0.05
# False discovery rate threshold for determining whether a pairwise loop detected by Fit-HiC2 is significant.

cluster_file=../data/GM12878/porec_GM12878_Nlalll_grch38_chr22.clusters
# Specify the location of the input file for the multi-way clusters.

chromosome_size_file=../data/${assembly}.chrom.sizes.txt
# Specify the size of each chromosome.

##################################################################################################
# The following variables are automatically generated and generally do not need to be configured.#
##################################################################################################

downweighting=none
# Accepts: none only, this parameter is only passed to get_sprite_contacts.py script. Users can not modify it.

result_dir=../result/hyperloops_results_${resolution}_$task_suffix
# Specify the location and format of the output folder.

loop_dir=${result_dir}/loops
# Specify the output folder for Fit-Hi-C2.

hyperloop_dir=${result_dir}/hyperloops
# Specify the output folder for hyperloops.


echo result dir is $result_dir
if [ ! -d $result_dir ]; then
    mkdir $result_dir
fi

echo loop dir is $loop_dir
if [ ! -d $loop_dir ]; then
    mkdir $loop_dir
fi

echo hyperloop dir is $hyperloop_dir
if [ ! -d $hyperloop_dir ]; then
    mkdir $hyperloop_dir
fi

cp $0 $result_dir/hyperloopfinder.sh

time_start1=$SECONDS

if [ $start_step -eq 0 ] && { [ $filter_by_loop -eq 1 ] || [ $correct_bias -eq 1 ] ;} 
then
    parallel --ungroup --verbose -I% -j $num_thread ./detect_loops.sh % $loop_dir $resolution $assembly $downweighting $cluster_file $max_cluster_size $loop_lowerbound $loop_upperbound $fithic_pass $filter_by_loop $correct_bias ::: ${chroms[*]}
fi
time_end1=$SECONDS
time_spend1=$(($time_end1-$time_start1))
echo $time_spend1 >> $hyperloop_dir/times.txt

python  hyperloopfinder.py       --cluster_file	$cluster_file\
                                 --loop_dir	$loop_dir\
                                 --chromosome_size_file $chromosome_size_file\
                                 --resolution	$resolution\
                                 --min_support_count	$min_support_count\
                                 --min_cluster_size	$min_cluster_size\
                                 --max_cluster_size	$max_cluster_size\
                                 --max_hyperloop_size $max_hyperloop_size\
                                 --loop_lowerbound $loop_lowerbound\
                                 --loop_upperbound $loop_upperbound\
                                 --loop_fdr_cut $loop_fdr_cut\
                                 --output_path	$hyperloop_dir\
                                 --chromosome ${chroms[*]}\
                                 --correct_bias $correct_bias\
                                 --filter_by_loop $filter_by_loop\
                                 --num_thread $num_thread\
                                 --start_step $start_step

# Save the executed script so you can record the parameters.
cp $0 $result_dir/hyperloopfinder.sh
