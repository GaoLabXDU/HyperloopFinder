#!/bin/bash
task_suffix=test

start_step=0
filter_by_loop=1
correct_bias=1
fithic_pass=2
# 0: detect loop
# 1: bin and split cluster file
# 2: mine frequent pattern
# 3: test connectivity by significant loops
# 4: test the significance of hyperloops
# 5: joint all results from different chromosomes

#chroms=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
chroms=('chr21')

assembly=hg38
downweighting=none
num_thread=6

min_cluster_size=3
max_cluster_size=100
resolution=25000
min_support_count=5
max_hyperloop_size=15
loop_lowerbound=25000
loop_upperbound=10000000
loop_fdr_cut=0.05

cluster_file=../data/GM12878/random_shuffled_porec_clusters/porec_GM12878_Nlalll_grch38_random_shuffled_repeat0.clusters
result_dir=../result/GM12878/hyperloops_results_${resolution}_${downweighting}_$task_suffix
loop_dir=${result_dir}/loops
hyperloop_dir=${result_dir}/hyperloops
chromosome_size_file=../data/hg38.chrom.sizes.txt


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

cp $0 $result_dir/hyperloopfinder.sh
