#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import networkx as nx
from scipy.special import comb
from itertools import combinations
import scipy as sp
from multiprocessing import Pool as ThreadPool
import time
import argparse
from scipy.stats import binomtest


# ============================
# parse input arguments
# ============================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Detect hyperloops from multi-way chromatin interaction data.')

    parser.add_argument('--cluster_file',
                        metavar="FILE",
                        action="store",
                        help="File path of cluster file.")

    parser.add_argument('--loop_dir',
                        metavar="FILE",
                        action="store",
                        help="File path of chromatin loops")

    parser.add_argument('--chromosome_size_file',
                        metavar="FILE",
                        action="store",
                        help="File path of chromosome size")

    parser.add_argument('--external_loop_path',
                        metavar="FILE",
                        action="store",
                        default='void',
                        help="File path of external loops")

    parser.add_argument('--filter_by_loop',
                        metavar="INT",
                        action="store",
                        type=int,
                        default=1,
                        help="Whether to filter hyperloops by significant loops? 1 is yes, 0 is no")

    parser.add_argument('--correct_bias',
                        metavar="INT",
                        action="store",
                        type=int,
                        default=1,
                        help="Whether to correct bias? 1 is yes, 0 is no")

    parser.add_argument('-o', '--output_path',
                        metavar="FILE",
                        action="store",
                        help="Path of outputs.")

    parser.add_argument('--min_support_count',
                        metavar="INT",
                        type=int,
                        action="store",
                        default=5,
                        help="Minimum support count for FPGrowth algorithm.")

    parser.add_argument('--max_cluster_size',
                        metavar='MAX',
                        type=int,
                        action='store',
                        default=1000,
                        help="Maximum cluster size of multi-way chromatin interaction data.")

    parser.add_argument('--max_hyperloop_size',
                        metavar='MAX',
                        type=int,
                        action='store',
                        default=10,
                        help="Maximum hyperloop size that user defined.")

    parser.add_argument('--loop_lowerbound',
                        metavar='MAX',
                        type=int,
                        action='store',
                        default=0,
                        help="Minimum hyperloop distance that user defined.")

    parser.add_argument('--loop_upperbound',
                        metavar='INT',
                        type=int,
                        action='store',
                        default=-1,
                        help="Maximum hyperloop distance that user defined.")

    parser.add_argument('--loop_fdr_cut',
                        metavar='FLOAT',
                        type=float,
                        action='store',
                        default=0.1,
                        help="Q value cut for significant pairwise loops.")

    parser.add_argument('--min_cluster_size',
                        metavar='INT',
                        type=int,
                        action='store',
                        default=3,
                        help="Minimum cluster size of multi-way chromatin interaction data.")

    parser.add_argument('--resolution',
                        metavar='INT',
                        type=int,
                        action='store',
                        default=40000,
                        help="Contact matrix resolution in bp.")

    parser.add_argument('--num_thread',
                        metavar='INT',
                        type=int,
                        action='store',
                        default=1,
                        help="Number of CPU threads.")

    parser.add_argument('--chromosome',
                        metavar='List',
                        nargs='+',
                        default='all',
                        help="Chromosomes.")

    parser.add_argument('--start_step',
                        metavar='INT',
                        type=int,
                        action='store',
                        default=1,
                        help="From which step does the program start?")

    return parser.parse_args()


def main():
    args = parse_arguments()
    resolution = args.resolution
    min_cluster_size = args.min_cluster_size
    max_cluster_size = args.max_cluster_size
    max_hyper_loop = args.max_hyperloop_size
    loop_lowerbound = args.loop_lowerbound
    loop_upperbound = args.loop_upperbound
    loop_fdr_cut = args.loop_fdr_cut
    min_sup = args.min_support_count
    cluster_file = args.cluster_file
    output_path = args.output_path
    loop_dir = args.loop_dir
    external_loop_path = args.external_loop_path
    num_thread = args.num_thread
    start_step = args.start_step
    filter_by_loop = args.filter_by_loop
    correct_bias = args.correct_bias
    chromosome_size_file = args.chromosome_size_file

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if args.chromosome[0] == 'all':
        chroms = ['chr'+str(i)for i in range(1, 23)]+['chrX']
    else:
        chroms = args.chromosome

    chromosome_size = pd.read_csv(chromosome_size_file, sep='\t', names=[
                                  'chrom', 'length'], index_col='chrom')

    params = [(chrom, chromosome_size.at[chrom, 'length'], min_sup, output_path, loop_dir, external_loop_path,
               resolution, max_hyper_loop, loop_lowerbound, loop_upperbound, loop_fdr_cut, start_step, filter_by_loop, correct_bias) for chrom in chroms]

    time_start = time.time()
    if start_step <= 1:
        print('Step 1, bin and split cluster file................')
        bin_clusters(cluster_file, output_path, min_cluster_size,
                     max_cluster_size, resolution, chroms)
        time_get_bins = time.time()
        print('Step 1 done.')
        print("Time cost in getting bins = {}s".format(time_get_bins-time_start))

    time_end1=time.time()
    with open(output_path+'/times.txt',mode='a') as tf:
        tf.write('{}'.format(time_end1-time_start))

    if start_step <= 4:
        pool = ThreadPool(processes=num_thread)
        results = pool.starmap(detect_by_chrom, params)
        pool.close()
        pool.join()

    # Joint all results from different chromosomes.
    if start_step <= 5:
        print('Step 5, joint all results from different chromosomes................')
        join_all_result_file(output_path, chroms, filter_by_loop)
        print('Step 5 done.')

    time_end = time.time()
    # print(results)
    print("Thread number = {},Time cost = {}s".format(
        num_thread, time_end-time_start))

# # Binning
def bin_cluster_on_line(line, clique_dict, min_cluster_size, max_cluster_size, resolution):
    items = line.strip().split('\t')
    if not min_cluster_size <= len(items)-1 <= max_cluster_size:
        return
    items_split = [(item.split(':')[0], int(item.split(':')[1]))
                   for item in items[1:]]
    chr_dict = dict()
    for item in items_split:
        bin_index = item[1]//resolution
        mid_point = bin_index*resolution+resolution//2
        if item[0] in chr_dict:
            chr_dict[item[0]].append(mid_point)
        else:
            chr_dict[item[0]] = [mid_point]
    for key in chr_dict:
        data = list(set(chr_dict[key]))
        if len(data) >= 3:
            if key in clique_dict:
                clique_dict[key].append(data)
            else:
                clique_dict[key] = [data]


def bin_clusters(cluster_file, output_path, min_cluster_size, max_cluster_size, resolution, chroms):
    clique_dict = dict()
    file_object = open(cluster_file, 'r')
    try:
        i = 0
        for line in file_object:
            bin_cluster_on_line(line, clique_dict,
                                min_cluster_size, max_cluster_size, resolution)
            # i+=1
            # if i==1000:
            #    break
    finally:
        file_object.close()

    for key in chroms:
        with open(output_path+'/'+key+'.clusters', 'w') as file_output:
            for cluster in clique_dict[key]:
                str_out = [str(item) for item in cluster]
                file_output.write('\t'.join(str_out)+'\n')


# # Mining frequent pattern
def detect_frequent_pattern(min_sup, max_hyper_loop, output_path, chrom, filter_by_loop):
    if filter_by_loop:
        node_file = output_path+'/{}_nodes.txt'.format(chrom)
        os.system('../utils/fpgrowth -ts -m3 -n{} -s-{} -q0 -Z -k " " -v " %a" -P {}/{}_pattern_spectrum_screen.txt -R {} {}/{}.clusters {}/{}_screen.fp'.format(
            max_hyper_loop, min_sup, output_path, chrom, node_file, output_path, chrom, output_path, chrom))
    else:
        os.system('../utils/fpgrowth -ts -m3 -n{} -s-{} -q0 -Z -k " " -v " %a" -P {}/{}_pattern_spectrum.txt {}/{}.clusters {}/{}.fp'.format(
            max_hyper_loop, min_sup, output_path, chrom, output_path, chrom, output_path, chrom))


# # Testing conectivity
def is_connected(bins, loop_set):
    g = nx.Graph()
    g.add_nodes_from(bins)
    for i in range(len(bins)):
        for j in range(len(bins)):
            edge = (bins[i], bins[j])
            if edge in loop_set:
                g.add_edge(*edge)
    return nx.is_connected(g)


def detect_connected_pattern(loop_dir, external_loop_path, chrom, resolution, loop_fdr_cut, output_path):
    if external_loop_path == 'void':
        loop_file = '{}/{}_fithic/FitHiC.spline_pass2.res{}.significances.txt.gz'.format(
            loop_dir, chrom, resolution)
        loop_file_new='{}/{}_fithic/FitHiC.spline_pass2.res{}.significances.new.txt.gz'.format(
            loop_dir, chrom, resolution)
        loop_df = pd.read_csv(loop_file, sep='\t')

        bias_file = '{}/{}_bias_{}.tsv.gz'.format(loop_dir, chrom, resolution)
        bias_pd = pd.read_csv(bias_file, sep='\t', header=None, names=[
                            'chrom', 'loc', 'bias_index'])
        bin_bias = bias_pd.bias_index.values
        group_num = 50

        loop_df = loop_df.sort_values(by='p-value')
        loop_df['distance'] = loop_df.fragmentMid2-loop_df.fragmentMid1

        loop_df['q_value_new'] = 1.0
        hp_dist_group, cuts = quantile_cut(loop_df.distance, group_num)
        cuts_bins = (cuts+resolution)//resolution
        grouped_test_nums_k_list = []
        for dist_limit in cuts_bins:
            grouped_test_nums_k_list.append(
                get_test_num(bin_bias, dist_limit, 2, -1))

        grouped_test_nums_k = np.array(grouped_test_nums_k_list)
        grouped_test_nums_k = grouped_test_nums_k[1:]-grouped_test_nums_k[:-1]
        for dist_group_index, df_k_d in loop_df.groupby(hp_dist_group):
            qvals = get_qvalues(df_k_d['p-value'],
                                grouped_test_nums_k[dist_group_index])
            loop_df.loc[df_k_d.index, 'q_value_new'] = qvals
        loop_df['q-value']=loop_df['q_value_new']

        #loop_significant=loop_df.query("q_value <= 0.05")
        loop_df=loop_df.drop(columns=['q_value_new'])
        loop_df.to_csv(loop_file_new, sep='\t', index=False, compression='gzip')
    else:
        loop_df = pd.read_csv(external_loop_path, sep='\t')
        loop_df = loop_df[loop_df.chr1 == chrom]
        #loop_significant = loop_df.query("q-value <= {}".format(loop_fdr_cut))
    
    loop_significant = loop_df[loop_df['q-value'] <= loop_fdr_cut]

    hyperloop_file = '{}/{}_screen.fp'.format(output_path, chrom)
    hyperloop_connected_file = '{}/{}_connected_hyperloop.csv'.format(
        output_path, chrom)

    loop_set = set(loop_significant.apply(lambda row: (
        row.fragmentMid1, row.fragmentMid2), axis=1))

    scf = open(hyperloop_connected_file, 'w')
    i = 0
    with open(hyperloop_file) as sf:
        for line in sf:
            items = line.strip().split()
            count = int(items[-1])
            bins = items[:-1]
            bins = [int(x)for x in bins]
            if is_connected(bins, loop_set):
                scf.write(line)
    scf.close()


def quantile_cut(x, q):
    x_np = np.array(x)
    # Ensure that there are at least 100 samples in each group.
    q_upper = x_np.shape[0]//100
    if q_upper == 0:
        return np.zeros_like(x_np), np.array([np.min(x_np), np.max(x_np)])
    q_last = min(q, q_upper)
    quantiles = np.linspace(0, 1, q_last+1)
    bins = np.quantile(x_np, quantiles, interpolation='higher')
    bins = np.unique(bins)
    ids = bins.searchsorted(x_np, side='right')
    q_real = bins.shape[0]-1
    ids[ids == q_real+1] = q_real
    ids = ids-1
    return ids, bins

# # Testing significance
def detect_significant_pattern(output_path, chrom, chr_len, loop_dir, external_loop_path, resolution, filter_by_loop, max_hyper_loop, loop_lowerbound, loop_upperbound, correct_bias):
    cluster_file = '{}/{}.clusters'.format(output_path, chrom)

    if filter_by_loop:
        hyperloop_file = '{}/{}_connected_hyperloop.csv'.format(
            output_path, chrom)
    else:
        hyperloop_file = '{}/{}.fp'.format(output_path, chrom)

    contact_file = '{}/{}_contacts_{}.tsv.gz'.format(
        loop_dir, chrom, resolution)

    if filter_by_loop:
        hyperloop_out_file = '{}/{}_connected_significant_hyperloop.csv'.format(
            output_path, chrom)
    else:
        hyperloop_out_file = '{}/{}_significant_hyperloop.csv'.format(
            output_path, chrom)

    chr_len_bin = chr_len//resolution+1

    bin_count, num_trials_array, loop_distribs_array = get_pairwise_loop_distribution(
        chr_len_bin, max_hyper_loop, cluster_file, resolution)

    cut_value = max(np.mean(bin_count)-2*np.std(bin_count), 0)
    bin_mask = (bin_count > cut_value)

    # Load the bias file
    if correct_bias:
        bias_file = loop_dir+'/{}_bias_{}.tsv.gz'.format(chrom, resolution)
        bias_df = pd.read_csv(bias_file, sep='\t', header=None, names=[
                              'chr', 'mid_point', 'bias'], index_col='mid_point')
        bias_series = bias_df.bias.copy()
        # bias_upper_bound=2
        # bias_lower_bound=0.5
        # bias_series.loc[bias_series<bias_lower_bound]=bias_lower_bound
        # bias_series.loc[bias_series>bias_upper_bound]=bias_upper_bound
        # bias=${loop_dir}/${chrom}_bias_${resolution}.tsv.gz
    else:
        bias_series = None

    processed_lines = 0
    hyperloops = []
    if loop_upperbound == -1:
        loop_upperbound = chr_len
    with open(hyperloop_file) as hf:
        for line in hf:
            items = line.strip().split()
            count = int(items[-1])
            bins = [int(item) for item in items[:-1]]
            bins = sorted(bins)
            if bins[-1] > chr_len:
                continue
            size = len(bins)
            distance = bins[-1]-bins[0]
            if distance < loop_lowerbound or distance > loop_upperbound:
                continue
            # loc_num=np.sum(bin_count>0)-(distance//resolution+1)
            # loc_num=chr_len_bin-(distance//resolution+1)
            relative_dist_bin = [(bin-bins[0])//resolution for bin in bins]
            loop_len_bin = distance//resolution+1
            candidates = np.arange(0, chr_len_bin-loop_len_bin+1).reshape(-1,
                                                                          1)+np.array(relative_dist_bin).reshape((1, -1))
            candidates_mask = np.array(
                [bin_mask[candidates[i, :]] for i in range(candidates.shape[0])])
            loc_num = np.sum(np.product(candidates_mask, axis=1))
            if loc_num == 0:
                continue

            # is_valid=True
            # for item in distance:
                # if item >=min_dist and item <=max_dist:
                # continue
                # else:
                # is_valid=False
            # if not is_valid:
                # continue
            distrib = loop_distribs_array[:, len(bins)-3].ravel()
            num_trial = int(num_trials_array[len(bins)])
            prob, expectation, total_bias, is_valid = get_expectation(
                distrib, bins, loc_num, correct_bias, bias_series, resolution, num_trial)
            if is_valid:
                fold_enrichment = count/expectation
                hyperloops.append([chrom, bins, count, size, distance, prob,
                                  num_trial, total_bias, expectation, fold_enrichment])
            processed_lines = processed_lines+1
            if processed_lines % 100000 == 0:
                print('{} lines have processed...'.format(processed_lines))
                # break

    hyperloop_df = pd.DataFrame(hyperloops, columns=[
                                'chrom', 'hyperloop', 'count', 'size', 'distance', 'prob', 'num_trial', 'total_bias', 'expectation', 'fold_enrichment'])

    group_num = 50

    def get_pvalue(row):
        k = row['count']
        n = row['num_trial']
        p = row['prob']
        return binomtest(k, n, p, alternative='greater').pvalue

    pvalues = hyperloop_df.apply(get_pvalue, axis=1)
    hyperloop_df['p_value'] = pvalues
    hyperloop_df = hyperloop_df.sort_values(by='p_value')
    hyperloop_df.reset_index(drop=True, inplace=True)

    hyperloop_df['q_value'] = 1.0
    for k, df_k in hyperloop_df.groupby('size'):
        hp_dist_group, cuts = quantile_cut(df_k.distance, group_num)
        cuts_bins = (cuts+resolution)//resolution
        # Correct the maximum and minimum values
        cuts_bins[0] = loop_lowerbound//resolution
        cuts_bins[-1] = chr_len_bin
        grouped_test_nums_k_list = []
        for dist_limit in cuts_bins:
            grouped_test_nums_k_list.append(
                get_test_num(bin_count, dist_limit, k, 0))

        grouped_test_nums_k = np.array(grouped_test_nums_k_list)
        grouped_test_nums_k = grouped_test_nums_k[1:]-grouped_test_nums_k[:-1]

        for dist_group_index, df_k_d in df_k.groupby(hp_dist_group):
            qvals = get_qvalues(
                df_k_d.p_value, grouped_test_nums_k[dist_group_index])
            hyperloop_df.loc[df_k_d.index, 'q_value'] = qvals

    hyperloop_df.to_csv(hyperloop_out_file, sep=',', index=False)


def get_pairwise_loop_distribution(chr_len_bin, max_hyper_loop, cluster_file, resolution):
    comb_counts = np.zeros(max_hyper_loop+1)
    bin_count = np.zeros((chr_len_bin,))
    D3 = np.zeros((chr_len_bin, max_hyper_loop-3+1))
    i = 0
    with open(cluster_file) as clusters:
        for line in clusters:
            bins = line.strip().split()
            bins = [int(bin) for bin in bins]
            bins = sorted(bins)
            size = len(bins)
            bin_count[np.array(bins)//resolution] += 1

            # for k in range(3,min(max_hyper_loop,size)+1):
            #     comb_counts[k]+=comb(size,k)
            #     #comb_counts[k]+=enum_kmer(bins,k,min_dist,max_dist)
            k_vec = np.arange(3, min(max_hyper_loop, size)+1)
            comb_counts[k_vec] += comb(size, k_vec)

            pairs_index = np.array([*combinations(range(size), 2)])
            pairs = np.array([*combinations(bins, 2)])
            # pairs=np.choose(pairs_index,bins)
            d_index = pairs_index[:, 1].ravel()-pairs_index[:, 0].ravel()
            d_genome = (pairs[:, 1].ravel()-pairs[:, 0].ravel())//resolution

            # loop version
            # for k in range(3,min(15,size)+1):
            #     weight=comb(size-(d_index+1),k-2)
            #     weight_series=pd.Series(data=weight,index=d_genome)
            #     weight_sum_series=weight_series.groupby(level=0).sum()
            #     D3[weight_sum_series.index,k-3]+=weight_sum_series.values
            N_comb = (size-(d_index+1)).reshape((-1, 1))
            k_comb = k_vec-2
            weight = comb(N_comb, k_comb)
            weight_pd = pd.DataFrame(weight, index=d_genome)
            weight_sum_pd = weight_pd.groupby(level=0).sum()
            try:
                D3[weight_sum_pd.index, :k_vec.shape[0]] += weight_sum_pd.values
            except IndexError as ie:
                print(ie)
    D3 = D3/D3.sum(axis=0)
    return bin_count, comb_counts, D3


def get_expectation(distrib, bins, loc_num, correct_bias, bias_series, resolution, num_trial):
    dist = [(bins[i]-bins[i-1])//resolution for i in range(1, len(bins))]
    probs = [distrib[item] for item in dist]
    # prob=np.product(probs)/loc_num

    # 实现1
    # prob=np.product(probs)/(loc_num-range_bin+1)
    # range_bin=(bins[-1]-bins[0])//resolution

    # 实现2
    prob = np.product(probs)/(loc_num)

    # print(num_trial)

    if correct_bias:
        bias_list = bias_series.loc[bins]
        total_bias = np.product(bias_list)
    else:
        total_bias = 1
    if total_bias <= 0:
        is_valid = False
    else:
        is_valid = True

    prob = prob*total_bias
    expectation = prob*num_trial
    # p_val=bdtrc(k-1,num_trial,prob)
    return prob, expectation, total_bias, is_valid


def get_test_num(bin_count, dist_limit, k_mer, baseline):
    chr_len_bin = bin_count.shape[0]
    test_num = 0
    for i in range(chr_len_bin):
        index_upper_bound = int(min(i+dist_limit, chr_len_bin))
        bin_count_window = bin_count[i:index_upper_bound]
        n_window = np.sum(bin_count_window > baseline)
        test_num += comb(n_window-1, k_mer-1)
    return test_num


def get_qvalues(pvalues, num_test):
    # Assuming it has been sorted
    pvalues = np.array(pvalues)
    qvalues = pvalues*num_test/np.arange(1, pvalues.shape[0]+1)
    # Guaranteed to increase monotonically
    qvalues = np.maximum.accumulate(qvalues)
    qvalues[qvalues > 1] = 1
    return qvalues


def detect_by_chrom(chrom, chromosome_size, min_sup, output_path, loop_dir, external_loop_path, resolution, max_hyper_loop, loop_lowerbound, loop_upperbound, loop_fdr_cut, start_step, filter_by_loop, correct_bias):
    print('{} start'.format(chrom))

    # # Mining frequent interaction patterns
    time_start = time.time()
    if start_step <= 2:
        print('Step 2, mine frequent pattern................')

        if filter_by_loop:
            if external_loop_path == 'void':
                loop_file = '{}/{}_fithic/FitHiC.spline_pass2.res{}.significances.txt.gz'.format(
                    loop_dir, chrom, resolution)
                loop_df = pd.read_csv(loop_file, sep='\t')
            else:
                loop_df = pd.read_csv(external_loop_path, sep='\t')
                loop_df = loop_df[loop_df.chr1 == chrom]
            loop_df = loop_df.rename(
                {'p-value': 'p_value', 'q-value': 'q_value'}, axis=1)
            loop_significant = loop_df.query("q_value <= 0.1")
            node_set = set(loop_significant.fragmentMid1.values) | set(
                loop_significant.fragmentMid2.values)
            node_series = pd.Series(list(node_set))
            node_file = output_path+'/{}_nodes.txt'.format(chrom)
            node_series.to_csv(node_file, index=False, header=False)
        detect_frequent_pattern(min_sup, max_hyper_loop,
                                output_path, chrom, filter_by_loop)
        print('Step 2 done.')
    
    time_end2=time.time()
    with open(output_path+'/times.txt',mode='a') as tf:
        tf.write('\n{}'.format(time_end2-time_start))

    # # Testing connectivity
    if start_step <= 3 and filter_by_loop:
        print('Step 3, test connectivity by significant loops................')
        detect_connected_pattern(
            loop_dir, external_loop_path, chrom, resolution, loop_fdr_cut, output_path)
        print('Step 3 done.')

    time_end3=time.time()
    with open(output_path+'/times.txt',mode='a') as tf:
        tf.write('\n{}'.format(time_end3-time_end2))

    # # Testing significance
    if start_step <= 4:
        print('Step 4, test the significance of hyperloops................')
        detect_significant_pattern(output_path, chrom, chromosome_size, loop_dir,
                                   external_loop_path, resolution, filter_by_loop, max_hyper_loop, loop_lowerbound, loop_upperbound, correct_bias)
        print('Step 4 done.')
    
    time_end4=time.time()
    with open(output_path+'/times.txt',mode='a') as tf:
        tf.write('\n{}'.format(time_end4-time_end3))
    print('{} end'.format(chrom))
    return '{} end'.format(chrom)


def join_all_result_file(output_path, chroms, filter_by_loop):
    result_list = []
    for chrom in chroms:
        if filter_by_loop:
            hyperloop_out_file = '{}/{}_connected_significant_hyperloop.csv'.format(
                output_path, chrom)
        else:
            hyperloop_out_file = '{}/{}_significant_hyperloop.csv'.format(
                output_path, chrom)
        hyperloop_df = pd.read_csv(hyperloop_out_file, dtype={
                                   'num_trial': np.float64})
        result_list.append(hyperloop_df)

    hyperloop_all = pd.concat(
        result_list, axis=0, ignore_index=True, sort=False)
    if filter_by_loop:
        hyperloop_all_out_file = '{}/all_connected_significant_hyperloop.csv'.format(
            output_path)
    else:
        hyperloop_all_out_file = '{}/all_significant_hyperloop.csv'.format(
            output_path)
    hyperloop_all.to_csv(hyperloop_all_out_file, index=False)


if __name__ == "__main__":
    main()
