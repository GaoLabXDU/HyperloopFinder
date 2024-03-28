# HyperloopFinder
An efficient and effective Python tool for detecting regulatory multi-way chromatin contacts

## Overview
Recent advances in chromatin conformation capture technologies, such as SPRITE and Pore-C, have enabled the detection of simultaneous contacts among multiple chromatin loci. This has made it possible to investigate the cooperative transcriptional regulation involving multiple genes and regulatory elements at the resolution of a single molecule. However, these technologies are unavoidably subject to the random polymer looping effect and technical biases, making it challenging to distinguish genuine regulatory relationships directly from random polymer interactions.

Therefore, we developed HyperloopFinder, a efficient and effective tool for identifying regulatory multi-way chromatin contacts (hyperloops) by jointly modeling the random polymer looping effect and technical biases to estimate the statistical significance of multi-way contacts.

## Required Package

HyperloopFinder could be installed in a linux-like system. The HyperloopFinder requires the following dependencies. We recommend using [Anaconda python distribution](https://www.anaconda.com/) to install the below python packages. Generally, the first six dependencies are directly supported in the default Anaconda distribution.

1. Python (tested on 3.12.+)
2. numpy (tested on 1.26.+)
3. pandas (tested on 2.2.+)
4. matplotlib (tested on 3.8.+)
5. networkx (tested on 3.2.+)
6. scipy (tested on 1.12.+)
7. fithic (tested on 2.0.+) （For detailed installation steps, please refer to: [https://github.com/ay-lab/fithic](https://github.com/ay-lab/fithic)）

## Installation

Download HyperloopFinder by

```shell
git clone https://github.com/gaolabXDU/HyperloopFinder
```
HyperloopFinder also require fpgrowth shell tool, please download it according to your OS from [https://borgelt.net/fpgrowth.html](https://borgelt.net/fpgrowth.html) and put it in `HyperloopFinder/utils` directory. If there is no version supported by your system, you can download the source code from this URL to compile it.

HyperloopFinder also require GNU parallel shell tools ([https://www.gnu.org/software/parallel/](https://www.gnu.org/software/parallel/)) for executing jobs in parallel using one or more CPUs. you can install it by `yum install parallel` command on CentOS/RHEL operating systems or `sudo apt install parallel` command on Ubuntu / Debian operating systems.

## Testing installation
```
cd ./code
bash hyperloopfinder_test.sh
```

## Usage
```shell
usage: hyperloopfinder.py [-h] [--cluster_file FILE] [--loop_dir FILE]
                          [--chromosome_size_file FILE]
                          [--external_loop_path FILE] [--filter_by_loop INT]
                          [--correct_bias INT] [-o FILE]
                          [--min_support_count INT] [--max_cluster_size MAX]
                          [--max_hyperloop_size MAX] [--loop_lowerbound MAX]
                          [--loop_upperbound INT] [--loop_fdr_cut FLOAT]
                          [--min_cluster_size INT] [--resolution INT]
                          [--num_thread INT] [--chromosome List [List ...]]
                          [--start_step INT]

Detect hyperloops from multi-way chromatin interaction data. This program contains the core functions of hyperloopFinder, which needs to input the cluster data file of multi-way chromatin interaction and the pair-wise chromatin loops detected by other methods such as Fit-HiC2, and finally output significant multi-way chromatin interactions.

optional arguments:
  -h, --help
                        Show this help message and exit.
  --cluster_file FILE
                        File path of cluster file.
  --loop_dir FILE
                        File path of chromatin loops.
  --chromosome_size_file FILE
                        File path of chromosome size.
  --external_loop_path FILE
                        File path of external loops.
  --filter_by_loop INT
                        Whether to filter hyperloops by significant loops? 1
                        is yes, 0 is no.
  --correct_bias INT
                        Whether to correct bias? 1 is yes, 0 is no.
  -o FILE, --output_path FILE
                        Path of outputs.
  --min_support_count INT
                        Minimum support count for FPGrowth algorithm.
  --max_cluster_size MAX
                        Maximum cluster size of multi-way chromatin
                        interaction data.
  --max_hyperloop_size MAX
                        Maximum hyperloop size that user defined.
  --loop_lowerbound MAX
                        Minimum hyperloop distance that user defined.
  --loop_upperbound INT
                        Maximum hyperloop distance that user defined.
  --loop_fdr_cut FLOAT
                        Q value cut for significant pairwise loops.
  --min_cluster_size INT
                        Minimum cluster size of multi-way chromatin
                        interaction data.
  --resolution INT
                        Contact matrix resolution in bp.
  --num_thread INT
                        Number of CPU threads.
  --chromosome List [List ...]
                        Chromosome list to be run by this tool.
  --start_step INT
                        From which step does the program start?

```
To facilitate the use of the software, the user can configure the `code/hyperloopfinder_test.sh` script in the code folder with the following example. This bash script implements all processes for detecting hyperloops, mainly including calling Fit-HiC2 to detect pairwise loops and calling hyperloopFinder to detect mult-way chromatin loops. Users can modify these parameters according to their own needs.

```
task_suffix=time_test
start_step=1
filter_by_loop=0
correct_bias=1
fithic_pass=2

chroms=('chr21')
assembly=hg38
downweighting=none
num_thread=6

min_cluster_size=3
max_cluster_size=100
resolution=25000
min_support_count=15
max_hyperloop_size=3
loop_lowerbound=25000
loop_upperbound=10000000
loop_fdr_cut=0.05

cluster_file=../data/GM12878/random_shuffled_porec_clusters/porec_GM12878_Nlalll_grch38_random_shuffled_repeat0.clusters
result_dir=../result/GM12878/hyperloops_results_${resolution}_${downweighting}_$task_suffix
loop_dir=${result_dir}/loops
hyperloop_dir=${result_dir}/hyperloops
chromosome_size_file=../data/hg38.chrom.sizes.txt
```

The following is a detailed explanation of each parameter.

### task_suffix
Task suffix for naming the folder of output. Users can set a friendly name prefix to distinguish other tasks. This prefix will be combined with other important parameters to set the name of the output directory.

### start_step
The task can be started at a user-specified step. This parameter allows the user to start executing the pipeline in the middle step. When the program runs incorrectly, on the premise of ensuring that the previous steps are executed correctly, the user can directly start to execute the error step after modifying the parameters, instead of starting the execution again.

Accepts:
0: detect loop
1: bin and split cluster file
2: mine frequent pattern
3: test connectivity by significant loops
4: test the significance of hyperloops
5: joint all results from different chromosomes

### filter_by_loop
Determine whether to use pairwise loops for filtering multi-way chromatin contacts. Our method supports the detection of significant multi-way interactions in all candidate sets, but this is usually not a good choice, because it can cause huge time and resource consumption and even cause computer crashes.

Accepts: 0 for not, 1 for yes.

### correct_bias
Determine whether to correct technical biases. It is recommended to correct the technical biases because the technical biases usually have a great impact on the significance estimation of hyperloops.

Accepts: 0 for no, 1 for yes

### fithic_pass
Determine the number of spline passes of Fit-HiC2 software. This parameter is only used by Fit-HiC2 software.

Accepts: integer values greater than or equal to 1.


### chroms
Specify the list of chromosomes to be tested for hyperloopFinder.

Accepts: the names of the chromosomes must have appeared in the input data.

### assembly
Specify the reference genome version of cluster data.

### downweighting
Accepts: none only, this parameter is only passed to get_sprite_contacts.py script.

### num_thread
Specify the number of threads that used by HyperloopFinder. Using one thread for one chromosome usually increases the speed of the software operation dramatically. But be careful, it may exceed the limit of the user's machine memory.

Accepts: integer values greater than or equal to 1.

### min_cluster_size

Specify the minimum cluster size to be used in hyperloop mining.

### max_cluster_size
Specify the maximum cluster size to be used in hyperloop mining. Very large clusters may be caused by technical bias. Sometimes it is necessary to filter out large clusters, especially for SPRTIE data.

### resolution
Specify the resolution of bin size for binning interaction loci of multi-way contacts.

### min_support_count
Specify the minimum support count parameter for the fpgrowth algorithm. The FP-growth algorithm finds all frequent itemsets whose support is greater than the minimum support count. Here the minimum support count refers to the minimum frequency of multi-way interactions. The minimum support count is greater than or equal to 1. Smaller support thresholds usually produce a large number of results, so be careful with setting small support counts.

### max_hyperloop_size
Specify the maximum number of interaction loci included in the hyperloop. This limits the number of sites contained in the maximum hyperloop and can be set as needed.

### loop_lowerbound

Integer representing lower bound for the intrachromosomal interactions to be considered in base pairs.

### loop_upperbound
Integer representing upper bound for the intrachromosomal interactions to be considered in base pairs.

### loop_fdr_cut
False discovery rate threshold for determining whether a pairwise loop detected by Fit-HiC2 is significant.

### cluster_file
Specify the location of the input file for the multi-way clusters.

### result_dir
Specify the location and format of the output folder.

### loop_dir
Specify the output folder for Fit-HiC2.

### hyperloop_dir
Specify the output folder for hyperloops.

### chromosome_size_file
Specify the size of each chromosome.


## Input Files

### Cluster file
This file contains interaction clusters of SPRITE and Pore-C and other multi-way chromatin interaction data. Each row in the data defines a multi-way interaction, containing the interaction ID and any number of chromatin sites, separated by tabs. Each chromatin site is represented by a chromosome number and chromosome coordinate, separated by semicolons. For example:
￼￼￼￼

```
112	chr6:154476803	chr1:151928081	
123	chr1:188725947	chr1:188749586	
157	chr8:13258973	chr8:28751057	chr8:42960052	
174	chr12:54361769	chr1:89477384	chr16:19915193	
188	chr10:27292831	chr10:28340389	
191	chr4:101389687	chr1:70768664	
197	chr3:86702014	chr3:86625869	chr3:86712997	chr12:54188132	
225	chr5:127766055	chr5:127765741	
231	chr1:123400235	chr1:123237507	chr1:122505444	
251	chr5:99594247	chr2:19412612	
258	chr19:57052551	chr19:6954693	chr14:21200296	chr14:47435112	chr9:132498916
273	chr7:2761780	chr12:15005233	
300	chr18:34539428	chr18:14598397	chr18:15021724	chr18:14598397	chr18:34566026
310	chr11:40133505	chr8:4597171	
333	chr8:123489010	chr8:120935941	chr8:126931236	chr5:179207858	chr5:178776180
```
### Chromosome size file
This file gives the length of each chromosome in the reference genome used. Each line of the file contains the name of the chromosome, the length of the chromosome, separated by the tab symbol. This file is usually available for download at [UCSC](https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=1996329094_gAaxsxTotF9BU2SuFgAzOrjvAuOm&redirect=manual&source=genome.ucsc.edu). For example, this file for hg19 reference genome can be download at [http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes). An example of hg19 chromosome size file is given below:
```
chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276
chr5	180915260
chr6	171115067
chr7	159138663
chrX	155270560
chr8	146364022
chr9	141213431
chr10	135534747
chr11	135006516
chr12	133851895
chr13	115169878
chr14	107349540
chr15	102531392
chr16	90354753
chr17	81195210
chr18	78077248
chr20	63025520
chrY	59373566
chr19	59128983
chr22	51304566
chr21	48129895
```

## Output file

### Hyperloop file

An example of hyperloop file as below:

| chrom | hyperloop                                                            | count | size | distance | prob                   | num_trial    | total_bias         | expectation           | fold_enrichment    | p_value                | q_value                |
|-------|----------------------------------------------------------------------|-------|------|----------|------------------------|--------------|--------------------|-----------------------|--------------------|------------------------|------------------------|
| chr1  | "[107765000, 107805000, 107825000, 107925000, 107935000, 108135000]" | 3     | 6    | 370000   | 1.229433219060843e-15  | 3980443175.0 | 0.2726566613257882 | 4.893689065929011e-06 | 613034.4530645991  | 1.9532429735507844e-17 | 1.0                    |
| chr1  | "[28135000, 28145000, 28155000]"                                     | 17    | 3    | 20000    | 7.259880058401419e-09  | 127855656.0  | 0.4122307105374009 | 0.9282167273482318    | 18.314688260969294 | 3.3020078275386484e-16 | 2.280069424993712e-11  |
| chr1  | "[156275000, 156295000, 156335000]"                                  | 13    | 3    | 60000    | 3.046907567919454e-09  | 127855656.0  | 0.4434136928875694 | 0.3895643658677064    | 33.3706086567854   | 5.324712421637251e-16  | 7.350339521076494e-11  |
| chr1  | "[1825000, 1835000, 1845000]"                                        | 14    | 3    | 20000    | 4.2683424031189644e-09 | 127855656.0  | 0.2423651365449814 | 0.5457317179833915    | 25.65363078351636  | 1.4333380536412846e-15 | 4.948671297099217e-11  |
| chr1  | "[154955000, 155025000, 155045000]"                                  | 11    | 3    | 90000    | 1.794342676011362e-09  | 127855656.0  | 0.3493851247502347 | 0.2294168599302281    | 47.94765303363229  | 1.881303729680033e-15  | 3.8948819245938686e-10 |
| chr1  | "[156255000, 156265000, 156275000]"                                  | 16    | 3    | 20000    | 6.816525793369158e-09  | 127855656.0  | 0.3870561563817115 | 0.8715313769521342    | 18.35848992144633  | 2.3344972079176424e-15 | 5.373312223464037e-11  |
| chr1  | "[30795000, 30815000, 30915000]"                                     | 11    | 3    | 120000   | 2.011330443916793e-09  | 127855656.0  | 0.4702191293530641 | 0.2571599733397528    | 42.77493055059193  | 6.438581833264574e-15  | 1.776990638744523e-09  |
| chr1  | "[9365000, 9425000, 9495000, 9575000, 9735000]"                      | 3     | 5    | 370000   | 3.510481768351902e-14  | 1110644975.0 | 0.244640028896073  | 3.898898935849154e-05 | 76944.80029774406  | 9.877839895108843e-15  | 1.0                    |
| chr1  | "[66275000, 66325000, 66335000, 68305000]"                           | 4     | 4    | 2030000  | 2.269549066766562e-12  | 349928876.0  | 0.4234131464664166 | 0.0007941807539604    | 5036.636785835644  | 1.6564955014060523e-14 | 1.0                    |
| chr1  | "[155865000, 155935000, 155995000, 156015000, 156235000]"            | 3     | 5    | 370000   | 4.597526906585736e-14  | 1110644975.0 | 0.2053268593849361 | 5.106220156226742e-05 | 58751.87336647976  | 2.2188642006410583e-14 | 1.0                    |
| chr1  | "[203025000, 203035000, 203355000, 203425000, 203565000]"            | 3     | 5    | 540000   | 5.434408821660472e-14  | 1110644975.0 | 0.3637038099858185 | 6.035698849872874e-05 | 49704.2691264026   | 3.66447511446979e-14   | 1.0                    |
| chr1  | "[109785000, 109955000, 109965000]"                                  | 10    | 3    | 180000   | 1.643588785801699e-09  | 127855656.0  | 0.3756925613286495 | 0.2101421224029197    | 47.58684211262663  | 3.8235582838716993e-14 | 1.5822266534489482e-08 |
| chr1  | "[42355000, 42525000, 42535000, 42555000, 42815000]"                 | 3     | 5    | 460000   | 5.81160441681998e-14   | 1110644975.0 | 0.21624132608695   | 6.454629242228914e-05 | 46478.26989616586  | 4.4816880910127164e-14 | 1.0                    |

The file uses the csv (Delimiter-separated values) file format. Each row of the file represents a hyperloop, and each column of the file represents an attribute of the hyperloop. The meanings of all attribute names are as follows:

| Column name           |Explanation|
|-----------------|----------------------------------------------------------------------|
| chrom           | Chromosome number.|
| hyperloop       |List of interaction sites contained in hyperloop.|
| count           |The number of times hyperloop appears in the reads of the cluster file.|
| size            |The number of interaction sites contained in the hyperloop.|
| distance        |The distance from the first interaction site to the last interaction site in the hyperloop.|
| prob            |Parameters of the binomial distribution model: interaction probabilities.|
| num_trial       |Parameters of the binomial distribution model: number of experiments.|
| total_bias      |The product of the biases of all interaction sites.|
| expectation     |The expected number of interactions of the hyperloop calculated through the binomial distribution model.|
| fold_enrichment |The observed frequency of multi-way chromatin interactions divided by the expected frequency of interactions.|
| p_value         |P value calculated from the binomial distribution model.|
| q_value         |The P value is corrected for multiple hypothesis testing to obtain the Q value.|