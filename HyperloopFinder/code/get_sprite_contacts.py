#!/usr/bin/env python
# coding: utf-8

# In[1]:


from enum import Enum
from itertools import combinations
from random import random
import numpy as np
import subprocess
import scipy as sp
import pandas as pd 
from scipy import sparse
import argparse


# # Contacts

# In[2]:


class Downweighting(Enum):
    """An enumeration of downweighting schemes.

    NONE -- No downweighting. Each contact has a value of 1.
    N_MINUS_ONE -- A contact from a cluster of n reads has a value of
                   1 /(n - 1).
    TWO_OVER_N -- A contact form a cluster of n reads has a value of
                  2 / n.
    UNKNOWN -- A default downweighing scheme for error checking.
    """

    NONE = 1
    N_MINUS_ONE = 2
    TWO_OVER_N = 3
    UNKNOWN = 4


class Contacts:
    """A class for making heatmaps from chromosomal conformation data.

    This class primarily deals with the cluster files from the Guttman Lab
    SPRITE workflow, but also contains some auxiliary methods for other Hi-C
    file formats.
    """


    def __init__(self, chromosome, build = "mm9", resolution = 1000000,
                 downweighting = "none"):
        """Constructs an instance of the Contacts class.

        Args:
            chromosome (str): The chromosome to visualize, or "genome" for an
                interchromosomal heatmap.
            build (str): The genome assembly.
            resolution (int): The resolution of the heatmap, in bp.
            downweighting (str): The downweighting scheme.
        """

        self._chromosome = chromosome
        self._resolution = resolution
        self._assembly = Build(build, resolution)

        if downweighting == "none":
            self._downweighting = Downweighting.NONE
        elif downweighting == "n_minus_one":
            self._downweighting = Downweighting.N_MINUS_ONE
        elif downweighting == "two_over_n":
            self._downweighting = Downweighting.TWO_OVER_N
        else:
            self._downweighting = Downweighting.UNKNOWN

        if self._chromosome == "genome":
            self.init_genome_matrix()
        else:
            self.init_chromosome_matrix()


    def get_raw_contacts_from_sprite_file(self, clusters_file,
                min_cluster_size = 2, max_cluster_size = 1000):
        """Parses a SPRITE clusters file and stores the relevant contacts."""

        if self._chromosome.startswith("chr"):
            self.get_raw_intrachromosomal_contacts_from_sprite_file(clusters_file,
                    min_cluster_size, max_cluster_size)
        elif self._chromosome == "genome":
            self.get_raw_interchromosomal_contacts_from_sprite_file(clusters_file,
                    min_cluster_size, max_cluster_size)
        else:
            raise Exception("Chromosome ID must start with 'chr' or " + 
                              "equal 'genome'")


    def get_raw_interchromosomal_contacts_from_sprite_file(self, clusters_file,
            min_cluster_size, max_cluster_size):
        """Parses a SPRITE clusters file and stores all contacts.

        This method is used for generating a genome-wide interchromosomal
        heatmap.

        Args:
            clusters_file (str): The path to the SPRITE clusters file
            min_cluster_size (int): Ignore clusters with a smaller size
            max_cluster_size (int): Ignore clusters with a larger size
        """

        with open(clusters_file, 'r') as f:

            for line in f:
                reads = line.split()[1:]
                if not min_cluster_size <= len(reads) <= max_cluster_size:
                    continue
                bins = set()

                for read in reads:
                    _, coord = read.split('_')
                    chrom, start, end = coord.replace('-', ':').split(':')
                    genome_pos = self.get_genomic_position(chrom, start)
                    if genome_pos is not None:  # genome_pos == None if chrom not in dict
                        bins.add(genome_pos)

                self.add_bins_to_contacts(bins)


    def get_raw_intrachromosomal_contacts_from_sprite_file(self, clusters_file,
                min_cluster_size, max_cluster_size):
        """Parses a SPRITE clusters file and stores contacts from within a
        single chromosome.

        Note:
            This method does not take a chromosome as an argument. The
            chromosome is specified when the Contacts object is constructed. If
            you need a heatmap for each of multiple chromosomes, you'll need to
            make a new Contacts object for each.

        Args:
            clusters_file (str): The path to the SPRITE clusters file
            min_cluster_size (int): Ignore clusters with a smaller size
            max_cluster_size (int): Ignore clusters with a larger size
        """

        with open(clusters_file, 'r') as f:

            for line in f:
                reads = line.split()[1:]
                if not min_cluster_size <= len(reads) <= max_cluster_size:
                    continue
                bins = set()

                for read in reads:
                    chrom, start = read.split(':')
                    if chrom == self._chromosome:
                        read_bin = int(start) // self._resolution
                        bins.add(read_bin)

                self.add_bins_to_contacts(bins)                


    def get_raw_intrachromosomal_contacts_from_aiden_hic_file(self, hic_file):
        """Parses a Hi-C file from the Aiden lab and stores the contacts.

        Erez's Hi-C files have three columns. Columns one and two contain
        positions on a chromosome that interact. Column three is the number
        of interactions. As an example:

        1213142 1213143 10
        1213142 1213144 4
        1213142 1213145 1

        Note:
            This method only handles intrachromosomal Hi-C files, e.g.,
            chr2-against-chr2.
        """

        with open(hic_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                pos1, pos2, count = line.split()
                pos1 = int(pos1) // self._resolution
                pos2 = int(pos2) // self._resolution
                count = int(float(count))
                self._contacts[pos1][pos2] = count
                self._contacts[pos2][pos1] = count


    def get_raw_contacts_from_ren_hic_file(self, hic_file):
        """Parses a Hi-C file from the Ren lab and stores the contacts.

        The Ren lab's Hi-C files have seven columns. Columns two and three
        contain the chromosome and position of a read in a contact. Columns
        five and six contain the chromosome and position of the other read.
        Each line corresponds to a single contact. As an example:

        HWI-ST216_0305:5:1104:16545:105833#AGTAAG chr1  3000000 - chr1  5404761 +
        HWI-ST216_0305:4:2208:8611:50989#AAATGA chr1  3000001 + chr1  3000252 -
        HWI-ST216_0305:5:2307:15998:173700#GGTTGT chr1  3000001 + chr1  3000218 -
        """

        if self._chromosome.startswith("chr"):
            self.get_raw_intrachromosomal_contacts_from_ren_hic_file(hic_file)
        elif self._chromosome == "genome":
            self.get_raw_interchromosomal_contacts_from_ren_hic_file(hic_file)
        else:
            raise Exception("Chromosome ID must start with 'chr' or " + 
                              "equal 'genome'")


    def get_raw_intrachromosomal_contacts_from_ren_hic_file(self, hic_file):
        """Parses a Hi-C file from the Ren lab and stores the intrachromosomal
        contacts on one chromosome.
        """

        assert self._chromosome.startswith("chr")

        with open(hic_file, 'r') as f:
            for line in f:
                _, chrom1, pos1, _, chrom2, pos2, _ = line.split()
                if self._chromosome == chrom1 == chrom2:
                    pos1 = int(pos1) // self._resolution
                    pos2 = int(pos2) // self._resolution
                    self._contacts[pos1][pos2] += 1
                    self._contacts[pos2][pos1] += 1


    def get_raw_interchromosomal_contacts_from_ren_hic_file(self, hic_file):
        """Parses a Hi-C file from the Ren lab and stores all contacts for a
        genome-wide interchromosomal heatmap.
        """

        assert self._chromosome == "genome"

        with open(hic_file, 'r') as f:

            for line in f:
                _, chrom1, pos1, _, chrom2, pos2, _ = line.split()
                bin1 = self.get_genomic_position(chrom1, pos1)
                bin2 = self.get_genomic_position(chrom2, pos2)

                # Bins == None if get_genomic_position passed unknown chrom
                if bin1 is not None and bin2 is not None:
                    self._contacts[bin1][bin2] += 1
                    self._contacts[bin2][bin1] += 1


    def get_genomic_position(self, chromosome, position):
        """Converts a chromosome and position to the appropriate heatmap index
        for a genome-wide interchromosomal heatmap

        For example, passing "chr3" and 10000 will return a value which
        (conceptually) is (10000 + length(chr2) + length(chr1)) / resolution.

        Args:
            chromosome (str): The chromosome, e.g., "chr1".
            position (int): The position along the chromosome.
        """

        read_bin = int(position) // self._resolution
        offset = self._assembly.get_offset(chromosome)
        if offset is not None:
            return read_bin + offset


    def add_bins_to_contacts(self, bins):
        """Stores all pairwise contacts implied by one SPRITE cluster."""

        if len(bins) > 1:
            if self._downweighting == Downweighting.TWO_OVER_N:
                inc = 2.0 / len(bins)
            elif self._downweighting == Downweighting.N_MINUS_ONE:
                inc = 1.0 / (len(bins) - 1)
            else:
                assert self._downweighting == Downweighting.NONE
                inc = 1.0

            for bin1, bin2 in combinations(bins, 2):
                try:
                    self._contacts[bin1][bin2] += inc
                    self._contacts[bin2][bin1] += inc
                except IndexError as ie:
                    print(ie)


    def zero_diagonal_entries(self):
        """Sets all diagonal entries in the internal heatmap matrix to zero."""

        for i in range(len(self._contacts)):
            self._contacts[i][i] = 0


    def init_chromosome_matrix(self):
        """Initializes an internal heatmap matrix with a number of rows
        determined by this object's assembly, chromosome and resolution (e.g.,
        chr1 on mm9 at a 100 Mb resolution.

        This method is used for single intrachromosomal heatmaps only.
        """

        chromosome_size = self._assembly.get_size(self._chromosome)
        num_bins = -(-chromosome_size // self._resolution)
        self._contacts = np.zeros((num_bins, num_bins))


    def init_genome_matrix(self):
        """Initializes an internal heatmap matrix with a number of rows
        determined by this object's assembly and resolution (e.g., mm9 at
        100 Mb resolution.

        This method is used for genome-wide interchromosomal heatmaps only.
        """

        num_bins = 0
        for chromosome_size in self._assembly._chromsizes.values():
            num_bins += -(-chromosome_size // self._resolution)
        self._contacts = np.zeros((num_bins, num_bins))


    def write_contacts_to_file(self, outfile, fmt):
        """Writes the internal heatmap matrix to file.

        Args:
            outfile (str): The path to write to.
            fmt (str): The numerical format to write.
        """

        np.savetxt(outfile, self._contacts, delimiter = "\t", fmt = fmt)


    def ice_raw_contacts(self, raw_contacts_file, bias_file, iterations,
                hicorrector_path):
        """Calls Hi-Corrector to apply IC normalization to the internal heatmap
        matrix.

        This method generates a file of biases by calling the ic executable,
        then scales each cell of the internal heatmap matrix by the two
        appropriate factors in that file.
 
        Args:
            raw_contacts_file (str): The file containing a raw contacts heatmap.
            bias_file (str): The path to write the Hi-Corrector output to.
            iterations (int): The number of Hi-Corrector iterations to perform.
            hicorrector_path (str): The path to the Hi-Corrector ic program.
        """

        biases = self.calculate_bias_factors(raw_contacts_file = raw_contacts_file,
                bias_file = bias_file, hicorrector = hicorrector_path,
                iterations = iterations)

        median_diagonal_value = self.get_median_diagonal_value()

        for row in range(self._contacts.shape[0]):
            for col in range(self._contacts.shape[1]):
                val = self._contacts[row][col]
                if val > 0:
                    val /= (biases[row] * biases[col])
                    self._contacts[row][col] = val


    def truncate_to_median_diagonal_value(self):
        """Scales all values in the internal heatmap matrix relative to the
        median value of the +1 and the -1 diagonals (offset from main diagonal
        by +/- 1).

        If the median value is 10, a value of 3 will be scaled to 0.3 and a
        value of 9 will be scaled to 0.9. A value of 11 would be set to 1
        (since 11 > 10) rather than being set to 1.1.
        """

        median_diagonal_value = self.get_median_diagonal_value()

        for row in range(self._contacts.shape[0]):
            for col in range(self._contacts.shape[1]):
                val = self._contacts[row][col]
                val = (1 if val >= median_diagonal_value
                       else val / median_diagonal_value)
                self._contacts[row][col] = val

 
    def calculate_bias_factors(self, raw_contacts_file, bias_file, hicorrector,
                iterations):
        """Runs Hi-Corrector on the raw contacts heatmap.

        The bias file that Hi-Corrector outputs is subsequently read and
        returned as a list of floats.

        Note:
            Hi-Corrector cannot access the internal numpy matrix of this
            object. The matrix needs to be written to disk first, then read by
            Hi-Corrector.

        Args:
            raw_contacts_file (str): The file containing a raw contacts heatmap.
            bias_file (str): The path to write the Hi-Corrector output to.
            iterations (int): The number of Hi-Corrector iterations to perform.
            hicorrector_path (str): The path to the Hi-Corrector ic program.
        """

        skip_first_row = "0"    # 0 == don't skip
        skip_first_column = "0"
        num_lines = self._contacts.shape[0]
        subprocess.check_call([hicorrector, raw_contacts_file, str(num_lines),
                str(iterations), skip_first_row, skip_first_column, bias_file])
        return self.parse_bias_file(bias_file)


    def parse_bias_file(self, bias_file):
        """Parses a Hi-Corrector ic output file and returns the values as a
        list of floats.
        """

        biases = []

        with open(bias_file) as f:
            for line in f:
                biases.append(float(line.strip()))
        return biases


    def get_median_diagonal_value(self):
        """Returns the median diagonal value of this object's internal heatmap
        matrix."

        Note:
            The median diagonal value is actually the median of the two
            diagonals offset by +1 and -1.
        """

        diagonal_values = []

        for i in range(self._contacts.shape[0] - 1):
            diagonal_values.append(self._contacts[i + 1][i])
            diagonal_values.append(self._contacts[i][i + 1])

        return np.median(diagonal_values)


    def downsample(self, target):
        """Downsample the internal heatmap matrix to a target number of
        contacts.

        Args:
            target (int): The number of contacts to downsample to.
        """

        dim = len(self._contacts)

        total_contacts = 0

        # Only sum contacts from diagonal and upper-triangle
        # Otherwise, will double-count.
        for i in range(dim):
            for j in range(i + 1):
                total_contacts += self._contacts[i][j]

        downsample_ratio = float(target) / total_contacts

        for i in range(dim):
            for j in range(i + 1):
                num_contacts = self._contacts[i][j]
                for contact in range(num_contacts):
                    if random() < downsample_ratio:
                        num_contacts -= 1
                self._contacts[i][j] = num_contacts
                self._contacts[j][i] = num_contacts


# # Assembly

# In[3]:


from collections import OrderedDict

__author__ = "Noah Ollikainen, Charlotte A Lai, Peter Chovanec"

class Assembly(object):
    _chromsizes = None
    _resolution = None

    def __init__(self, resolution):
        self._resolution = resolution
        self.init_offsets()

    def init_offsets(self):
        count = 0
        self._offsets = OrderedDict()
        for (chrom, size) in self._chromsizes.items():
            self._offsets[chrom] = count
            count += -(-size // self._resolution)


    def get_size(self, chrom):
        return self._chromsizes.get(chrom)


    def get_offset(self, chrom):
        return self._offsets.get(chrom)


    def get_position(self, n):
        for (chrom, offset) in reversed(self._offsets.items()):
            if n - offset >=  0:
                return (chrom, (n - offset) * self._resolution)
        raise ValueError


    def get_index(self, chrom, pos):
        read_bin = pos // self._resolution
        offset = self.get_offset(chrom)
        if offset is not None:
            return read_bin + offset


class Mm9(Assembly):

    def __init__(self, resolution):
        self._chromsizes = OrderedDict([
            ('chr1', 197195432),
            ('chr2', 181748087),
            ('chr3', 159599783),
            ('chr4', 155630120),
            ('chr5', 152537259),
            ('chr6', 149517037),
            ('chr7', 152524553),
            ('chr8', 131738871),
            ('chr9', 124076172),
            ('chr10', 129993255),
            ('chr11', 121843856),
            ('chr12', 121257530),
            ('chr13', 120284312),
            ('chr14', 125194864),
            ('chr15', 103494974),
            ('chr16', 98319150),
            ('chr17', 95272651),
            ('chr18', 90772031),
            ('chr19', 61342430),
            ('chrX', 166650296),
            ('chrY', 15902555),
            ('chrM', 16299)])
        super(Mm9, self).__init__(resolution)

    @classmethod
    def is_named(cls, name):
        return name == "mm9"


class Mm10(Assembly):

    def __init__(self, resolution):
        self._chromsizes = OrderedDict([
            ('chr1', 195471971),
            ('chr2', 182113224),
            ('chr3', 160039680),
            ('chr4', 156508116),
            ('chr5', 151834684),
            ('chr6', 149736546),
            ('chr7', 145441459),
            ('chr8', 129401213),
            ('chr9', 124595110),
            ('chr10', 130694993),
            ('chr11', 122082543),
            ('chr12', 120129022),
            ('chr13', 120421639),
            ('chr14', 124902244),
            ('chr15', 104043685),
            ('chr16', 98207768),
            ('chr17', 94987271),
            ('chr18', 90702639),
            ('chr19', 61431566),
            ('chrX', 171031299),
            ('chrY', 91744698),
            ('chrM', 16299)])
        super(Mm10, self).__init__(resolution)

    @classmethod
    def is_named(cls, name):
        return name == "mm10"


class Hg19(Assembly):

    def __init__(self, resolution):
        self._chromsizes = OrderedDict([
            ('chr1', 249250621),
            ('chr2', 243199373),
            ('chr3', 198022430),
            ('chr4', 191154276),
            ('chr5', 180915260),
            ('chr6', 171115067),
            ('chr7', 159138663),
            ('chr8', 146364022),
            ('chr9', 141213431),
            ('chr10', 135534747),
            ('chr11', 135006516),
            ('chr12', 133851895),
            ('chr13', 115169878),
            ('chr14', 107349540),
            ('chr15', 102531392),
            ('chr16', 90354753),
            ('chr17', 81195210),
            ('chr18', 78077248),
            ('chr19', 59128983),
            ('chr20', 63025520),
            ('chr21', 48129895),
            ('chr22', 51304566),
            ('chrX', 155270560),
            ('chrY', 59373566),
            ('chrM', 16571)])
        super(Hg19, self).__init__(resolution)

    @classmethod
    def is_named(cls, name):
        return name == "hg19"


class Hg38(Assembly):

    def __init__(self, resolution):
        self._chromsizes = OrderedDict([
            ('chr1', 248956422),
            ('chr2', 242193529),
            ('chr3', 198295559),
            ('chr4', 190214555),
            ('chr5', 181538259),
            ('chr6', 170805979),
            ('chr7', 159345973),
            ('chr8', 145138636),
            ('chr9', 138394717),
            ('chr10', 133797422),
            ('chr11', 135086622),
            ('chr12', 133275309),
            ('chr13', 114364328),
            ('chr14', 107043718),
            ('chr15', 101991189),
            ('chr16', 90338345),
            ('chr17', 83257441),
            ('chr18', 80373285),
            ('chr19', 58617616),
            ('chr20', 64444167),
            ('chr21', 46709983),
            ('chr22', 50818468),
            ('chrX', 156040895),
            ('chrY', 57227415),
            ('chrM', 16569)])
        super(Hg38, self).__init__(resolution)

    @classmethod
    def is_named(cls, name):
        return name == "hg38"


def Build(name, resolution):
    for cls in Assembly.__subclasses__():
        if cls.is_named(name):
            return cls(resolution)
    raise ValueError


# # main

# In[17]:

def main():
    args = parse_arguments()

    contacts = Contacts(build = args.assembly,
                        chromosome = args.chromosome,
                        resolution = args.resolution,
                        downweighting = args.downweighting)

    contacts.get_raw_contacts_from_sprite_file(
            clusters_file = args.clusters,
            min_cluster_size = args.min_cluster_size,
            max_cluster_size = args.max_cluster_size)

    #contacts.write_contacts_to_file(outfile = args.raw_contacts, fmt = "%1f")

    data=contacts._contacts
    data=np.triu(data,1)

    data_coo=sp.sparse.coo_matrix(data)

    data_df=pd.DataFrame(data={'chr1':args.chromosome,'cord1':data_coo.row*args.resolution+args.resolution//2,'chr2':args.chromosome,'cord2':data_coo.col*args.resolution+args.resolution//2,'frequency':data_coo.data})

    data_df.to_csv(args.raw_contacts,sep='\t',header=None,index=False,compression='gzip')

    '''
    import seaborn as sns
    import matplotlib.pyplot as plt

    step=100
    for i in range(data.shape[0]//step):
        start=step*i
        plt.figure(figsize=(14,14))
        end=start+step
        sns.heatmap(data[start:end,start:end],square=True)
        plt.show()
    '''

    mid_list=[i*args.resolution+args.resolution//2 for i in range(data.shape[0])]

    row_sum=data.sum(axis=1).astype(int)
    fragment_df=pd.DataFrame({'chr':args.chromosome,'extra':1,'mid':mid_list,'frequency':row_sum,'map':1})

    fragment_df.to_csv(args.fragments,sep='\t',header=None,index=False,compression='gzip')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates an ICEd matrix of intra- or ' +
                'interchromosomal contacts from a SPRITE clusters file.')

    parser.add_argument('--clusters',
                        metavar = "FILE",
                        action = "store",
                        help = "The input clusters file.")

    parser.add_argument('--raw_contacts',
                        metavar = "FILE",
                        action = "store",
                        help = "An intermediate file of raw (un-ICEd) " + 
                               "contacts.")

    parser.add_argument('--biases',
                        metavar = "FILE",
                        action = "store",
                        help = "An intermediate file of biases generated by " +
                               "hicorrector.")

    parser.add_argument('--iced',
                        metavar = "FILE",
                        action = "store",
                        help = "An intermediate file of ICEd contacts.")

    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        help = "The output ICEd and scaled contacts file.")

    parser.add_argument('--assembly',
                        metavar = "ASSEMBLY",
                        action = "store",
                        choices = ["mm9", "mm10", "hg19", "hg38"],
                        default = "mm10",
                        help = "The genome assembly. (default mm10)")

    parser.add_argument('--chromosome',
                        metavar = "CHROM",
                        action = "store",
                        help = "The chromosome of interest. Input 'genome' " +
                               "for an interchromosomal matrix.")

    parser.add_argument('--max_cluster_size',
                        metavar = 'MAX',
                        type = int,
                        action = 'store',
                        default = 1000,
                        help = "Skip read-clusters with more reads than MAX. (default 1000)")

    parser.add_argument('--min_cluster_size',
                        metavar = 'MIN',
                        type = int,
                        action = 'store',
                        default = 2,
                        help = "Skip read-clusters with fewer reads than MIN. (default 2)")

    parser.add_argument('--resolution',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        default = 1000000,
                        help = "Contact matrix resolution in bp. (default 1000000 = 1 Mb)")

    parser.add_argument('--downweighting',
                        metavar = 'DW',
                        action = 'store',
                        default = "none",
                        choices = ["none", "n_minus_one", "two_over_n"],
                        help = "Downweighting strategy")

    parser.add_argument('--hicorrector',
                        metavar = "FILE",
                        action = 'store',
                        default = "/groups/guttman/software/hicorrector/1.2/bin/ic",
                        help = "Path to hicorrector ic. (default /mnt/new-storage/Software/hicorrector/1.2/bin/ic)")

    parser.add_argument('--iterations',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        default = 100,
                        help = "Number of ICE iterations (default 100)")
    parser.add_argument('--fragments',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'The output fragments file')

    return parser.parse_args()


if __name__ == "__main__":
    main()
