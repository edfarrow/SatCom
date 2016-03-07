__author__ = 'ewjf201'
import itertools
import regex
import json
import os
import multiprocessing as mp
import cPickle
import time
import string


# Generates regex for the searches

class MicrosatelliteRegex:
    def __init__(self, min_seq_length=[18, 9, 6, 4], bases=("A", "C", "T", "G"), motif_lengths=[18],
                 fuzzy_substitution=1):
        """
        :param min_seq_length:
        :param bases:
        :param motif_lengths:
        :param fuzzy_substitution:
        :return:
        """
        self.regex_template = ""
        self.fuzzy_substitution = fuzzy_substitution
        self.min_seq_length = min_seq_length
        self.bases = bases
        self.motif_lengths = motif_lengths
        self.motifs = []
        self.regex = []
        self.regex_template = []
        self.matches = []
        self.output = []

    def _generate_motifs(self):
        # Generates all motifs for the parameters given to the class
        for motif_len in self.motif_lengths:
            # Generates reverse combinations too, so no need to worry about them
            iter_motifs = [''.join(motif) for motif in itertools.product(self.bases, repeat=motif_len)]
            self.motifs.append(iter_motifs)

    def _generate_regex(self):

        def _internal_regex(min_seq_len, motif, subs):
            if len(motif) is 1:
                motif_len = min_seq_len[0]
            elif len(motif) is 2:
                motif_len = min_seq_len[1]
            elif len(motif) is 3:
                motif_len = min_seq_len[2]
            elif len(motif) is 4:
                motif_len = min_seq_len[3]
            elif len(motif) is 5:
                motif_len = min_seq_len[4]
            else:
                motif_len = min_seq_len[0]

            base = "[%s]" % "".join(self.bases)
            if self.fuzzy_substitution == 0:
                motif_1 = r"(?P<whole>(?P<start>%s{,15})(?P<msat>(%s){%d,})(?P<end>%s{,15}))" % (
                    base, motif, motif_len, base)
            else:
                motif_1 = r"(?P<whole>(?P<start>%s{,15})(?P<msat>(%s){%d,}){s<=%d}(?P<end>%s{,15}))" % (
                    base, motif, motif_len, subs, base)
            return motif_1

        regex_templates = [_internal_regex(self.min_seq_length, motifs, self.fuzzy_substitution)
                           for motifs in self.motifs]

        regex_template_full = "|".join(regex_templates)
        self.regex = regex.compile(regex_template_full, flags=(regex.V1 | regex.BESTMATCH))

    def save_regex(self, output_1):
        """
        Outputs a pickle
        :param output_path: The full path to output a pickled file to.
        :param output_name: The name of the file, EXCLUDING any format (eg. "output" rather than "output.pkl").
        :return: True or False
        """
        with open(output_1, "wb") as out_file:
            cPickle.dump(self.regex, out_file, -1)

    def load_regex(self, file_path):
        """
        Loads a pkl-encoded file containing pre-compiled regex
        :param file_path:
        :param file_name:
        :return:
        """
        with open(file_path, 'rb') as in_file:
            regex_raw = cPickle.load(in_file)
            return regex_raw

    def create_regex(self, preloaded=False):
        if preloaded is False:
            self._generate_motifs()
        else:
            self.motifs = ['CCGCG', 'CCCGG',
                           'CCCCG', 'ATGCC', 'ATCGC', 'ATCCG', 'ATCCC', 'ATATC', 'AGGGG', 'AGGGC', 'AGGCG', 'AGGCC',
                           'AGGAT', 'AGCTC', 'AGCGG', 'AGCGC', 'AGCCT', 'AGCCG', 'AGCCC', 'AGCAT', 'AGATG', 'AGATC',
                           'AGAGG', 'AGAGC', 'ACTGG', 'ACTGC', 'ACTCT', 'ACTCG', 'ACTCC', 'ACTAT', 'ACTAG', 'ACGTC',
                           'ACGGG', 'ACGGC', 'ACGCT', 'ACGCG', 'ACGCC', 'ACGAT', 'ACGAG', 'ACCTG', 'ACCTC', 'ACCGT',
                           'ACCGG', 'ACCGC', 'ACCCT', 'ACCCG', 'ACCCC', 'ACCAT', 'ACCAG', 'ACATG', 'ACATC', 'ACAGT',
                           'ACAGG', 'ACAGC', 'ACACT', 'ACACG', 'ACACC', 'AATTC', 'AATGT', 'AATGG', 'AATGC', 'AATCT',
                           'AATCG', 'AATCC', 'AATAT', 'AATAG', 'AATAC', 'AAGTG', 'AAGTC', 'AAGGT', 'AAGGG', 'AAGGC',
                           'AAGCT', 'AAGCG', 'AAGCC', 'AAGAT', 'AAGAG', 'AAGAC', 'AACTT', 'AACTG', 'AACTC', 'AACGT',
                           'AACGG', 'AACGC', 'AACCT', 'AACCG', 'AACCC', 'AACAT', 'AACAG', 'AACAC', 'AAATT', 'AAATG',
                           'AAATC', 'AAAGT', 'AAAGG', 'AAAGC', 'AAACT', 'AAACG', 'AAACC', 'AAAAT', 'AAAAG', 'AAAAC',
                           'CGCGG', 'CCGGG', 'CGGGG', 'GGCAT', 'GCGAT', 'CGGAT', 'GGGAT', 'GATAT', 'CCCCT', 'GCCCT',
                           'CGCCT', 'GGCCT', 'ATCCT', 'GAGCT', 'CCGCT', 'GCGCT', 'AGGCT', 'CGGCT', 'GGGCT', 'ATGCT',
                           'CATCT', 'GATCT', 'CCTCT', 'GCTCT', 'CCAGT', 'GCAGT', 'AGAGT', 'CGAGT', 'GGAGT', 'ATAGT',
                           'CTAGT', 'GACGT', 'CCCGT', 'GCCGT', 'AGCGT', 'CGCGT', 'GGCGT', 'ATCGT', 'CTCGT', 'CAGGT',
                           'GAGGT', 'ACGGT', 'CCGGT', 'GCGGT', 'AGGGT', 'CGGGT', 'GGGGT', 'ATGGT', 'CTGGT', 'CATGT',
                           'GATGT', 'ACTGT', 'CCTGT', 'GCTGT', 'AGTGT', 'CGTGT', 'GGTGT', 'GAATT', 'ACATT', 'CCATT',
                           'GCATT', 'AGATT', 'CGATT', 'GGATT', 'ATATT', 'CTATT', 'GTATT', 'CACTT', 'GACTT', 'ACCTT',
                           'CCCTT', 'GCCTT', 'AGCTT', 'CGCTT', 'GGCTT', 'ATCTT', 'CTCTT', 'GTCTT', 'AAGTT', 'CAGTT',
                           'GAGTT', 'ACGTT', 'CCGTT', 'GCGTT', 'AGGTT', 'CGGTT', 'GGGTT', 'ATGTT', 'CTGTT', 'GTGTT',
                           'AATTT', 'CATTT', 'GATTT', 'ACTTT', 'CCTTT', 'GCTTT', 'AGTTT', 'CGTTT', 'GGTTT', 'ATTTT',
                           'CTTTT', 'GTTTT', 'CCGG', 'CCCG', 'ATGC', 'ATCG', 'ATCC', 'AGGG', 'AGGC', 'AGCG', 'AGCC',
                           'ACGG', 'AATT', 'AATG', 'AATC', 'AGAT', 'ACTG', 'ACTC', 'ACGT', 'ACGC', 'ACCT', 'ACCG',
                           'ACCC', 'ACAT', 'ACAG', 'AAGT', 'AAGG', 'AAGC', 'AACT', 'AACG', 'AACC', 'AAAT', 'AAAG',
                           'AAAC', 'CCGG', 'CGGG', 'GCAT', 'CGAT', 'GGAT', 'CCCT', 'GCCT', 'CGCT', 'GGCT', 'CCGT',
                           'AATT', 'CATT', 'GATT', 'ATCT', 'CAGT', 'GAGT', 'ACGT', 'GCGT', 'AGGT', 'CGGT', 'GGGT',
                           'ATGT', 'CTGT', 'ACTT', 'CCTT', 'GCTT', 'AGTT', 'CGTT', 'GGTT', 'ATTT', 'CTTT', 'GTTT',
                           'CCG', 'ATC', 'AGG', 'AGC', 'ACT', 'ACG', 'ACC', 'AAT', 'AAG', 'AAC', 'CGG', 'GAT', 'CCT',
                           'GCT', 'AGT', 'CGT', 'GGT', 'ATT', 'CTT', 'GTT', 'CG', 'AT', 'AG', 'AC', 'CG',
                           'AT', 'CT', 'GT', 'A', 'C', 'T', 'G']
        self._generate_regex()
        return self.regex
