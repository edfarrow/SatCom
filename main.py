from __future__ import print_function
from satcom.file_loader import SatcomFile
from satcom.MicrosatelliteRegex import MicrosatelliteRegex
import os
import argparse
from argparse import Namespace
import logging
from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import Queue
import regex
import time
import sys
import itertools
import traceback
from collections import OrderedDict

parser = argparse.ArgumentParser(description='Analyses FASTQ files for microsatellites.')
parser.add_argument('-m', '--multithread', action="store_true", help="Multithread the search (may perform worse)")
parser.add_argument('-o', '--output-file', nargs=1, help="Specify where to output the results to.")
parser.add_argument('-f', '--fuzzy', action="store_true", help="If specified, fuzzes matches in Regex.")
parser.add_argument('-s', '--substitutions', nargs=1, help="If -f is specified, sets the fuzziness. "
                                                           "Default is 2 substitutions allowed.", default=2)
parser.add_argument('-c', '--cpus', nargs=1, type=int, help="The number of CPUs to use")
parser.add_argument('-v', '--verbose', action="store_true", help="Do you want debug output?")
parser.add_argument('FILE', help="FastQ input file (including path)")
logging.info("Parsing arguements provided:")
args = parser.parse_args()
logging.info(args)
if args.verbose is True:
    print(args)

logging.info("Starting timer")
start_time = time.time()


def _check_match(match):
    if len(match) is not 0:
        matches_1 = match
        return list(matches_1)

def regex_worker(sequence):
    """
    :param sequence:
    :return:
    """
    matches = regex.match(regex_pregen, str(sequence))
    if matches is not None:
        match_capt = matches.groupdict()
        matches = [match_capt['whole'], match_capt['start'], match_capt['msat'], match_capt['end']]
    return matches


def regex_worker_multithread(sequence):
    matches = regex.match(regex_pregen, str(sequence))
    if matches is not None:
        match_capt = matches.groupdict()
        matches = [match_capt['whole'], match_capt['start'], match_capt['msat'], match_capt['end']]
    return matches

# Exec search
logging.info("Loading FASTQ")
fastq_file = SatcomFile(args.FILE, 'fastq')


logging.info("Loading Regex")
start_regex_time = time.time()

min_seq_len = [14, 7, 5, 5, 3]

if args.fuzzy is not None:
    if args.fuzzy is True:
        regex_search = MicrosatelliteRegex(fuzzy_substitution=args.substitutions, min_seq_length=min_seq_len)
        fuzzy = True
    else:
        regex_search = MicrosatelliteRegex(fuzzy_substitution=0,  min_seq_length=min_seq_len)
        fuzzy = False
else:
    regex_search = MicrosatelliteRegex(fuzzy_substitution=0,  min_seq_length=min_seq_len)
    fuzzy = False

regex_pregen = regex_search.create_regex(preloaded=True)
logging.info("Starting searches")

if __name__ == '__main__':
    result_1 = []
    output_1 = []

    if args.multithread is True:
        p = Pool(nodes=args.cpus[0])
        res = p.amap(regex_worker_multithread, fastq_file.sequences)
        count = 0
        while not res.ready():
            count += 2
            print("\rWaiting. Timer: %d" % count, end='')
            time.sleep(2)
        result_1 = res.get()
        print("\nDone searching")

    elif args.multithread is False:
        num_seq = 0
        for seq in fastq_file.sequences:
            num_seq += 1
            rexex = regex_worker(seq.seq)
            if rexex:
                result_1.append(rexex)
            print("\rCompleted %d searches!" % num_seq, end='')
    else:
        logging.critical("No valid multithread arguement!")
        exit()
    output = result_1

    print("\nFinished Regex search")

logging.info("Finished searches")
if args.verbose is True:
    print("*** Microsatellites Found ***")
    for microsatellite in output:
        print("Full sequence: %s" % microsatellite[0], end=" ")
        print("Start sequence: %s" % microsatellite[1], end=" ")
        print("Microsat sequence: %s" % microsatellite[2], end=" ")
        print("End sequence: %s" % microsatellite[3])
if args.output_file is not None and output != []:
    with open(args.output_file[0], "wb") as out_file:
        out_file.write("*** SatCom output file ***\n")
        out_file.write("*** Authored by Ed Farrow, available at github.com\n")
        for microsatellite in output:
            out_file.write("Full sequence: %s\n" % microsatellite[0])
            out_file.write("Start sequence: %s\n" % microsatellite[1])
            out_file.write("Microsat sequence: %s\n" % microsatellite[2])
            out_file.write("End sequence: %s\n" % microsatellite[3])
            out_file.write("\n")
if args.verbose is True:
    print("*** End Microsatellites ***\n")

print("*** Found %d microsatellites ***" % len(output))
print("Finished analysing microsatellites in %.3f seconds!" % (time.time() - start_time))
logging.info("Finished!")
