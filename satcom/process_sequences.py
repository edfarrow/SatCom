from __future__ import print_function


class SequenceProcessor:
    import regex
    from satcom import MicrosatelliteRegex
    from satcom.file_loader import SatcomFile
    import logging
    import time
    import sys
    import itertools
    from pathos.multiprocessing import ProcessingPool as Pool
    def __init__(self, fastq_file, fuzzy=True, multithread=True, cpus=8, ):




def regex_worker(sequence):
    """
    Uses a loop through the pregenerated regex instead of through sequences because it's easier to keep track of
    :param sequence:
    :return:
    """
    # matches = []
    # for regex_indiv in regex_pregen:
    #     matches_1 = regex.findall(regex_indiv, str(sequence))
    #     if len(matches_1) is not 0:
    #         if matches_1[0] is not "" or " ":
    #             for match in matches_1:
    #                 matches.append([match[0], match[1], match[2], match[3], regex_indiv])
    # return matches

    def _check_match(match):
        if len(match) is not 0:
            if match[0] is not "" or " ":
                return [match[0], match[1], match[2], match[3]]
    matches = [[_check_match(m) for m in regex.findall(r, str(sequence))] for r in regex_pregen]
    return matches


def regex_worker_multithread(sequence):
    # try:
    def _check_match(match):
        if len(match) is not 0:
            if match[0] is not "" or " ":
                return [match[0], match[1], match[2], match[3]]
    matches = [[_check_match(m) for m in regex.findall(r, str(sequence.seq))] for r in regex_pregen]
    # for regex_indiv in regex_pregen:
    #     matches_1 = regex.findall(regex_indiv, str(sequence.seq))
    #     if len(matches_1) is not 0:
    #         if matches_1[0] is not "" or " ":
    #             for match in matches_1:
    #                 matches.append([match[0], match[1], match[2], match[3]])
    return matches
    # except Exception as e:
    #     for frame in traceback.extract_tb(sys.exc_info()[2]):
    #         fname,lineno,fn,text = frame
    #         errors.put(["Error in %s on line %d" % (fname, lineno)])
    #     errors.put(e)


# Exec search
logging.info("Loading FASTQ")
fastq_file = SatcomFile(args.FILE)

if args.regex is not None:
    logging.info("Loading Regex")
    start_regex_time = time.time()

    if args.fuzzy is not None:
        if args.fuzzy is True:
            regex_search = MicrosatelliteRegex(fuzzy=True)
            fuzzy = True
        else:
            regex_search = MicrosatelliteRegex(fuzzy=False)
            fuzzy = False
    else:
        regex_search = MicrosatelliteRegex(fuzzy=False)
        fuzzy = False

    regex_pregen = regex_search.create_regex(preloaded=True)
    if args.save_regex is not None:
        regex_search.save_regex(args.save_regex[0])
    elif args.existing_regex is not None:
        regex_search.load_regex(args.existing_regex[0])
    print("Loaded %d regex in %d seconds!" % (len(regex_search.regex), time.time() - start_regex_time))
    start_regex_time = None
    logging.info("Starting searches")

    #results = search_regex(regex_worker, regex_res, fastq_file, args.cpus)
    if __name__ == '__main__':
        result_1 = []
        output_1 = []

        if args.multithread is True:
            # errors = Queue()
            p = Pool(nodes=args.cpus[0])
            res = p.amap(regex_worker_multithread, fastq_file.sequences)
            count = 0
            while not res.ready():
                count += 2
                print("\rWaiting. Timer: %d" % count, end='')
                time.sleep(2)
            # for e in errors.get():
            #     print(e)
            # p.terminate()
            result_1 = res.get()
            print("\nDone searching")

        elif args.multithread is False:
            num_seq = 0
            for seq in fastq_file.sequences:
                num_seq += 1
                rexex = regex_worker(seq.seq)
                result_1.append(rexex)
                print("\rCompleted %d searches!" % num_seq, end='')
        else:
            logging.critical("No valid multithread arguement!")
            exit()

        for result in result_1:
            if result is not []:
                for res in result:
                    output_1.append([res[1], res[2], res[3]])
        # remove duplicates
        output = [k for k, _ in itertools.groupby(output_1.sort())]


        print("\nFinished Regex search")

    logging.info("Finished searches")
    if args.verbose is True:
        print("*** Microsatellites Found ***")
        for microsatellite in output:
            print("Start sequence: " + microsatellite[0], end=" ")
            print("Microsatellite sequence: " + microsatellite[1], end=" ")
            print("End sequence: " + microsatellite[2])
    if args.output_file is not None and output[0] is not None:
        with open(args.output_file[0], "wb") as out_file:
                out_file.write("Start sequence: " + microsatellite[0] + "\n")
                out_file.write("Microsatellite sequence: " + microsatellite[1] + "\n")
                out_file.write("End sequence: " + microsatellite[2] + "\n")
                out_file.write("\n")
    if args.verbose is True:
        print("*** End Microsatellites ***")

else:
    logging.info("Loading Regex")
    start_regex_time = time.time()

    if args.fuzzy is not None:
        if args.fuzzy is True:
            regex_search = MicrosatelliteRegex(fuzzy=True)
            fuzzy = True
        else:
            regex_search = MicrosatelliteRegex(fuzzy=False)
            fuzzy = False
    else:
        regex_search = MicrosatelliteRegex(fuzzy=False)
        fuzzy = False