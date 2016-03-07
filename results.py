__author__ = 'ewjf201'
import os
import subprocess

test_file_1 = "/path/to/1250-fastq.fastq"
test_file_2 = "/path/to/10000-fastq.fastq"
test_file_3 = "/path/to/100000-fastq.fastq"
process = [0, 1, 2, 3, 4, 5, 6, 7, 8]
fuzzy = [True, False]
test_files = [test_file_1]

# Run tests sequentially
curdir = os.path.dirname(os.path.realpath(__file__))
os.chdir(curdir)
files_1 = test_files
for files in files_1:
    print("*** Testing on file %s ***" % files)
    for processes in process:
        print("\n\n*** Test cores: %d ***" % processes)
        for repeat in range(1, 5):
            print("\n\nRepeat count: %d" % repeat)
            if processes is not 0:
                for f in fuzzy:
                    if f is True:
                        print("\nFuzzy results:")
                        proc = subprocess.Popen(["/path/to/python", "main.py", "-m",
                                                 "-c", str(processes), "-f", files])
                        proc.communicate()
                    else:
                        print("\nNon-fuzzy results:")
                        proc = subprocess.Popen(["/path/to/python", "main.py", "-m",
                                                 "-c", str(processes), files])
                        proc.communicate()
            else:
                for f in fuzzy:
                    if f is True:
                        print("\nFuzzy results:")
                        proc = subprocess.Popen(["/path/to/python", "main.py", "-f", files])
                        proc.communicate()
                    else:
                        print("\nNon-fuzzy results:")
                        proc = subprocess.Popen(["/path/to/python", "main.py", files])
                        proc.communicate()