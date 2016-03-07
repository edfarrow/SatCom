__author__ = 'ewjf201'
import sys
import os


def is_odd(num):
    return num & 0x1

print(os.getcwd())
#output_file_name = raw_input("Please enter the filename to output to: ")
output_file_name = "test-fastqgen3.fastq"
try:
    output_file = open(output_file_name, mode="wb")
except IOError as e:
    print("There was an error opening the file. Exiting!")
    sys.exit()

output = []
seq_length = 100
insert_point = ["none", "start", "middle", "end"]

while True:
    sequence = ""
    sequence_chars = raw_input("Enter the sequence you like to repeat, or False to end: ")
    if sequence_chars == 'False':
        break
    sequence_chars = sequence_chars.upper()
    substitution_raw = raw_input("Enter a sequence you would like to substitute: ").upper()

    repeats = len(sequence_chars)
    for i in range(0, seq_length):
        if len(sequence) < 100:
            sequence += sequence_chars
        elif len(sequence) is 100:
            break
        elif len(sequence) > 100:
            cut_len = len(sequence) - 100
            sequence = sequence[:-cut_len]

    sequence_pres = sequence

    if substitution_raw != "":
        for i in range(1,10):

            substitution = substitution_raw * i
            for sub in insert_point:
                sequence = sequence_pres
                substitution_len = len(substitution)
                if sub == "middle":

                    if is_odd(substitution_len):
                        if substitution_len != 1:
                            sub_len_1 = float(substitution_len) / 2
                            sub_len_1 += 0.5
                            sub_len_1 = int(sub_len_1)
                            sub_len_2 = float(substitution_len) / 2
                            sub_len_2 -= 0.5
                            sub_len_2 = int(sub_len_2)
                        else:
                            sub_len_1 = 1
                            sub_len_2 = 0

                        sequence = sequence[:50-sub_len_1] + substitution + sequence[50+sub_len_2:]
                    else:
                        sub_len_1 = substitution_len / 2
                        sub_len_2 = sub_len_1

                        sequence = sequence[:50-sub_len_1] + substitution + sequence_pres[50+sub_len_2:]

                elif sub == "start":
                    sequence = substitution + sequence[substitution_len:]

                elif sub == "end":
                    sequence = sequence[:100-substitution_len] + substitution

                if substitution_len > 1 and sub != "none":
                    print(sequence)
                    output.append(sequence)
                elif substitution_len == 1:
                    print(sequence)
                    output.append(sequence)
    else:
        print(sequence)
        output.append(sequence)
test_i = 0
for out in output:
    test_i += 1
    buffer_out = []
    output_file.write("@SAMPLESEQUENCE %d\n" % test_i)
    output_file.write(out + "\n")
    output_file.write("+\n")
    for i in range(0, 100):
        buffer_out.append("~")
    output_file.write(''.join(buffer_out))
    output_file.write("\n")

output_file.close()



