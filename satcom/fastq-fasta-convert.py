__author__ = 'ewjf201'
from Bio import SeqIO

with open("/Users/Ed/PycharmProjects/Bio-Project-GUI/data/100k-fastq.fasta", "w") as output:
    with open("/Users/Ed/PycharmProjects/Bio-Project-GUI/data/100k-fastq.fastq", "r") as input_fastq:
        sequences = SeqIO.parse(input_fastq, "fastq")
        SeqIO.write(sequences, output, "fasta")