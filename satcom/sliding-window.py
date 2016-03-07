__author__ = 'ewjf201'
import numpy as np

class MicrosatelliteString:
    def __init__(self, sequence_length=18, bases=("A", "C", "T", "G"), min_microsat = 1, max_microsat = 6, save=False):
        self.seq_len = sequence_length
        self.bases = bases
        self.save = save
        self.min_microsat = min_microsat
        self.max_microsat = max_microsat

    def _generate_strings(self):
        self.generated_strings = np.array([], np.string_)
        min_repeats = self.seq_len /