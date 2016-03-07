__author__ = 'ewjf201'
from Bio import SeqIO


class SatcomFile:
    def __init__(self, file_name, file_type):
        """

        :param file_name:
        :param file_type:
        :return:
        """
        self.file_name = file_name
        self.file_type = file_type
        self.sequences = self._load_file()

    def _open_file(self):
        try:
            file_handle = open(self.file_name)
        except IOError as e:
            return False
        return file_handle

    def _load_seq(self, file_handle, file_type):
        try:
            sequences = SeqIO.parse(file_handle, file_type)
        except ValueError as e:
            return False
        return sequences

    def _load_file(self):
        for i in range(0,3):
            file_handle = self._open_file()
            if file_handle is not False:
                break
            if i >= 3:
                return False

        for i in range(0,3):
            seq_gen = self._load_seq(file_handle, self.file_type)
            if file_handle is not False:
                return seq_gen
            if i >= 3:
                return False
