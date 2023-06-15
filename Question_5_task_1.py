import os
from Levenshtein import distance
import logging

class LevenshteinDistanceCalculator:
    def __init__(self, file_path):
        self.file_path = file_path
        self.let7_sequences = []
        self.total_let7_miRNA = 0
        self.total_distance = 0
        self.total_pairs = 0
        self.logger = logging.getLogger('LevenshteinDistanceCalculator')
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def calculate_levenshtein_distance(self):
        if not os.path.isfile(self.file_path):
            self.logger.error("The specified file does not exist.")
            return

        with open(self.file_path, 'r') as file:
            species = ""
            sequence = ""
            for line in file:
                if line.startswith('>'):
                    header = line[1:].strip()
                    species = self.extract_species(header)
                    sequence = ""
                else:
                    sequence = line.strip()

                if 'let-7a' in header and species and sequence:
                    self.let7_sequences.append((species, sequence))
                    self.total_let7_miRNA += 1

        for i in range(len(self.let7_sequences) - 1):
            for j in range(i + 1, len(self.let7_sequences)):
                species1, seq1 = self.let7_sequences[i]
                species2, seq2 = self.let7_sequences[j]
                levenshtein_distance = distance(seq1, seq2)
                self.total_distance += levenshtein_distance
                self.total_pairs += 1
                self.logger.info("Species: %s - %s | Levenshtein Distance: %d", species1, species2, levenshtein_distance)

        if self.total_pairs > 0:
            average_distance = self.total_distance / self.total_pairs
            self.logger.info("Total 'let-7a' miRNAs: %d", self.total_let7_miRNA)
            self.logger.info("Average Levenshtein Distance of 'let-7a' miRNAs: %.2f", average_distance)

    def extract_species(self, header):
        species = header.split(' ')[0]
        return species


# Set up the logging configuration
logging.basicConfig(level=logging.INFO)

# File path
file_path = r"C:\Users\melak\Advancedprogramming\exam\mature.fa"

# Create an instance of the LevenshteinDistanceCalculator class and calculate the distances
calculator = LevenshteinDistanceCalculator(file_path)
calculator.calculate_levenshtein_distance()
