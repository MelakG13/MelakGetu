import os
import Levenshtein
import matplotlib.pyplot as plt
import logging

class LevenshteinDistanceAnalyzer:
    def __init__(self, file_path):
        self.file_path = file_path
        self.let7_sequences = {}
        self.let7_frequency = {}
        self.avg_distances = []
        self.frequencies = []
        self.miRNA_families = []
        self.logger = logging.getLogger('LevenshteinDistanceAnalyzer')
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def average_levenshtein_distance(self):
        if not os.path.isfile(self.file_path):
            self.logger.error("The specified file does not exist.")
            return

        with open(self.file_path, 'r') as file:
            let7_code = ""
            sequence = ""
            for line in file:
                if line.startswith('>'):
                    header = line[1:].strip()
                    let7_code = self.extract_let7_code(header)
                    sequence = ""
                else:
                    sequence = line.strip()

                if let7_code and sequence:
                    self.let7_sequences.setdefault(let7_code, []).append(sequence)
                    self.let7_frequency[let7_code] = self.let7_frequency.get(let7_code, 0) + 1

        total_sequences_count = 0

        for let7_code, sequences in self.let7_sequences.items():
            total_distance = 0
            total_pairs = 0

            if len(sequences) < 2:
                continue

            for i in range(len(sequences) - 1):
                for j in range(i + 1, len(sequences)):
                    total_distance += Levenshtein.distance(sequences[i], sequences[j])
                    total_pairs += 1

            average_distance = total_distance / total_pairs
            frequency = self.let7_frequency.get(let7_code, 0)
            self.avg_distances.append(average_distance)
            self.frequencies.append(frequency)
            self.miRNA_families.append(let7_code)
            self.logger.info("The Average Levenshtein distance among all pairs for miRNA family %s: %.2f", let7_code, average_distance)
            self.logger.info("The frequency of miRNA family %s: %d", let7_code, frequency)
            total_sequences_count += frequency

        self.logger.info("Total sequences count: %d", total_sequences_count)

        self.plot_results()

    def extract_let7_code(self, header):
        if 'let-7' in header:
            code = header.split('let-7')[1].strip()[0]
            return f"let-7{code}"
        return ""

    def plot_results(self):
        x_pos = range(len(self.miRNA_families))
        fig, ax1 = plt.subplots()
        ax1.bar(x_pos, self.avg_distances, align='center', alpha=0.5)
        ax1.set_ylabel('Average Levenshtein Distance')
        ax1.set_title('Average Levenshtein Distance and Frequency for miRNA Families')

        ax2 = ax1.twinx()
        ax2.plot(x_pos, self.frequencies, 'r')
        ax2.set_ylabel('Frequency')

        plt.xticks(x_pos, self.miRNA_families)
        plt.xlabel('miRNA Family')

        plt.show()


# Set up the logging configuration
logging.basicConfig(level=logging.INFO)

# File path
file_path = r"C:\Users\melak\Advancedprogramming\exam\mature.fa"

# Create an instance of the LevenshteinDistanceAnalyzer class and calculate the average distances
analyzer = LevenshteinDistanceAnalyzer(file_path)
analyzer.average_levenshtein_distance()
