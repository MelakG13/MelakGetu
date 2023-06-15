import Levenshtein
import matplotlib.pyplot as plt
import logging

class LevenshteinDistanceAnalyzer:
    def __init__(self, file_path):
        self.file_path = file_path
        self.let7a_mirnas = []
        self.distances = []
        self.logger = logging.getLogger('LevenshteinDistanceAnalyzer')
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def extract_let7a_mirnas(self):
        with open(self.file_path, "r") as file:
            miRNA_name = ""
            sequence = ""
            for line in file:
                if line.startswith(">"):
                    if miRNA_name != "":
                        if miRNA_name.startswith("hsa-let-7"):
                            self.let7a_mirnas.append((miRNA_name, sequence))
                    miRNA_name = line.strip()[1:]
                    sequence = ""
                else:
                    sequence += line.strip()

            if miRNA_name.startswith("-let-7"):
                self.let7a_mirnas.append((miRNA_name, sequence))

    def calculate_pairwise_distances(self):
        for mirna, seq in self.let7a_mirnas:
            self.logger.info("Pairwise Levenshtein distances for %s", mirna)
            for other_mirna, other_seq in self.let7a_mirnas:
                if mirna != other_mirna:
                    distance = Levenshtein.distance(seq, other_seq)
                    self.distances.append(distance)
                    self.logger.info("Levenshtein distance between %s and %s: %d", mirna, other_mirna, distance)

    def generate_histogram_plot(self):
        plt.hist(self.distances, bins=10)  # Adjust the number of bins as needed
        plt.xlabel("Levenshtein Distance")
        plt.ylabel("Frequency")
        plt.title("Distribution of Levenshtein Distances")
        plt.show()

    def run_analysis(self):
        self.extract_let7a_mirnas()
        self.calculate_pairwise_distances()
        self.generate_histogram_plot()


# Set up the logging configuration
logging.basicConfig(level=logging.INFO)

# Specify the file path to the mature.fa file on your computer
file_path = r"C:\Users\melak\Advancedprogramming\exam\mature.fa"

# Create an instance of the LevenshteinDistanceAnalyzer class and run the analysis
analyzer = LevenshteinDistanceAnalyzer(file_path)
analyzer.run_analysis()
