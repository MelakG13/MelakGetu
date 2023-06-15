import logging
import matplotlib.pyplot as plt
import numpy as np
import os

class MiRNAAnalyzer:
    def __init__(self, file_path):
        self.file_path = file_path
        self.let7_species = {}
        self.logger = logging.getLogger('MiRNAAnalyzer')
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def extract_let7_code(self, header):
        if 'let-7' in header:
            code = header.split('let-7')[1].split()[0]
            return f"let-7{code[0]}"
        return ""

    def extract_species_code(self, header):
        start_index = header.find('>') + 1
        end_index = header.find('-', start_index)
        if start_index < end_index:
            return header[start_index:end_index]
        return ""

    def process_file(self):
        if not os.path.isfile(self.file_path):
            self.logger.error("The specified file does not exist.")
            return

        current_species = ""

        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    header = line[1:].strip()
                    current_species = self.extract_species_code(header)
                else:
                    let7_code = self.extract_let7_code(header)
                    if let7_code:
                        self.let7_species.setdefault(current_species, {}).setdefault(let7_code, 0)
                        self.let7_species[current_species][let7_code] += 1

    def plot_let7_family_presence(self):
        let7_species_filtered = {species: counts for species, counts in self.let7_species.items() if sum(counts.values()) >= 5}

        species_list = list(let7_species_filtered.keys())
        let7_codes = list(set().union(*[d.keys() for d in let7_species_filtered.values()]))

        presence_matrix = np.zeros((len(species_list), len(let7_codes)))

        for i, species in enumerate(species_list):
            for j, let7_code in enumerate(let7_codes):
                presence_matrix[i, j] = let7_species_filtered[species].get(let7_code, 0)

        x = np.arange(len(species_list))
        bar_width = 0.35

        fig, ax = plt.subplots(figsize=(10, 6))

        bottom = np.zeros(len(species_list))

        for i, let7_code in enumerate(let7_codes):
            ax.bar(x, presence_matrix[:, i], bottom=bottom, label=let7_code)
            bottom += presence_matrix[:, i]

        ax.set_xlabel('Species')
        ax.set_ylabel('Presence Count')
        ax.set_title('Presence of let-7 Family per Species')
        ax.set_xticks(x)
        ax.set_xticklabels(species_list, rotation=45, ha='right')
        ax.legend()

        plt.tight_layout()
        plt.show()

    def run_analysis(self):
        self.process_file()
        self.logger.info("Analysis completed successfully.")
        self.plot_let7_family_presence()


# Set up the logging configuration
logging.basicConfig(level=logging.INFO)

# Specify the file path to the mature.fa file on your computer
file_path = r"C:\Users\melak\Advancedprogramming\exam\mature.fa"

# Create an instance of the MiRNAAnalyzer class and run the analysis
analyzer = MiRNAAnalyzer(file_path)
analyzer.run_analysis()
