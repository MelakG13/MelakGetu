import logging

class MiRNACounter:
    def __init__(self, file_path):
        self.file_path = file_path
        self.let7_count = 0
        self.logger = logging.getLogger('MiRNACounter')
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def count_let7_miRNAs(self):
        with open(self.file_path, "r") as file:
            for line in file:
                if line.startswith(">"):
                    miRNA_name = line.strip()[1:]
                    if "let-7" in miRNA_name:
                        self.let7_count += 1

    def run_analysis(self):
        self.count_let7_miRNAs()
        self.logger.info("Total number of let-7 miRNAs across all species: %d", self.let7_count)

# Set up the logging configuration
logging.basicConfig(level=logging.INFO)

# Specify the file path to the mature.fa file on your computer
file_path = r"C:\Users\melak\Advancedprogramming\exam\mature.fa"

# Create an instance of the MiRNACounter class and run the analysis
counter = MiRNACounter(file_path)
counter.run_analysis()
