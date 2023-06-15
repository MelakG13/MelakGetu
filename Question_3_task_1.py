import requests
import re
import logging

class miRBaseVersionFetcher:
    def __init__(self):
        self.url = "https://www.mirbase.org/ftp/CURRENT/README"
        self.logger = logging.getLogger('miRBaseVersionFetcher')
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def get_mirbase_version(self):
        response = requests.get(self.url)
        readme_text = response.text
        # Search for the version information in the README file
        version_match = re.search(r"Release (\d+\.\d+)", readme_text)
        if version_match:
            miRBase_version = version_match.group(1)
            return miRBase_version
        else:
            return None

    def print_mirbase_version(self):
        mirbase_version = self.get_mirbase_version()
        if mirbase_version:
            self.logger.info("Current version of miRBase: %s", mirbase_version)
        else:
            self.logger.info("Unable to retrieve the current version of miRBase.")

# Set up the logging configuration
logging.basicConfig(level=logging.INFO)

# Create an instance of the miRBaseVersionFetcher class
version_fetcher = miRBaseVersionFetcher()

# Call the method to print the miRBase version
version_fetcher.print_mirbase_version()
