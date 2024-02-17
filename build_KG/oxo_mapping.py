""" 
This script is modified from https://github.com/WengLab-InformaticsResearch/oxo_py/blob/master/OxO.py. Thanks WengLab-InformaticsResearch for providing this code.

Basic implementation of EMBL-EBI's OxO mappings (https://www.ebi.ac.uk/spot/oxo/index) in Python using mapping files provided by EMBL-EBI.

Notes:
    The mapping algorithm may be different from OxO's. A few test cases have been run and produced matching results.

Examples:
    OxO.find_mappings('DOID:162')
    OxO.find_mappings('UMLS:C0002199', distance=3)
    OxO.find_mappings('SNOMEDCT:136111001', distance=3, targets=['MeSH', 'UMLS'])
"""

## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
import csv
from collections import defaultdict
from typing import List, Dict, Tuple, Set, Union, Optional

## Import custom libraries
from utils import get_logger

class OxO:

    def __init__(self, mapping_file_dir: str):
        # Initialize logger
        self.logger = get_logger()

        # set up mapping file paths
        self._file_ols = os.path.join(mapping_file_dir, 'ols_mappings.csv')
        if not os.path.exists(self._file_ols):
            self.logger.error(f'Could not find OLS mapping file at {self._file_ols}')
            return
        
        self._file_umls = os.path.join(mapping_file_dir, 'umls_mappings.csv')
        if not os.path.exists(self._file_umls):
            self.logger.error(f'Could not find UMLS mapping file at {self._file_umls}')
            return
        
        self._file_ols_terms = os.path.join(mapping_file_dir, 'ols_terms.csv')
        if not os.path.exists(self._file_ols_terms):
            self.logger.error(f'Could not find OLS terms file at {self._file_ols_terms}')
            return
        
        self._file_umls_terms = os.path.join(mapping_file_dir, 'umls_terms.csv')
        if not os.path.exists(self._file_umls_terms):
            self.logger.error(f'Could not find UMLS terms file at {self._file_umls_terms}')
            return
        
        # Load mapping files
        self.logger.info(f'Loading OxO mapping files from {mapping_file_dir}')
        self._load_files()

    def _load_files(self):
        # Initialize
        self._terms = dict()
        self._mappings = defaultdict(set)

        # Read in the ols terms
        with open(self._file_ols_terms, 'r', newline='') as fh:
            reader = csv.reader(fh, delimiter=',', quotechar='"', doublequote=False, lineterminator='\r\n', escapechar='\\')
            # Skip the header line
            reader.__next__()
            # Read in term definitions
            for identifier, curie, label, uri, prefix in tqdm(reader, desc='Loading OLS terms'):
                self._terms[curie] = {'label': label, 'uri': uri, 'prefix': prefix, 'identifier': identifier}

        # Read in the umls terms
        with open(self._file_umls_terms, 'r', newline='') as fh:
            reader = csv.reader(fh, delimiter=',', quotechar='"', doublequote=False, lineterminator='\r\n', escapechar='\\')
            # Skip the header line
            reader.__next__()
            # Read in term definitions
            for identifier, curie, label, uri, prefix in tqdm(reader, desc='Loading UMLS terms'):
                self._terms[curie] = {'label': label, 'uri': uri, 'prefix': prefix, 'identifier': identifier}

        # Read in OLS mapping file
        with open(self._file_ols, 'r', newline='') as fh:
            reader = csv.reader(fh, delimiter=',', quotechar='"', doublequote=False, lineterminator='\r\n', escapechar='\\')
            # Skip the header line
            reader.__next__()
            # Read in all mappings
            for row in tqdm(reader, desc='Loading OLS mappings'):
                curie_from = row[0]
                curie_to = row[1]
                self._mappings[curie_from].add(curie_to)
                self._mappings[curie_to].add(curie_from)

        # Read in UMLS mapping file
        with open(self._file_umls, 'r', newline='') as fh:
            reader = csv.reader(fh, delimiter=',', quotechar='"', doublequote=False, lineterminator='\r\n', escapechar='\\')
            # Skip the header line
            reader.__next__()
            # Read in all mappings
            for row in tqdm(reader, desc='Loading UMLS mappings'):
                curie_from = row[0]
                curie_to = row[1]
                self._mappings[curie_from].add(curie_to)
                self._mappings[curie_to].add(curie_from)

    def find_mappings(self, curie_source: Union[List, str], distance: int = 1, targets: Union[List, str] = None):
        self.logger.info(f'Finding mappings for {curie_source} with distance {distance} and targets {targets}')

        if not hasattr(self, '_terms') or not hasattr(self, '_mappings'):
            self._load_files()

        found = dict()  # mapping results (key:curie, value:distance)
        visited = set()  # nodes already visited
        # nodes to visit on this iteration
        if isinstance(curie_source, str):
            searching = set([curie_source])
        elif type(curie_source) is list:
            searching = set(curie_source)
        prefix_source = [curie.split(':')[0] for curie in searching]

        # Convert targets to a set
        if targets is None:
            targets = set()
        elif type(targets) is str:
            targets = set([targets])
        else:
            targets = set(targets)

        for i in range(distance):
            search_add = set()  # nodes to search in the next iteration

            # Mark all nodes that we're about to visit as already visited
            visited = visited.union(searching)

            # Visit each new node
            for curie in tqdm(searching, desc=f'Searching with distance={i}'):
                curr_mappings = self._mappings[curie]

                # Add new mappings to the set to search in the next iteration if we have not already visited
                search_add = search_add.union([x for x in curr_mappings if x not in visited])

                # Add new mappings to the set of found mappings if it's in the target ontologies
                for m in curr_mappings:
                    prefix_curr = m.split(':')[0]
                    if m not in found and prefix_curr not in prefix_source and (len(targets) == 0 or prefix_curr in targets):
                        info = {'distance': i + 1, 'label': '', 'uri': ''}
                        if m in self._terms:
                            term = self._terms[m]
                            info['label'] = term['label']
                            info['uri'] = term['uri']
                        found[m] = info

            searching = search_add

        return found


if __name__ == '__main__':
    oxo = OxO('oxo_mapping_data/oxo-mappings-2020-02-04/')
    oxo.find_mappings('EFO:1001950', distance=1, targets=['MeSH','ICD9CM','ICD10CM','ICD11CM','OMIM'])