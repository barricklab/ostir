#!/usr/bin/env python
"""
This module contains an early rewrite of the OSTIRFactory class. It is not yet
implemented, but it will be used to replace the OSTIRFactory class in the future.

Parrellism is more straightforward to implement in function based code, so this
module turns the class based legacy code into a function based implementation.

Some of the codebase of this file originated from the Salis Lab's RBS calculator which is distributed under GPL3.
See <http://www.gnu.org/licenses/>.
Copyright 2008-2009 is owned by the University of California Regents. All rights reserved.
"""

raise NotImplementedError("This module is not yet implemented")

from dataclasses import dataclass

@dataclass
class OSTIRResult():
    """Class for storing the results of an OSTIR calculation"""
    # These are the results that most people care about
    name: str = ""
    start_codon: str = ""
    start_position: int = 0
    expression: float = 0.0
    RBS_distance_bp: int = 0
    dG_total: float = 0.0
    dG_rRNA_mRNA: float = 0.0
    dG_mRNA: float = 0.0
    dG_spacing: float = 0.0
    dG_standby: float = 0.0
    dG_start_codon: float = 0.0

    # These are the results that are not as important as the above
    mRNA_structure: str = ""
    mRNA_rRNA_uncorrected_structure: str = ""
    mRNA_rRNA_corrected_structure: str = ""
    most_5p_mRNA: str = ""
    longest_helix: str = ""
    dG_start_energy: float = 0.0
    helical_loop: str = ""
    bulge_loop: str = ""
    min_bp_prob: float = 0.0
    kinetic_score: float = 0.0

    def __getitem__(self, key):
        if key in self.__annotations__:
            return getattr(self, key)
        raise KeyError(key)

    def results(self):
        """Returns a dictionary of the core results"""
        return {
            'name': self.name,
            'start_codon': self.start_codon,
            'start_position': self.start_position,
            'expression': self.expression,
            'RBS_distance_bp': self.RBS_distance_bp,
            'dG_total': self.dG_total,
            'dG_rRNA:mRNA': self.dG_rRNA_mRNA,
            'dG_mRNA': self.dG_mRNA,
            'dG_spacing': self.dG_spacing,
            'dG_standby': self.dG_standby,
            'dG_start_codon': self.dG_start_codon,
        }

    def addional_results(self):
        """Returns a dictionary additional results"""
        return {
            'mRNA_structure': self.mRNA_structure,
            'mRNA_rRNA_uncorrected_structure': self.mRNA_rRNA_uncorrected_structure,
            'mRNA_rRNA_corrected_structure': self.mRNA_rRNA_corrected_structure,
            'most_5p_mRNA': self.most_5p_mRNA,
            'longest_helix': self.longest_helix,
            'dG_start_energy': self.dG_start_energy,
            'helical_loop': self.helical_loop,
            'bulge_loop': self.bulge_loop,
            'min_bp_prob': self.min_bp_prob,
            'kinetic_score': self.kinetic_score,
        }


class OstirManager():
    """This is the class that accepts sequences, prepares parallelization, and sends parts out for computation.
    Everything here is managed in python, while the actual computation is sent off to cython"""
    def __init__(self):
        self.start_range_1 = None
        self.start_codons = None

    @staticmethod
    def _job() -> OSTIRResult:
        
        pass

    def run(self) -> list:
        """Runs the OSTIR calculation on the submitted sequences"""
        
        
    def find_start_codons(self, sequence):
        """Finds all start codons in an mRNA sequence. Creates a list."""

        self.start_position_list = []

        seq_len = len(sequence)

        #Switch to zero-indexed positions happens here
        end_0 = min(self.start_range_1[1]-1, seq_len - 2)
        begin_0 = min(self.start_range_1[0]-1, end_0)

        for i in range(begin_0, end_0 + 1):
            codon = sequence[i:i + 3]
            if codon.upper() in self.start_codons:
                self.start_position_list.append(i)
                yield (i, codon)
            else:
                pass
    
    def submit_sequence_str(self):
        pass

    def submit_sequence_file(self):
        pass
