#!/usr/bin/env python
"""
This module contains the OSTIRFactory class, which is the main script used to call the ostir calculation functions
and compile a result. There should be no need to interact with this module directly except in advanced use cases.

Some of the codebase of this file originated from the Salis Lab's RBS calculator which is distributed under GPL3.
See <http://www.gnu.org/licenses/>.
Copyright 2008-2009 is owned by the University of California Regents. All rights reserved.
"""


import re
import math
import os
import concurrent.futures
from dataclasses import dataclass
from .ostir_calculations import calc_longest_loop_bulge, calc_longest_helix, calc_dG_standby_site, find_start_codons, calc_dG_mRNA, calc_dG_mRNA_rRNA, calc_expression_level

class CalcError(Exception):
    """Base class for exceptions in this module."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

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
    longest_helix: str = ""
    dG_start_energy: float = 0.0
    helical_loop: str = ""
    bulge_loop: str = ""
    min_bp_prob: float = 0.0
    kinetic_score: float = 0.0

    def __getitem__(self, key):
        return getattr(self, key)

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
            'longest_helix': self.longest_helix,
            'dG_start_energy': self.dG_start_energy,
            'helical_loop': self.helical_loop,
            'bulge_loop': self.bulge_loop,
            'min_bp_prob': self.min_bp_prob,
            'kinetic_score': self.kinetic_score,
        }


class OSTIRFactory:
    """Class for calculating the rate of translation initiation of a given mRNA sequence"""
    def __init__(self, mRNA, start_range, rRNA, constraints, verbose=False, decimal_places=4, circular=False, name="unnamed"):
        """
        Initializes the RBS Calculator class with the mRNA sequence and the range of start codon positions considered.
        start_range is a pair of 1-indexed positions
        """

        # Sets deaults
        # From OSTIR calibration using Salis2009 data. See calibration directory for procedure
        self.Beta = 0.40002512
        self.RT_eff = 1/self.Beta
        self.logK = 7.279194329
        self.K = math.exp(self.logK)

        # Global parameters -- constants
        self.infinity = 1e12  # For all practical purposes, here.
        self.RNA_model = 'rna2004'
        self.start_codon_energies = {"ATG": -1.194, "AUG": -1.194, "GTG": -0.0748, "GUG": -0.0748, "TTG": -0.0435,
                                        "UUG": -0.0435, "CTG": -0.03406, "CUG": -0.03406}  # hybridization to CAT
        self.auto_dangles = True
        self.dangles_default = "all"
        self.temp = 37.0
        self.optimal_spacing = 5  # aligned spacing

        # From OSTIR calibration using Salis2009 data. See calibration directory for procedure
        self.dG_spacing_constant_push = [17.20965071, 3.46341492, 1.790848365, 3.0]
        self.dG_spacing_constant_pull = [0.06422042, 0.275640836, 0.0]
        self.cutoff = 35  # number of nt +- start codon considering for folding
        self.standby_site_length = 4  # Number of nt before SD sequence that must be unpaired for ribosome binding
        self.energy_cutoff = 3.0
        self.start_codons = ["ATG", # substituted U for T in actual calcs. Ignores CTG/CUG
                             "AUG",
                             "GTG",
                             "GUG",
                             "TTG",
                             "UUG"]
        self.footprint = 1000
        """Footprint of the 30S complex that prevents formation of secondary structures
        downstream of the start codon. Here, we assume that the entire post-start RNA
        sequence does not form secondary structures once the 30S complex has bound.
        """
        self.circular = circular

        # NuPACK.__init__(self,sequences,self.RNA_model)
        exp = re.compile('[ATGCU._]', re.IGNORECASE)
        if exp.match(mRNA) is None:
            raise ValueError(f"Invalid letters found in sequence {mRNA}. Only ATGCU accepted.")
        mRNA = mRNA.replace('.', '')
        mRNA = mRNA.replace('_', '')

        if start_range[0] < 1:
            start_range[0] = 1
        if start_range[1] > len(mRNA):
            start_range[1] = len(mRNA)


        self.install_location = os.path.dirname(os.path.realpath(__file__))
        self.name = name
        self.mRNA_input = mRNA.upper()
        self.rRNA = rRNA
        self.constraints = constraints
        self.rRNA_len = len(self.rRNA)
        self.mRNA_len = len(self.mRNA_input)
        self.total_sequence_length = len(mRNA) + len(self.rRNA)
        self.run = 0
        self.start_range = start_range
        self.verbose = verbose
        self.threads = 1
        self.decimal_places = decimal_places
        self.results = []


    def calc_dG(self):
        """Calculates each dG term in the free energy model and sums them together to create dG_total"""

        mRNA_length = len(self.mRNA_input)
        if self.circular:
            if mRNA_length > 200:
                self.mRNA_input = self.mRNA_input + self.mRNA_input[:200]
            else:
                self.mRNA_input = self.mRNA_input + self.mRNA_input[:mRNA_length]

        parallelizer_arguments = [[], [], []]
        for _, (start_pos, codon) in enumerate(find_start_codons(self.mRNA_input[:mRNA_length], self.start_range)):
            parallelizer_arguments[0].append(self)
            parallelizer_arguments[1].append(start_pos)
            parallelizer_arguments[2].append(codon)
        if self.threads > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as multiprocessor:
                parallel_output = multiprocessor.map(self._parallel_dG, *parallelizer_arguments)
        else:
            parallel_output = map(self._parallel_dG, *parallelizer_arguments)

        parallel_output = [x for x in parallel_output if x is not None]
        self.results = parallel_output

        self.run = 1


    @staticmethod
    def _parallel_dG(factory_obj, start_pos, codon):
        try:

            # print "Top of calc_dG here"

            # Set dangles based on length between 5' end of mRNA and start codon
            if factory_obj.auto_dangles:
                cutoff = int(factory_obj.cutoff)
                if start_pos > cutoff:
                    dangles = "none"

                else:
                    dangles = "all"

            else:
                dangles = factory_obj.dangles_default
                # print "Auto Dangles set to ", self.dangles

            # Start codon energy
            dG_start_codon = factory_obj.start_codon_energies[codon]

            # Energy of mRNA folding. This also gets us the kinetic score
            dG_mRNA, mRNA_structure, kinetic_score, min_bp_prob = calc_dG_mRNA(factory_obj.mRNA_input, start_pos, dangles, factory_obj.constraints)

            # Energy of mRNA:rRNA hybridization & folding
            dg_mRNA_rRNA_output = calc_dG_mRNA_rRNA(factory_obj.mRNA_input, factory_obj.rRNA, start_pos, dangles, factory_obj.constraints)
            dG_mRNA_rRNA_withspacing = dg_mRNA_rRNA_output[0]
            mRNA_rRNA_structure = dg_mRNA_rRNA_output[1]
            spacing_value = dg_mRNA_rRNA_output[2]
            if not dG_mRNA_rRNA_withspacing:
                return None

            dG_mRNA_rRNA_withspacing -= 2.481  # Modifying hybridization penalty to match NUPACK

            dG_mRNA_rRNA_nospacing = mRNA_rRNA_structure["dG_mRNA_rRNA"]
            dG_mRNA_rRNA_nospacing -= 2.481  # Modifying hybridization penalty to match NUPACK

            # Standby site correction:
            dG_standby_site, corrected_structure = calc_dG_standby_site(mRNA_rRNA_structure,
                                                                                    dangles,
                                                                                    factory_obj.standby_site_length,
                                                                                    factory_obj.constraints,
                                                                                    factory_obj.rRNA)

            # Total energy is mRNA:rRNA + start - rRNA - mRNA - standby_site
            dG_total = dG_mRNA_rRNA_withspacing + dG_start_codon - dG_mRNA - dG_standby_site

            # Calculate dG to open SD sequence
            # ddG_SD_open = self.calc_dG_SDopen(mRNA_structure, mRNA_rRNA_structure)

            loop_result = calc_longest_loop_bulge(structure=mRNA_structure)
            helical_loop_list, bulge_loop_list = loop_result[0], loop_result[1]

            parallel_output = OSTIRResult(
                name= factory_obj.name,
                start_codon=codon,
                start_position=start_pos+1,
                expression=round(calc_expression_level(dG_total),
                                 factory_obj.decimal_places),
                RBS_distance_bp=int(spacing_value),
                dG_total=round(float(dG_total), factory_obj.decimal_places),
                dG_rRNA_mRNA=round(float(dG_mRNA_rRNA_nospacing), factory_obj.decimal_places),
                dG_mRNA=round(float(dG_mRNA), factory_obj.decimal_places),
                dG_spacing=round(float(mRNA_rRNA_structure["dG_spacing"]), factory_obj.decimal_places),
                dG_standby=round(float(dG_standby_site), factory_obj.decimal_places),
                dG_start_codon=round(float(dG_start_codon), factory_obj.decimal_places),

                mRNA_structure=mRNA_structure,
                mRNA_rRNA_uncorrected_structure=mRNA_rRNA_structure,
                mRNA_rRNA_corrected_structure=corrected_structure,
                longest_helix=calc_longest_helix(mRNA_structure),
                dG_start_energy=dG_start_codon,
                helical_loop=helical_loop_list,
                min_bp_prob=min_bp_prob,
                bulge_loop=bulge_loop_list,
                kinetic_score=kinetic_score
            )

        except ValueError as msg:
            if "leaderless start codon" not in str(msg):
                raise ValueError(msg) from msg
            print(msg)
            parallel_output = None

        return parallel_output


# ----------------------------------------------------------------------------------------------------------
# End RBS_Calculator class
# ----------------------------------------------------------------------------------------------------------
