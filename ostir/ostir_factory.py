#!/usr/bin/env python
# Given an mRNA sequence, this Python class predicts the dG_total and translation initiation rate.

# This file is part of the Ribosome Binding Site Calculator.

# The Ribosome Binding Site Calculator is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# The Ribosome Binding Site Calculator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.
# Copyright 2008-2009 is owned by the University of California Regents. All rights reserved.


try:
    from ostir.ViennaRNA import *
except ModuleNotFoundError:
    from ViennaRNA import *
import re
import math
import os
import concurrent.futures



class CalcError(Exception):
    """Base class for exceptions in this module."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class OSTIRFactory:

    # From OSTIR calibration using Salis2009 data. See calibration directory for procedure.
    Beta = 0.40002512
    RT_eff = 1/Beta
    logK = 7.279194329
    K = math.exp(logK)

    # Global parameters -- constants
    infinity = 1e12  # For all practical purposes, here.
    RNA_model = 'rna2004'
    start_codon_energies = {"ATG": -1.194, "AUG": -1.194, "GTG": -0.0748, "GUG": -0.0748, "TTG": -0.0435,
                            "UUG": -0.0435, "CTG": -0.03406, "CUG": -0.03406}  # hybridization to CAT
    auto_dangles = True
    dangles_default = "all"
    temp = 37.0
    optimal_spacing = 5  # aligned spacing

    # From OSTIR calibration using Salis2009 data. See calibration directory for procedure
    dG_spacing_constant_push = [17.20965071, 3.46341492, 1.790848365, 3.0]
    dG_spacing_constant_pull = [0.06422042, 0.275640836, 0.0]
    cutoff = 35  # number of nt +- start codon considering for folding
    standby_site_length = 4  # Number of nt before SD sequence that must be unpaired for ribosome binding
    energy_cutoff = 3.0
    start_codons = ["ATG", "AUG", "GTG", "GUG", "TTG", "UUG"]  # substituted U for T in actual calcs. Ignores CTG/CUG
    # rRNA = "acctcctta" #These are the last 9 nt (3' end) of the 16S rRNA in E. coli

    footprint = 1000  # Footprint of the 30S complex that prevents formation of secondary structures downstream of the start codon. Here, we assume that the entire post-start RNA sequence does not form secondary structures once the 30S complex has bound.

    def __init__(self, mRNA, start_range_1, rRNA, verbose=False):
        """
        Initializes the RBS Calculator class with the mRNA sequence and the range of start codon positions considered.
        start_range_1 is a pair of 1-indexed positions
        """

        # NuPACK.__init__(self,sequences,self.RNA_model)
        exp = re.compile('[ATGCU._]', re.IGNORECASE)
        if exp.match(mRNA) == None:
            raise ValueError(f"Invalid letters found in sequence {mRNA}. Only ATGCU accepted.")
        mRNA = mRNA.replace('.', '')
        mRNA = mRNA.replace('_', '')

        if start_range_1[0] < 1: start_range_1[0] = 1
        if start_range_1[1] > len(mRNA): start_range_1[1] = len(mRNA)


        self.install_location = os.path.dirname(os.path.realpath(__file__))
        self.name = mRNA
        self.mRNA_input = mRNA.upper()
        self.rRNA = rRNA
        self.rRNA_len = len(self.rRNA)
        self.mRNA_len = len(self.mRNA_input)
        self.total_sequence_length = len(mRNA) + len(self.rRNA)
        self.dG_rRNA = self.calc_dG_rRNA()
        self.run = 0
        self.start_range_1 = start_range_1
        self.verbose = verbose
        self.threads = 1

    def find_min(self, input_list):
        """Finds the minimum of a list of numbers."""

        min_item = self.infinity
        min_index = 0

        for i, item in enumerate(input_list):
            if item < min_item:
                min_item = item
                min_index = i

        return (min_item, min_index)

    def find_start_codons(self, sequence):
        """Finds all start codons in an mRNA sequence. Creates a list."""

        self.start_position_list = []
        self.start_codon_list = []

        seq_len = len(sequence)

        #Switch to zero-indexed positions happens here
        end_0 = min(self.start_range_1[1]-1, seq_len - 2)
        begin_0 = min(self.start_range_1[0]-1, end_0)

        for i in range(begin_0, end_0 + 1):
            codon = sequence[i:i + 3]
            if codon.upper() in self.start_codons:
                self.start_position_list.append(i)
                self.start_codon_list.append(codon)
                yield (i, codon)
            else:
                pass

    def calc_aligned_spacing(self, mRNA, start_pos, bp_x, bp_y):
        """Calculates the aligned spacing between the 16S rRNA binding site and the start codon."""

        # rRNA is the concatenated at the end of the sequence in 5' to 3' direction
        # first: identify the farthest 3' nt in the rRNA that binds to the mRNA and return its mRNA base pairer

        Ok = False
        seq_len = len(mRNA) + self.rRNA_len
        for (rRNA_nt) in range(seq_len, seq_len - self.rRNA_len, -1):

            if rRNA_nt in bp_y:
                rRNA_pos = bp_y.index(rRNA_nt)
                if bp_x[rRNA_pos] < start_pos:
                    Ok = True
                    farthest_3_prime_rRNA = rRNA_nt - len(mRNA)

                    mRNA_nt = bp_x[rRNA_pos]
                    distance_to_start = start_pos - mRNA_nt + 1  # start_pos is counting starting from 0 (python)

                    break
                else:
                    break

        if Ok:
            aligned_spacing = distance_to_start - farthest_3_prime_rRNA
        else:
            aligned_spacing = self.infinity

        return aligned_spacing

    def calc_dG_spacing(self, aligned_spacing):
        """Calculates the dG_spacing according to the value of the aligned spacing. This relationship was determined through experiments."""

        if (aligned_spacing < self.optimal_spacing):
            ds = aligned_spacing - self.optimal_spacing

            dG_spacing_penalty = self.dG_spacing_constant_push[0] / (
                        1.0 + math.exp(self.dG_spacing_constant_push[1] * (ds + self.dG_spacing_constant_push[2]))) ** \
                                 self.dG_spacing_constant_push[3]

        else:
            ds = aligned_spacing - self.optimal_spacing
            dG_spacing_penalty = self.dG_spacing_constant_pull[0] * ds * ds + self.dG_spacing_constant_pull[1] * ds + \
                                 self.dG_spacing_constant_pull[2]


        return dG_spacing_penalty

    def calc_dG_mRNA_rRNA(self, start_pos, dangles):
        """Calculates the dG_mRNA_rRNA from the mRNA and rRNA sequence. Considers all feasible 16S rRNA binding sites and includes the effects of non-optimal spacing."""

        begin = max(0, start_pos - self.cutoff)
        mRNA_len = min(len(self.mRNA_input), start_pos + self.cutoff)
        start_pos_in_subsequence = min(start_pos, self.cutoff)
        startpos_to_end_len = mRNA_len - start_pos_in_subsequence - begin

        # 1. identify a list of rRNA-binding sites. Binding sites are hybridizations between the mRNA and rRNA and can include mismatches, bulges, etc. Intra-molecular folding is also allowed within the mRNA. The subopt program is used to generate a list of optimal & suboptimal binding sites.
        # Constraints: the entire rRNA-binding site must be upstream of the start codon
        mRNA = self.mRNA_input[begin:start_pos]

        if len(mRNA) == 0:
            raise CalcError("Warning: There is a leaderless start codon, which is being ignored.")

        # print "After exception"

        fold = ViennaRNA([mRNA, self.rRNA], material=self.RNA_model)
        # print(fold)
        fold.subopt([1, 2], None, self.energy_cutoff, dangles=dangles, Temp=self.temp)
        if len(fold["subopt_basepairing_x"]) == 0:
            return None, None, None

        # 2. Calculate dG_spacing for each 16S rRNA binding site

        # Calculate the aligned spacing for each binding site in the list
        aligned_spacing = []
        proposed_pairings = list(zip(fold["subopt_basepairing_x"], fold["subopt_basepairing_y"]))
        for (bp_x, bp_y) in proposed_pairings:
            aligned_spacing.append(self.calc_aligned_spacing(mRNA, start_pos_in_subsequence, bp_x, bp_y))

        dG_spacing_list = []
        dG_mRNA_rRNA = []
        dG_mRNA_rRNA_withspacing = []

        # Calculate dG_spacing using aligned spacing value. Add it to dG_mRNA_rRNA.
        for (counter) in range(len(fold["subopt_basepairing_x"])):
            dG_mRNA_rRNA.append(fold["subopt_energy"][counter])
            val = self.calc_dG_spacing(aligned_spacing[counter])
            dG_spacing_list.append(val)
            dG_mRNA_rRNA_withspacing.append(val + fold["subopt_energy"][counter])


        # 3. Find 16S rRNA binding site that minimizes dG_spacing+dG_mRNA_rRNA.
        [dG_mRNA_rRNA_folding, index] = self.find_min(dG_mRNA_rRNA_withspacing)
        dG_spacing_final = dG_spacing_list[index]
        spacing_value = aligned_spacing[index]



        dG_mRNA_rRNA_nospacing = dG_mRNA_rRNA[index]

        # Check: Is the dG spacing large compared to the energy gap? If so, this means the list of suboptimal 16S rRNA binding sites generated by subopt is too short.
        if dG_spacing_final > self.energy_cutoff:
            if self.verbose:
                print("Warning: The spacing penalty is greater than the energy gap. dG (spacing) = ", dG_spacing_final)

        # 4. Identify the 5' and 3' ends of the identified 16S rRNA binding site. Create a base pair list.

        most_5p_mRNA = self.infinity
        most_3p_mRNA = -self.infinity

        # Generate a list of rRNA-mRNA base paired nucleotides
        bp_x_target = []
        bp_y_target = []

        bp_x = fold["subopt_basepairing_x"][index]
        bp_y = fold["subopt_basepairing_y"][index]


        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_y > len(mRNA):  # nt is rRNA
                most_5p_mRNA = min(most_5p_mRNA, bp_x[bp_y.index(nt_y)])
                most_3p_mRNA = max(most_3p_mRNA, bp_x[bp_y.index(nt_y)])
                bp_x_target.append(nt_x)
                bp_y_target.append(nt_y)


        # The rRNA-binding site is between the nucleotides at positions most_5p_mRNA and most_3p_mRNA
        # Now, fold the pre-sequence, rRNA-binding-sequence and post-sequence separately. Take their base pairings and combine them together. Calculate the total energy. For secondary structures, this splitting operation is allowed.

        # We postulate that not all of the post-sequence can form secondary structures. Once the 30S complex binds to the mRNA, it prevents the formation of secondary structures that are mutually exclusive with ribosome binding. We define self.footprint to be the length of the 30S complex footprint. Here, we assume that the entire mRNA sequence downstream of the 16S rRNA binding site can not form secondary structures.
        begin = int(begin)
        most_5p_mRNA = int(most_5p_mRNA)
        most_3p_mRNA = int(most_3p_mRNA)

        mRNA_pre = self.mRNA_input[begin:begin + most_5p_mRNA - 1]
        post_window_end = mRNA_len + 1
        post_window_begin = min(start_pos + self.footprint, post_window_end)  # Footprint
        post_window_end = mRNA_len + 1
        mRNA_post = self.mRNA_input[post_window_begin:post_window_end]

        mRNA_pre_len = len(mRNA_pre)
        mRNA_post_len = len(mRNA_post)
        mRNA_rRNA_binding_len = most_3p_mRNA - most_5p_mRNA + 1
        total_folded_len = mRNA_pre_len + mRNA_post_len + mRNA_rRNA_binding_len

        total_bp_x = []
        total_bp_y = []

        # Calculate pre-sequence folding
        if len(mRNA_pre) > 0:
            fold_pre = ViennaRNA([mRNA_pre], material=self.RNA_model)
            fold_pre.mfe([1], None, Temp=self.temp, dangles=dangles, )
            bp_x_pre = fold_pre["mfe_basepairing_x"][0]
            bp_y_pre = fold_pre["mfe_basepairing_y"][0]

        else:
            bp_x_pre = []
            bp_y_pre = []

        # Add pre-sequence base pairings to total base pairings
        offset = 0  # Begins at 0
        for (nt_x, nt_y) in zip(bp_x_pre, bp_y_pre):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + offset)


        # Add rRNA-binding site base pairings to total base pairings
        offset = 0  # Begins at zero
        if startpos_to_end_len < self.cutoff:
            rRNA_offset = startpos_to_end_len
        else:
            rRNA_offset = startpos_to_end_len

        for (nt_x, nt_y) in zip(bp_x_target, bp_y_target):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + rRNA_offset)

        # Calculate post-sequence folding
        if len(mRNA_post) > 0:
            fold_post = ViennaRNA([mRNA_post], material=self.RNA_model)
            fold_post.mfe([1], None, Temp=self.temp, dangles=dangles)
            bp_x_post = fold_post["mfe_basepairing_x"][0]
            bp_y_post = fold_post["mfe_basepairing_y"][0]
        else:
            bp_x_post = []
            bp_y_post = []

        offset = post_window_begin - begin
        for (nt_x, nt_y) in zip(bp_x_post, bp_y_post):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + offset)


        mRNA = self.mRNA_input[begin:mRNA_len]
        fold = ViennaRNA([mRNA, self.rRNA], material=self.RNA_model)

        total_energy = fold.energy([1, 2], total_bp_x, total_bp_y, Temp=self.temp, dangles=dangles)

        energy_nowindows = dG_mRNA_rRNA_nospacing
        total_energy_withspacing = total_energy + dG_spacing_final

        structure = fold
        structure["program"] = "subopt"
        structure["mRNA"] = mRNA
        structure["MinStructureID"] = 0
        structure["dG_mRNA_rRNA"] = total_energy
        structure["dG_mRNA_rRNA_withspacing"] = total_energy_withspacing
        structure["dG_spacing"] = dG_spacing_final
        structure["subopt_energy"] = [total_energy_withspacing]
        structure["subopt_basepairing_x"] = [total_bp_x]
        structure["subopt_basepairing_y"] = [total_bp_y]
        structure["subopt_composition"] = [1, 2]
        structure["bp_x"] = total_bp_x
        structure["bp_y"] = total_bp_y

        return total_energy_withspacing, structure, spacing_value

    def calc_dG_standby_site(self, structure_old, dangles, rRNA_binding=True):
        """Calculates the dG_standby given the structure of the mRNA:rRNA complex"""

        # To calculate the mfe structure while disallowing base pairing at the standby site, we split the folded mRNA sequence into three parts: (i) a pre-sequence (before the standby site) that can fold; (ii) the standby site, which can not fold; (iii) the 16S rRNA binding site and downstream sequence, which has been previously folded.

        import copy
        structure = copy.deepcopy(structure_old)
        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]
        energy_before = structure["dG_mRNA_rRNA"]  # without spacing effects

        # Identify the most 5p mRNA nt that is bound to rRNA
        most_5p_mRNA = 0
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y > len(mRNA):  # nt_x is mRNA, nt_y is rRNA, they are bound.
                most_5p_mRNA = nt_x  # starts counting from 0
                break

        # Extract the base pairings that are 3' of the most_5p_mRNA base pairing
        bp_x_3p = []
        bp_y_3p = []
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x >= most_5p_mRNA:
                bp_x_3p.append(nt_x)
                bp_y_3p.append(nt_y)

        # Create the mRNA subsequence
        mRNA_subsequence = mRNA[0:max(0, most_5p_mRNA - self.standby_site_length - 1)]
        standby_site = mRNA[most_5p_mRNA - self.standby_site_length - 1:most_5p_mRNA]

        # Fold it and extract the base pairings
        if (len(mRNA_subsequence)) > 0:
            fold = ViennaRNA([mRNA_subsequence], material=self.RNA_model)
            fold.mfe([1], None, Temp=self.temp, dangles=dangles)
            energy_after_5p = fold["mfe_energy"][0]
            bp_x_5p = fold["mfe_basepairing_x"][0]  # [0] added 12/13/07
            bp_y_5p = fold["mfe_basepairing_y"][0]
        else:
            bp_x_5p = []
            bp_y_5p = []
            energy_after_5p = 0.0

        # Put the sets of base pairings together
        bp_x_after = []
        bp_y_after = []
        for (nt_x, nt_y) in zip(bp_x_5p, bp_y_5p):
            bp_x_after.append(nt_x)
            bp_y_after.append(nt_y)

        for (nt_x, nt_y) in zip(bp_x_3p, bp_y_3p):
            bp_x_after.append(nt_x)
            bp_y_after.append(nt_y)

        # Calculate its energy
        fold = ViennaRNA([mRNA, self.rRNA], material=self.RNA_model)
        energy_after = fold.energy([1, 2], bp_x_after, bp_y_after, dangles=dangles, Temp=self.temp)

        dG_standby_site = energy_before - energy_after

        if (dG_standby_site > 0.0): dG_standby_site = 0.0

        index = structure["MinStructureID"]
        structure["bp_x"] = bp_x_after
        structure["bp_y"] = bp_y_after
        structure["subopt_basepairing_x"][index] = bp_x_after
        structure["subopt_basepairing_y"][index] = bp_y_after
        structure["subopt_energy"][index] = energy_after
        structure["dG_mRNA_rRNA_corrected"] = energy_after

        return (dG_standby_site, structure)

    def calc_dG_mRNA(self, start_pos, dangles):
        """Calculates the dG_mRNA given the mRNA sequence."""
        mRNA = self.mRNA_input[max(0, start_pos - self.cutoff):min(len(self.mRNA_input), start_pos + self.cutoff)]
        fold = ViennaRNA([mRNA], self.RNA_model)
        fold.mfe([1], None, Temp=self.temp, dangles=dangles)

        structure = fold
        structure["mRNA"] = mRNA
        structure["bp_x"] = fold["mfe_basepairing_x"][0]
        structure["bp_y"] = fold["mfe_basepairing_y"][0]
        structure["dG_mRNA"] = fold["mfe_energy"][0]
        structure["MinStructureID"] = 0

        dG_mRNA_folding = fold["mfe_energy"][0]

        return (dG_mRNA_folding, structure)

    def calc_dG_rRNA(self):
        """Calculates the dG of folding for the last 9 nt of the 16S rRNA. Not used in the free energy model."""
        fold = ViennaRNA([self.rRNA], self.RNA_model)
        fold.mfe([1], constraints=None, Temp=self.temp, dangles="all")
        dG_rRNA_folding = fold["mfe_energy"][0]
        return dG_rRNA_folding

    def calc_kinetic_score(self, structure=None, mRNA_in=None, bp_x_in=None, bp_y_in=None):
        """Calculate a "kinetic score", a heuristic measure of the maximum time required for the mRNA secondary structure to form. This is related to the RNA polymer model by David et. al. This heuristic should not be used in any way to quantify the folding kinetics of an mRNA sequence because it completely ignores cooperative RNA folding mechanisms, such as zipping or strand displacement. Here, we use it to eliminate mRNA sequences that MAY fold slowly."""

        if not (structure is None):
            program = structure["program"]
            mRNA = structure["mRNA"]
            index = structure["MinStructureID"]
            bp_x = structure[program + "_basepairing_x"][index]
            bp_y = structure[program + "_basepairing_y"][index]

        if not (bp_x_in is None) and not (bp_y_in is None) and not (mRNA_in is None):
            mRNA = mRNA_in[:]
            bp_x = bp_x_in[:]
            bp_y = bp_y_in[:]

        largest_range_helix = 0
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y <= len(mRNA):
                val = nt_y - nt_x
                largest_range_helix = max(val, largest_range_helix)

        kinetic_score = float(largest_range_helix) / float(len(mRNA))
        if float(largest_range_helix) > 0:
            min_bp_prob = float(largest_range_helix) ** (-1.44)  # From David et. al.
        else:
            min_bp_prob = 1.0

        return (kinetic_score, min_bp_prob)

    def calc_most_5p_mRNA(self, structure_old):
        """Calculates the most 5' nucleotide in the 16S rRNA binding site."""

        import copy
        structure = copy.deepcopy(structure_old)
        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]

        # Identify the most 5p mRNA nt that is bound to rRNA
        most_5p_mRNA = 0
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y > len(mRNA):  # nt_x is mRNA, nt_y is rRNA, they are bound.
                most_5p_mRNA = nt_x
                break

        return most_5p_mRNA

    def calc_longest_helix(self, structure):
        """Calculate the longest helical structure (longest contiguous list of base pairings) in the secondary structure"""

        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]

        longest_helix = 0
        helix_length = 1

        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if (bp_x.count(nt_x + 1) > 0 and bp_y.count(nt_y - 1) > 0):
                helix_length += 1
            else:
                longest_helix = max(longest_helix, helix_length)
                helix_length = 1

        return longest_helix

    def calc_longest_loop_bulge(self, structure, output_start_end=False, InRBSOnly=False, RBS=None):
        """Calculate the longest helical loop and bulge structure (longest contiguous list of un-base paired nucleotides starting and ending with a helix (loop -> same helix, bulge -> different helix) in the secondary structure"""

        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]

        loop_length = 0
        begin_helix = 1

        bulge_loop_list = []
        helical_loop_list = []

        if output_start_end:
            bulge_loop_start_end = []
            helical_loop_start_end = []

        if InRBSOnly and RBS is not None:
            RBS_begin = mRNA.find(RBS)
            RBS_end = RBS_begin + len(RBS)
            nucleotide_range = list(range(RBS_begin, RBS_end + 1))

        else:
            nucleotide_range = list(range(1, len(mRNA) + 1))

        # Find loops. Find bulges.
        for n in nucleotide_range:
            if bp_x.count(n) == 0 and bp_y.count(n) == 0:  # nth nucleotide is not base-paired.

                # Determine if nearest neighbor nucleotides are base-paired
                (x1, x2, y1, y2) = (bp_x.count(n - 1), bp_x.count(n + 1), bp_y.count(n - 1), bp_y.count(n + 1))

                # print "#", n, (x1,x2,y1,y2)

                # middle unpaired nt
                if (x1, x2, y1, y2) == (0, 0, 0, 0):
                    loop_length += 1

                # single mismatch -- loop
                elif (x1, x2, y1, y2) == (1, 0, 0, 1) or (x1, x2, y1, y2) == (0, 1, 1, 0):
                    loop_length += 1
                    begin_helix = n - 1
                    end_helix = n + 1

                # single mismatch -- bulge
                elif (x1, x2, y1, y2) == (1, 1, 0, 0) or (x1, x2, y1, y2) == (0, 0, 1, 1):
                    loop_length += 1
                    begin_helix = n - 1
                    end_helix = n + 1

                # starting unpaired nt
                elif (x1, x2, y1, y2) == (1, 0, 0, 0) or (x1, x2, y1, y2) == (0, 0, 1, 0):
                    loop_length += 1
                    begin_helix = n - 1

                # ending unpaired nt
                elif (x1, x2, y1, y2) == (0, 1, 0, 0) or (x1, x2, y1, y2) == (0, 0, 0, 1):
                    loop_length += 1
                    end_helix = n + 1

                # 1,0,1,0 is impossible w/o psuedoknots
                # 0,1,0,1 is impossible w/o psuedoknots
                # Also, all binary combinations with 3 or 4 true are impossible (n-1 or n+1 can not be in both bp_x and bp_y)


            elif loop_length > 0:
                # Bulge or loop?
                # loop

                # print "begin = ", begin_helix
                # print "end = ", end_helix

                if (bp_x.count(begin_helix) > 0 and bp_y.count(end_helix) > 0 and bp_x.index(begin_helix) == bp_y.index(
                        end_helix)):
                    helical_loop_list.append(loop_length)
                    loop_length = 0

                    if output_start_end:
                        # Also return the starting and ending positions of each loop/bulge
                        helical_loop_start_end.append((begin_helix, end_helix))

                else:

                    bp_end = 0
                    bp_begin = 0

                    if (bp_x.count(end_helix) > 0): bp_begin = bp_y[bp_x.index(end_helix)]
                    if (bp_y.count(end_helix) > 0): bp_end = bp_x[bp_y.index(end_helix)]

                    if (bp_x.count(begin_helix) > 0): bp_end = bp_y[bp_x.index(begin_helix)]
                    if (bp_y.count(begin_helix) > 0): bp_begin = bp_x[bp_y.index(begin_helix)]

                    if bp_end > bp_begin:
                        bulge_loop_list.append(loop_length)
                        loop_length = 0

                        if output_start_end:
                            # Also return the starting and ending positions of each loop/bulge
                            bulge_loop_start_end.append((begin_helix, end_helix))

                    else:
                        loop_length = 0
        if output_start_end:
            return (helical_loop_list, bulge_loop_list, helical_loop_start_end, bulge_loop_start_end)
        else:
            return (helical_loop_list, bulge_loop_list)

    def calc_dG(self):
        """Calculates each dG term in the free energy model and sums them together to create dG_total"""

        # Initialization of data structures
        self.start_pos_list = []
        self.dG_total_list = []
        self.dG_details = []
        self.dG_mRNA_list = []
        self.dG_mRNA_rRNA_list = []
        self.fold_x_list = []
        self.fold_y_list = []
        self.dG_start_energy_list = []
        self.dG_spacing_list = []
        self.mRNA_structure_list = []
        self.mRNA_rRNA_uncorrected_structure_list = []
        self.mRNA_rRNA_corrected_structure_list = []
        self.dG_standby_site_list = []
        self.dG_spacer_site_list = []
        self.kinetic_score_list = []
        self.min_bp_prob_list = []
        self.longest_helix_list = []
        self.three_state_indicator_list = []
        self.helical_loop_list_list = []
        self.bulge_loop_list_list = []
        self.start_codon_list = []

        self.dS1_list = []
        self.dS2_list = []
        self.most_5p_mRNA_list = []
        self.Expression_list = []
        self.dg_scores = []

        parallelizer_arguments = [[], [], []]
        for _, (start_pos, codon) in enumerate(self.find_start_codons(self.mRNA_input)):
            parallelizer_arguments[0].append(self)
            parallelizer_arguments[1].append(start_pos)
            parallelizer_arguments[2].append(codon)
            checked_codon = self.mRNA_input[start_pos:start_pos+3]
            if checked_codon != codon:
                raise ValueError(f'Codon mismatch at {start_pos}: {codon, checked_codon}')

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as multiprocessor:
            parallel_output = multiprocessor.map(self._parallel_dG, *parallelizer_arguments)

        parallel_output = [x for x in parallel_output if x is not None]
        for output in parallel_output:
            if not output['start_position_list'] or output['start_position_list'] == None:
                continue
            self.dg_scores.append(output['dG_scores'])
            self.mRNA_structure_list.append(output['mRNA_structure_list'])
            self.mRNA_rRNA_uncorrected_structure_list.append(output['mRNA_rRNA_uncorrected_structure_list'])
            self.mRNA_rRNA_corrected_structure_list.append(output['mRNA_rRNA_corrected_structure_list'])
            self.most_5p_mRNA_list.append(output['most_5p_mRNA_list'])
            self.longest_helix_list.append(output['longest_helix_list'])
            self.dG_start_energy_list.append(output['dG_start_energy_list'])
            self.dG_spacing_list.append(output['dG_spacing_list'])
            self.dG_total_list.append(output['dG_total_list'])
            self.dG_standby_site_list.append(output['dG_standby_site_list'])
            self.helical_loop_list_list.append(output['helical_loop_list_list'])
            self.min_bp_prob_list.append(output['min_bp_prob_list'])
            self.bulge_loop_list_list.append(output['bulge_loop_list_list'])
            self.kinetic_score_list.append(output['kinetic_score_list'])
            if output['start_position_list']:
                self.start_pos_list.append(output['start_position_list'])
            if output['Expression_list']:
                self.Expression_list.append(output['Expression_list'])
            self.dG_details.append(output['dG_details'])
            self.start_codon_list.append(output['codon'])

        self.run = 1

    @staticmethod
    def _parallel_dG(ostir_factory_object, start_pos, codon):
        try:

            # print "Top of calc_dG here"

            # Set dangles based on length between 5' end of mRNA and start codon
            if ostir_factory_object.auto_dangles:
                cutoff = int(ostir_factory_object.cutoff)
                if start_pos > cutoff:
                    dangles = "none"

                else:
                    dangles = "all"

            else:
                dangles = ostir_factory_object.dangles_default
                # print "Auto Dangles set to ", self.dangles

            # Start codon energy
            dG_start_codon = ostir_factory_object.start_codon_energies[codon]

            # Energy of mRNA folding
            [dG_mRNA, mRNA_structure] = ostir_factory_object.calc_dG_mRNA(start_pos, dangles)

            # Energy of mRNA:rRNA hybridization & folding
            [dG_mRNA_rRNA_withspacing, mRNA_rRNA_structure, spacing_value] = ostir_factory_object.calc_dG_mRNA_rRNA(start_pos, dangles)
            if not dG_mRNA_rRNA_withspacing:
                return None

            dG_mRNA_rRNA_withspacing -= 2.481  # Modifying hybridization penalty to match NUPACK

            dG_mRNA_rRNA_nospacing = mRNA_rRNA_structure["dG_mRNA_rRNA"]
            dG_mRNA_rRNA_nospacing -= 2.481  # Modifying hybridization penalty to match NUPACK


            # Standby site correction:
            [dG_standby_site, corrected_structure] = ostir_factory_object.calc_dG_standby_site(mRNA_rRNA_structure,
                                                                                               rRNA_binding=True, dangles=dangles)

            # Total energy is mRNA:rRNA + start - rRNA - mRNA - standby_site
            dG_total = dG_mRNA_rRNA_withspacing + dG_start_codon - dG_mRNA - dG_standby_site

            #print(f"{dG_total} = {dG_mRNA_rRNA_withspacing} + {dG_start_codon} - {dG_mRNA} - {dG_standby_site} - 2.481")

            # Calculate 'kinetic score': directly related to probability of base pair formation
            (kinetic_score, min_bp_prob) = ostir_factory_object.calc_kinetic_score(mRNA_structure)

            # Calculate dG to open SD sequence
            # ddG_SD_open = self.calc_dG_SDopen(mRNA_structure, mRNA_rRNA_structure)


            (helical_loop_list, bulge_loop_list) = ostir_factory_object.calc_longest_loop_bulge(structure=mRNA_structure)

            parallel_output = {
                'dG_scores': [dG_mRNA_rRNA_withspacing, dG_start_codon, dG_mRNA, dG_standby_site],
                'mRNA_structure_list': mRNA_structure,
                'mRNA_rRNA_uncorrected_structure_list': mRNA_rRNA_structure,
                'mRNA_rRNA_corrected_structure_list': corrected_structure,
                'most_5p_mRNA_list': ostir_factory_object.calc_most_5p_mRNA(mRNA_rRNA_structure),
                'longest_helix_list': ostir_factory_object.calc_longest_helix(structure=mRNA_structure),
                'dG_start_energy_list': dG_start_codon,
                'dG_mRNA_list': dG_mRNA,
                'dG_mRNA_rRNA_list': dG_mRNA_rRNA_nospacing,
                'dG_spacing_list': mRNA_rRNA_structure["dG_spacing"],
                'dG_standby_site_list': dG_standby_site,
                'dG_total_list': dG_total,
                'helical_loop_list_list': helical_loop_list,
                'bulge_loop_list_list': bulge_loop_list,
                'min_bp_prob_list': min_bp_prob,
                'kinetic_score_list': kinetic_score,
                'start_position_list': start_pos,
                'Expression_list': ostir_factory_object.calc_expression_level(dG_total),
                'dG_details': [dG_mRNA_rRNA_nospacing, dG_start_codon, dG_mRNA, dG_standby_site,
                               mRNA_rRNA_structure["dG_spacing"], spacing_value],
                'codon': codon
            }

        except CalcError as msg:
            print(msg)
            parallel_output = {
                'dG_scores': [],
                'mRNA_structure_list': [],
                'mRNA_rRNA_uncorrected_structure_list': [],
                'mRNA_rRNA_corrected_structure_list': [],
                'most_5p_mRNA_list': ostir_factory_object.infinity,
                'longest_helix_list': ostir_factory_object.infinity,
                'dG_start_energy_list': ostir_factory_object.infinity,
                'dG_mRNA_list': ostir_factory_object.infinity,
                'dG_mRNA_rRNA_list': ostir_factory_object.infinity,
                'dG_spacing_list': ostir_factory_object.infinity,
                'dG_standby_site_list': ostir_factory_object.infinity,
                'dG_total_list': ostir_factory_object.infinity,
                'helical_loop_list_list': ostir_factory_object.infinity,
                'bulge_loop_list_list': ostir_factory_object.infinity,
                'min_bp_prob_list': ostir_factory_object.infinity,
                'kinetic_score_list': ostir_factory_object.infinity,
                'start_position_list': None,
                'Expression_list': None,
                'dG_details': None,
                'codon': None
            }
        return parallel_output


    def calc_expression_level(self, dG):
        import math
        return self.K * math.exp(-dG / self.RT_eff)

# ----------------------------------------------------------------------------------------------------------
# End RBS_Calculator class
# ----------------------------------------------------------------------------------------------------------
