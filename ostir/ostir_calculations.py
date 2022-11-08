# cython: profile=True

import math
from dataclasses import dataclass
from copy import deepcopy
from .ViennaRNA import ViennaRNA, subopt, mfe, energy
import numpy as np
from datetime import datetime

@dataclass
class OstirConstants():
    Beta: float = 0.40002512
    RT_eff: float = 1/Beta
    logK: float = 7.279194329
    K: float = math.exp(logK)
    infinity = np.inf
    RNA_model: str = "rna2004"
    auto_dangles: bool = True
    dangles_default: str = "all"
    temp: float = 37.0
    optimal_spacing: int = 5  # aligned spacing
    cutoff: int = 35

        # From OSTIR calibration using Salis2009 data. See calibration directory for procedure
    dG_spacing_constant_push = np.array([17.20965071, 3.46341492, 1.790848365, 3.0], dtype=np.float64)
    dG_spacing_constant_pull = np.array([0.06422042, 0.275640836, 0.0], dtype=np.float64)
    cutoff: int = 35  # number of nt +- start codon considering for folding
    standby_site_length: int = 4  # Number of nt before SD sequence that must be unpaired for ribosome binding
    start_codons: list = ("ATG", # substituted U for T in actual calcs. Ignores CTG/CUG
                             "AUG",
                             "GTG",
                             "GUG",
                             "TTG",
                             "UUG")
    footprint: int = 1000
    energy_cutoff: float = 3.0
    verbose: bool = False

@dataclass
class StartEnergies:
    """Class to hold the start codon energies"""
    ATG: float = -1.194
    AUG: float = -1.194
    GTG: float = -0.0748
    GUG: float = -0.0748
    TTG: float = -0.0435
    UUG: float = -0.0435
    CTG: float = -0.03406
    CUG: float = -0.03406

ostir_constants = OstirConstants()
start_energies = StartEnergies()

# Timer decorator
def timer(func):
    def wrapper(*args, **kwargs):
        start = datetime.now()
        result = func(*args, **kwargs)
        end = datetime.now()
        print(f"[{datetime.now()}] {func.__name__} took {end - start}")
        return result
    return wrapper


def calc_longest_loop_bulge(structure, output_start_end=False, InRBSOnly=False, RBS=None):
    """Calculate the longest helical loop and bulge structure
    (longest contiguous list of un-base paired nucleotides starting and
    ending with a helix (loop -> same helix, bulge -> different helix) in the secondary structure"""

    mRNA = structure["mRNA"]

    bp_x = structure["bp_x"]
    bp_y = structure["bp_y"]

    loop_length = 0
    begin_helix = 1

    bulge_loop_list = []
    helical_loop_list = []
    end_helix = 0

    if output_start_end:
        bulge_loop_start_end = []
        helical_loop_start_end = []

    if InRBSOnly and RBS is not None:
        RBS_begin = mRNA.find(RBS)
        RBS_end = RBS_begin + len(RBS)
        loop_start = RBS_begin
        loop_end = RBS_end +1

    else:
        loop_start = 1
        loop_end = len(mRNA) + 1


    # Find loops. Find bulges.
    for i in range(loop_start, loop_end):
        if np.count_nonzero(bp_x == i) == 0 and np.count_nonzero(bp_y == i) == 0:  # nth nucleotide is not base-paired.
            # Determine if nearest neighbor nucleotides are base-paired
            x_1, x_2, y_1, y_2 = np.count_nonzero(bp_x == i - 1), np.count_nonzero(bp_x == i + 1), np.count_nonzero(bp_y == i - 1), np.count_nonzero(bp_y == i + 1)

            # print "#", n, (x_1,x_2,y_1,y_2)

            # middle unpaired nt
            if (x_1, x_2, y_1, y_2) == (0, 0, 0, 0):
                loop_length += 1

            # single mismatch -- loop
            elif (x_1, x_2, y_1, y_2) == (1, 0, 0, 1) or (x_1, x_2, y_1, y_2) == (0, 1, 1, 0):
                loop_length += 1
                begin_helix = i - 1
                end_helix = i + 1

            # single mismatch -- bulge
            elif (x_1, x_2, y_1, y_2) == (1, 1, 0, 0) or (x_1, x_2, y_1, y_2) == (0, 0, 1, 1):
                loop_length += 1
                begin_helix = i - 1
                end_helix = i + 1

            # starting unpaired nt
            elif (x_1, x_2, y_1, y_2) == (1, 0, 0, 0) or (x_1, x_2, y_1, y_2) == (0, 0, 1, 0):
                loop_length += 1
                begin_helix = i - 1

            # ending unpaired nt
            elif (x_1, x_2, y_1, y_2) == (0, 1, 0, 0) or (x_1, x_2, y_1, y_2) == (0, 0, 0, 1):
                loop_length += 1
                end_helix = i + 1

            # 1,0,1,0 is impossible w/o psuedoknots
            # 0,1,0,1 is impossible w/o psuedoknots
            # Also, all binary combinations with 3 or 4 true are impossible
            # (n-1 or n+1 can not be in both bp_x and bp_y)


        elif loop_length > 0:
            # Bulge or loop?
            # loop

            # print "begin = ", begin_helix
            # print "end = ", end_helix

            if (np.count_nonzero(bp_x == begin_helix) > 0 and np.count_nonzero(bp_y == end_helix) > 0 and np.where(bp_x == begin_helix) == np.where(bp_y ==
                    end_helix)):
                helical_loop_list.append(loop_length)
                loop_length = 0

                if output_start_end:
                    # Also return the starting and ending positions of each loop/bulge
                    helical_loop_start_end.append((begin_helix, end_helix))

            else:

                bp_end = 0
                bp_begin = 0

                if np.count_nonzero(bp_x == end_helix) > 0:
                    bp_begin = bp_y[np.where(bp_x == end_helix)]
                if np.count_nonzero(bp_y == end_helix) > 0:
                    bp_end = bp_x[np.where(bp_y == end_helix)]

                if np.count_nonzero(bp_x == begin_helix) > 0:
                    bp_end = bp_y[np.where(bp_x == begin_helix)]
                if np.count_nonzero(bp_y == begin_helix) > 0:
                    bp_begin = bp_x[np.where(bp_y == begin_helix)]

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


def calc_longest_helix(structure: ViennaRNA):
    """Calculate the longest helical structure (longest contiguous list of base pairings)
    in the secondary structure"""

    bp_x = structure["bp_x"]
    bp_y = structure["bp_y"]

    longest_helix = 0
    helix_length = 1

    for (nt_x, nt_y) in zip(bp_x, bp_y):
        if (np.count_nonzero(bp_x == nt_x + 1) > 0 and np.count_nonzero(bp_y == nt_y - 1) > 0):
            helix_length += 1
        else:
            longest_helix = max(longest_helix, helix_length)
            helix_length = 1

    return longest_helix


def calc_kinetic_score(mRNA_in=None, bp_x_in=None, bp_y_in=None):
    """Calculate a "kinetic score", a heuristic measure of the maximum time required for
    the mRNA secondary structure to form. This is related to the RNA polymer model by David et. al.
    This heuristic should not be used in any way to quantify the folding kinetics of an mRNA sequence
    because it completely ignores cooperative RNA folding mechanisms, such as zipping or strand
    displacement. Here, we use it to eliminate mRNA sequences that MAY fold slowly."""


    mRNA = mRNA_in
    bp_x = bp_x_in
    bp_y = bp_y_in

    mrnalen = len(mRNA)
    largest_range_helix = 0

    for (nt_x, nt_y) in zip(bp_x, bp_y):
        if nt_x <= len(mRNA) and nt_y <= mrnalen:
            val = nt_y - nt_x
            largest_range_helix = max(val, largest_range_helix)

    kinetic_score = largest_range_helix / mrnalen
    if largest_range_helix > 0:
        min_bp_prob = largest_range_helix ** (-1.44)  # From David et. al.
    else:
        min_bp_prob = 1.0

    return kinetic_score, min_bp_prob


def calc_dG_standby_site(structure_old, dangles, standby_site_length, constraints, rRNA):
    """Calculates the dG_standby given the structure of the mRNA:rRNA complex

    To calculate the mfe structure while disallowing base pairing at the standby site,
    we split the folded mRNA sequence into three parts: (i) a pre-sequence (before the standby
    site) that can fold; (ii) the standby site, which can not fold; (iii) the 16S rRNA binding
    site and downstream sequence, which has been previously folded.
    """


    structure = deepcopy(structure_old)
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

    #Create the mRNA subsequence
    mRNA_subsequence = mRNA[0:max(0,most_5p_mRNA - standby_site_length - 1)]
    if constraints is None:
        constraint_subsequence = None
    else:
        constraint_subsequence = constraints[0:max(0,most_5p_mRNA - standby_site_length - 1)]

    # Fold it and extract the base pairings
    if (len(mRNA_subsequence)) > 0:
        mfe_basepairing_x, mfe_basepairing_y, _ = mfe([mRNA_subsequence], constraint_subsequence, temp=ostir_constants.temp, dangles=dangles)
        bp_x_5p = mfe_basepairing_x
        bp_y_5p = mfe_basepairing_y
    else:
        bp_x_5p = []
        bp_y_5p = []

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
    energy_after = energy([mRNA, rRNA], bp_x_after, bp_y_after, dangles=dangles, Temp=ostir_constants.temp)

    dG_standby_site = energy_before - energy_after
    if dG_standby_site > 0.0:
        dG_standby_site = 0.0
    # catch negative 0

    index = structure["MinStructureID"]


    structure["bp_x"] = np.array(bp_x_after, dtype=np.intc)
    structure["bp_y"] = np.array(bp_y_after, dtype=np.intc)
    structure["subopt_energy"][index] = energy_after
    structure["dG_mRNA_rRNA_corrected"] = energy_after

    return (dG_standby_site, structure)


def calc_dG_mRNA(mRNA, start_pos, dangles, constraints):
    """Calculates the dG_mRNA given the mRNA sequence."""
    mRNA = cutoff_mRNA(mRNA, start_pos)

    fold = ViennaRNA([mRNA])
    if constraints:
        constraints = constraints[max([0,start_pos-ostir_constants.cutoff]) : min([len(mRNA), start_pos+ostir_constants.cutoff])]
        mfe_basepairing_x, mfe_basepairing_y, mfe_energy = mfe([mRNA], constraints, temp=ostir_constants.temp, dangles=dangles)
    else:
        mfe_basepairing_x, mfe_basepairing_y, mfe_energy = mfe([mRNA], None, temp=ostir_constants.temp, dangles=dangles)

    structure = fold
    structure["mRNA"] = mRNA
    structure["bp_x"] = mfe_basepairing_x
    structure["bp_y"] = mfe_basepairing_y
    structure["dG_mRNA"] = mfe_energy
    structure["MinStructureID"] = 0

    dG_mRNA_folding = mfe_energy


    kinetic_score, min_bp_prob = calc_kinetic_score(mRNA, mfe_basepairing_x, mfe_basepairing_y)

    return dG_mRNA_folding, structure, kinetic_score, min_bp_prob



def calc_dG_mRNA_rRNA(mRNA_in, rRNA, start_pos, dangles, constraints):
    """Calculates the dG_mRNA_rRNA from the mRNA and rRNA sequence.
    Considers all feasible 16S rRNA binding sites and includes the effects of non-optimal spacing."""

    # Collect all constants
    cutoff = ostir_constants.cutoff
    temp = ostir_constants.temp
    energy_cutoff = ostir_constants.energy_cutoff

    begin = max(0, start_pos - cutoff)
    mRNA_len = min(len(mRNA_in), start_pos + cutoff)
    start_pos_in_subsequence = min(start_pos, cutoff)
    startpos_to_end_len = mRNA_len - start_pos_in_subsequence - begin


    # 1. identify a list of rRNA-binding sites. Binding sites are hybridizations between the mRNA and rRNA and can include mismatches, bulges, etc. Intra-molecular folding is also allowed within the mRNA. The subopt program is used to generate a list of optimal & suboptimal binding sites.
    # Constraints: the entire rRNA-binding site must be upstream of the start codon

    mRNA = mRNA_in[begin:start_pos]
    if begin == start_pos:
        raise ValueError("Warning: There is a leaderless start codon, which is being ignored.")

    #include viennaRNA folding constraints due to binding of global regulator


    if constraints and len(constraints) >= len(mRNA):
        constraints = constraints[begin:start_pos] #added so constraints file will be the same as the sequence
    elif constraints and len(constraints) < len(mRNA):
        constraints = constraints + "."*(len(mRNA)-len(constraints)) # Fill the rest with dots
    else:
        constraints = None
    subopt_energy, subopt_basepairing_x, subopt_basepairing_y = subopt([mRNA, rRNA], constraints , energy_cutoff, dangles=dangles, temp=temp)


    if len(subopt_basepairing_x) == 0:
        return None, None, None

    # 2. Calculate dG_spacing for each 16S rRNA binding site

    # Calculate the aligned spacing for each binding site in the list
    aligned_spacing = []
    loop_len = len(subopt_basepairing_x)
    for i in range(loop_len):
        aligned_spacing.append(calc_aligned_spacing(rRNA, mRNA, start_pos_in_subsequence, subopt_basepairing_x[i], subopt_basepairing_y[i]))

    dG_spacing_list = []
    dG_mRNA_rRNA = []
    dG_mRNA_rRNA_withspacing = []

    # Calculate dG_spacing using aligned spacing value. Add it to dG_mRNA_rRNA.
    for i in range(loop_len):
        dG_mRNA_rRNA.append(subopt_energy[i])
        val = calc_dG_spacing(aligned_spacing[i])
        dG_spacing_list.append(val)
        dG_mRNA_rRNA_withspacing.append(val + subopt_energy[i])

    # 3. Find 16S rRNA binding site that minimizes dG_spacing+dG_mRNA_rRNA.
    _, index = find_min(dG_mRNA_rRNA_withspacing)
    dG_spacing_final = dG_spacing_list[index]
    spacing_value = aligned_spacing[index]

    # Check: Is the dG spacing large compared to the energy gap? If so, this means the list of suboptimal 16S rRNA binding sites generated by subopt is too short.
    if dG_spacing_final > energy_cutoff:
        if ostir_constants.verbose:
            print("Warning: The spacing penalty is greater than the energy gap. dG (spacing) = ", dG_spacing_final)

    # 4. Identify the 5' and 3' ends of the identified 16S rRNA binding site. Create a base pair list.



    most_5p_mRNA = ostir_constants.infinity
    most_3p_mRNA = -ostir_constants.infinity

    # Generate a list of rRNA-mRNA base paired nucleotides
    bp_x_target = []
    bp_y_target = []

    bp_x = subopt_basepairing_x[index]
    bp_y = subopt_basepairing_y[index]
    loop_len = len(bp_x)


    for i in range(loop_len):
        nt_x = bp_x[i]
        nt_y = bp_y[i]
        if nt_y > len(mRNA):  # nt is rRNA
            index_match = bp_y.index(nt_y)
            most_5p_mRNA = min(most_5p_mRNA, bp_x[index_match])
            most_3p_mRNA = max(most_3p_mRNA, bp_x[index_match])
            bp_x_target.append(nt_x)
            bp_y_target.append(nt_y)


    """
    The rRNA-binding site is between the nucleotides at positions most_5p_mRNA and most_3p_mRNA
    Now, fold the pre-sequence, rRNA-binding-sequence and post-sequence separately. 
    Take their base pairings and combine them together. Calculate the total energy. 
    For secondary structures, this splitting operation is allowed.
    We postulate that not all of the post-sequence can form secondary structures. 
    Once the 30S complex binds to the mRNA, it prevents the formation of secondary 
    structures that are mutually exclusive with ribosome binding. We define self.footprint 
    to be the length of the 30S complex footprint. Here, we assume that the entire mRNA 
    sequence downstream of the 16S rRNA binding site can not form secondary structures.
    """


    begin = int(begin)
    most_5p_mRNA = int(most_5p_mRNA)
    most_3p_mRNA = int(most_3p_mRNA)

    mRNA_pre = mRNA_in[begin:begin + most_5p_mRNA - 1]
    post_window_end = mRNA_len + 1
    post_window_begin = min(start_pos + ostir_constants.footprint, post_window_end)  # Footprint
    post_window_end = mRNA_len + 1
    mRNA_post = mRNA_in[post_window_begin:post_window_end]

    total_bp_x = []
    total_bp_y = []

    if constraints:
        pre_constraints = constraints[begin:begin+most_5p_mRNA-1]
        post_constraints = constraints[post_window_begin:post_window_end]
    else:
        pre_constraints = None
        post_constraints = None

    # Calculate pre-sequence folding
    if len(mRNA_pre) > 0:
        mfe_basepairing_x, mfe_basepairing_y, _ = mfe([mRNA_pre], pre_constraints, temp=temp, dangles=dangles, )
        bp_x_pre = mfe_basepairing_x
        bp_y_pre = mfe_basepairing_y

    else:
        bp_x_pre = []
        bp_y_pre = []

    # Add pre-sequence base pairings to total base pairings
    total_bp_x.extend(bp_x_pre)
    total_bp_y.extend(bp_y_pre)


    # Add rRNA-binding site base pairings to total base pairings
    if startpos_to_end_len < cutoff:
        rRNA_offset = startpos_to_end_len
    else:
        rRNA_offset = startpos_to_end_len

    total_bp_x.extend(bp_x_target)
    for nt_y in bp_y_target:
        total_bp_y.append(nt_y + rRNA_offset)

    # Calculate post-sequence folding
    if len(mRNA_post) > 0:
        mfe_basepairing_x, mfe_basepairing_y, _ = mfe([mRNA_post], post_constraints, temp=temp, dangles=dangles)
        bp_x_post = mfe_basepairing_x
        bp_y_post = mfe_basepairing_y
    else:
        bp_x_post = []
        bp_y_post = []

    offset = 0  # Begins at zero
    offset = post_window_begin - begin
    for (nt_x, nt_y) in zip(bp_x_post, bp_y_post):
        total_bp_x.append(nt_x + offset)
        total_bp_y.append(nt_y + offset)


    mRNA = mRNA_in[begin:mRNA_len]
    fold = ViennaRNA([mRNA, rRNA])

    total_energy = energy([mRNA, rRNA], total_bp_x, total_bp_y, Temp=temp, dangles=dangles)

    total_energy_withspacing = total_energy + dG_spacing_final

    structure = fold
    #structure["program"] = "subopt"
    structure["mRNA"] = mRNA
    structure["MinStructureID"] = 0
    structure["dG_mRNA_rRNA"] = total_energy
    structure["dG_mRNA_rRNA_withspacing"] = total_energy_withspacing
    structure["dG_spacing"] = dG_spacing_final
    structure["subopt_energy"] = [total_energy_withspacing]
    structure["bp_x"] = np.array(total_bp_x, dtype=np.intc)
    structure["bp_y"] = np.array(total_bp_y, dtype=np.intc)

    return total_energy_withspacing, structure, spacing_value


def calc_dG_spacing(aligned_spacing):
    """Calculates the dG_spacing according to the value of the aligned spacing.
    This relationship was determined through experiments."""

    optimal_spacing = ostir_constants.optimal_spacing
    push_constants = ostir_constants.dG_spacing_constant_push
    pull_constants = ostir_constants.dG_spacing_constant_pull

    if aligned_spacing < optimal_spacing:
        ds = aligned_spacing - optimal_spacing

        dG_spacing_penalty = push_constants[0] / (
                    1.0 + math.exp(push_constants[1] * (ds + push_constants[2]))) ** \
                            push_constants[3]

    else:
        ds = aligned_spacing - optimal_spacing
        dG_spacing_penalty = pull_constants[0] * ds * ds + pull_constants[1] * ds + \
                            pull_constants[2]

    return dG_spacing_penalty


def calc_aligned_spacing(rRNA, mRNA, start_pos, bp_x, bp_y):
    """Calculates the aligned spacing between the 16S rRNA binding site and the start codon."""

    # rRNA is the concatenated at the end of the sequence in 5' to 3' direction
    # first: identify the farthest 3' nt in the rRNA that binds to the mRNA and return its mRNA base pairer

    rRNA_len = len(rRNA)
    Ok = False
    seq_len = len(mRNA) + rRNA_len

    loop_end = seq_len - rRNA_len

    for rRNA_nt in range(seq_len, loop_end, -1):

        if rRNA_nt in bp_y:
            rRNA_pos = bp_y.index(rRNA_nt)
            x_start = bp_x[rRNA_pos]
            if x_start < start_pos:
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
        aligned_spacing = ostir_constants.infinity

    return aligned_spacing


def calc_expression_level(dG):
    """Calculates the expression level of a given dG value."""
    K = ostir_constants.K
    RT_eff = ostir_constants.RT_eff
    return K * math.exp(-dG / RT_eff)


def find_min(input_list):
    """Finds the minimum of a list of numbers."""

    min_item = ostir_constants.infinity
    min_index = 0

    for i, item in enumerate(input_list):
        if item < min_item:
            min_item = item
            min_index = i

    return (min_item, min_index)


def find_start_codons(sequence, start_range):
    """Finds all start codons in an mRNA sequence. Creates a list."""

    seq_len = len(sequence)

    #Switch to zero-indexed positions happens here
    end_0 = min(start_range[1]-1, seq_len - 2)
    begin_0 = min(start_range[0]-1, end_0)

    for i in range(begin_0, end_0 + 1):
        codon = sequence[i:i + 3]
        if codon.upper() in ostir_constants.start_codons:
            yield (i, codon)
        else:
            pass

def cutoff_mRNA(mRNA, start_pos):
    mRNA = mRNA[max(0, start_pos - ostir_constants.cutoff):min(len(mRNA), start_pos + ostir_constants.cutoff)]
    return mRNA
