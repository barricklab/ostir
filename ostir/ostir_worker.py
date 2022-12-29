from dataclasses import dataclass
from .data_classes import OSTIRResult

from .ostir_calculations import calc_dG_standby_site, calc_longest_loop_bulge, calc_dG_mRNA_rRNA, calc_dG_mRNA, calc_expression_level, calc_longest_helix


# Default values. @TODO: Find a way to overwrite these
auto_dangles = True
dangles_default = "all"
start_codon_energies = {"ATG": -1.194, "AUG": -1.194, "GTG": -0.0748, "GUG": -0.0748, "TTG": -0.0435,
                                        "UUG": -0.0435, "CTG": -0.03406, "CUG": -0.03406}  # hybridization to CAT
standby_site_length = 4
decimal_places = 4
cutoff = 35

def ostir_worker(task, start_position, task_ID, verbosity=0):
    """ Run OSTIR on a single sequence. """

    # Get the start codon from the start position and the task
    sequence = task.sequence.upper()
    start_codon = start_position[1]
    start_position = start_position[0]

    # Set dangles based on length between 5' end of mRNA and start codon
    if auto_dangles:
        global cutoff
        cutoff = int(cutoff)
        if start_position > cutoff:
            dangles = "none"

        else:
            dangles = "all"

    else:
        dangles = dangles_default

    # Start codon energy
    dG_start_codon = start_codon_energies[start_codon]

    # Energy of mRNA folding. This also gets us the kinetic score
    dG_mRNA, mRNA_structure, kinetic_score, min_bp_prob = calc_dG_mRNA(sequence, start_position, dangles, task.constraints)

    # Energy of mRNA:rRNA hybridization & folding
    dg_mRNA_rRNA_output = calc_dG_mRNA_rRNA(sequence, task.rrna, start_position, dangles, task.constraints)
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
                                                                            standby_site_length,
                                                                            task.constraints,
                                                                            task.rrna)

    # Total energy is mRNA:rRNA + start - rRNA - mRNA - standby_site
    dG_total = dG_mRNA_rRNA_withspacing + dG_start_codon - dG_mRNA - dG_standby_site

    # Calculate dG to open SD sequence
    # ddG_SD_open = self.calc_dG_SDopen(mRNA_structure, mRNA_rRNA_structure)

    loop_result = calc_longest_loop_bulge(structure=mRNA_structure)
    helical_loop_list, bulge_loop_list = loop_result[0], loop_result[1]

    parallel_output = OSTIRResult(
        name= task.name,
        start_codon=start_codon,
        start_position=start_position+1,
        expression=round(calc_expression_level(dG_total),
                            decimal_places),
        RBS_distance_bp=int(spacing_value),
        dG_total=round(float(dG_total), decimal_places),
        dG_rRNA_mRNA=round(float(dG_mRNA_rRNA_nospacing), decimal_places),
        dG_mRNA=round(float(dG_mRNA), decimal_places),
        dG_spacing=round(float(mRNA_rRNA_structure["dG_spacing"]), decimal_places),
        dG_standby=round(float(dG_standby_site), decimal_places),
        dG_start_codon=round(float(dG_start_codon), decimal_places),

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


    return task_ID, parallel_output