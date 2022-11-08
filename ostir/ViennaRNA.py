"""

Some of the codebase of this file originated from the Salis Lab's RBS calculator which is distributed under GPL3.
See <http://www.gnu.org/licenses/>.
Copyright 2008-2009 is owned by the University of California Regents. All rights reserved.
"""

from shutil import which
from operator import itemgetter
from dataclasses import dataclass
from functools import cache
import warnings
import os
import array
import sys
import re
import RNA


#  On import check dependencies
dependencies = [which('RNAfold') is not None,
                which('RNAsubopt') is not None,
                which('RNAeval') is not None]
if False in dependencies:
    warnings.warn('RBS Calculator Vienna is missing dependency ViennaRNA!')
#  End dependency check

debug = 0

@dataclass
class ViennaConstants():
    debug = 0
    RT = 0.61597
    material = 'rna2004'

vienna_constants = ViennaConstants()


#Class that encapsulates all of the functions from ViennaRNA
class ViennaRNA(dict):


    def __init__(self, Sequence_List):
        self.RT = 0.61597 #gas constant times 310 Kelvin (in units of kcal/mol)

        for seq in Sequence_List:
            if re.compile('[ATGCU]', re.IGNORECASE).match(seq) is None:
                error_string = "Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\"."
                raise ValueError(error_string)

        self["sequences"] = Sequence_List


def mfe(sequences, constraints, temp , dangles):
    '''
    Calculate the MFE of a sequence using ViennaRNA as a module
    '''

    params = get_paramater_object(vienna_constants.material, temp, dangles)
    seq_string = "&".join(sequences).upper().replace("T", "U")

    rna = RNA.fold_compound(seq_string, params)

    if constraints:
        rna.hc_add_from_db(constraints)

    vienna_mfe = rna.mfe()
    bracket_string = vienna_mfe[0]
    vienna_energy = round(vienna_mfe[1], 2)

    bp_x, bp_y = convert_bracket_to_numbered_pairs(bracket_string)

    mfe_basepairing_x = array.array('i', bp_x)
    mfe_basepairing_y = array.array('i', bp_y)
    mfe_energy = vienna_energy

    return mfe_basepairing_x, mfe_basepairing_y, mfe_energy


def subopt(sequences, constraints, energy_gap, temp = 37.0, dangles = "some"):

    #self["subopt_composition"] = strands

    if temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

    seq_string = "&".join(sequences)
    #seq_string = self["sequences"][0]

    #print(f'SEQ STRING: {seq_string}')

    params = get_paramater_object(vienna_constants.material, temp, dangles)
    seq_string = "&".join(sequences).upper().replace("T", "U")

    rna = RNA.fold_compound(seq_string, params)

    if constraints:
        if len(sequences) > 1:
            constraints = constraints + "."*len(sequences[1])
        rna.hc_add_from_db(constraints)

    vienna_subopt = rna.subopt(int((energy_gap+2.481)*100))
    subopt_output = [[str(output.structure), str(round(output.energy, 3))] for output in vienna_subopt]

    subopt_energy = []
    subopt_basepairing_x = []
    subopt_basepairing_y = []

    if len(sequences) > 1:
        subopt_output = process_fold_outputs(subopt_output, multi_sequence=True)
    else:
        subopt_output = process_fold_outputs(subopt_output, multi_sequence=False)

    for finding in subopt_output:
        
        subopt_energy.append(float(finding[1]))
        subopt_basepairing_x.append(finding[2])
        subopt_basepairing_y.append(finding[3])


    return subopt_energy, subopt_basepairing_x, subopt_basepairing_y


def energy(sequences, base_pairing_x, base_pairing_y, Temp, dangles): 
    if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
    temp = Temp

    seq_string = "&".join(sequences)
    strands = [len(seq) for seq in sequences]
    bracket_string = convert_numbered_pairs_to_bracket(strands,base_pairing_x,base_pairing_y)

    params = get_paramater_object(vienna_constants.material, temp, dangles)
    seq_string = "&".join(sequences).upper().replace("T", "U")

    rna = RNA.fold_compound(seq_string, params)
    result = rna.eval_structure(bracket_string.replace("&", ""))

    vienna_energy = round(result, 2)
    return vienna_energy


def convert_numbered_pairs_to_bracket(strands, bp_x, bp_y):

    bp_x = [pos-1 for pos in bp_x[:]] #Shift so that 1st position is 0
    bp_y = [pos-1 for pos in bp_y[:]] #Shift so that 1st position is 0

    bracket_notation = []
    counter = 0
    for (strand_number, seq_len) in enumerate(strands):
        if strand_number > 0: bracket_notation.append("&")
        for pos in range(counter, seq_len+counter):
            if pos in bp_x:
                bracket_notation.append("(")
            elif pos in bp_y:
                bracket_notation.append(")")
            else:
                bracket_notation.append(".")
        counter += seq_len

    return "".join(bracket_notation)

def convert_bracket_to_numbered_pairs(bracket_string):

    #print(all_seq_len)
    bp_x = []
    bp_y = []

    for _ in range(bracket_string.count(")")):
        bp_y.append([])

    last_nt_x_list = []
    num_strands=0
    #print(bracket_string)
    for (pos,letter) in enumerate(bracket_string[:]):
        if letter == ".":
            pass

        elif letter == "(":
            bp_x.append(pos-num_strands)
            last_nt_x_list.append(pos-num_strands)

        elif letter == ")":
            nt_x = last_nt_x_list.pop() #nt_x is list of "(" except last entry
            #print('this is the last_nt_x_list ' + str(last_nt_x_list.pop()))
            nt_x_pos = bp_x.index(nt_x)
            bp_y[nt_x_pos] = pos-num_strands

        elif letter == "&":
            num_strands += 1

        else:
            print("Error! Invalid character in bracket notation.")

    if len(last_nt_x_list) > 0:
        print("Error! Leftover unpaired nucleotides when converting from bracket notation to numbered base pairs.")

    if len(bp_y) > 1:
        bp_x = [int(pos+1) for pos in bp_x[:]] #Shift so that 1st position is 1
        bp_y = [int(pos+1) for pos in bp_y[:]] #Shift so that 1st position is 1
        #print("subopt bp_x " + str(bp_x))
        #print("subopt bp_y " + str(bp_y))

    return bp_x, bp_y

def process_fold_outputs(findings, multi_sequence = False):
    filtered_findings = []
    # Earlier vienna versions would report 'folds' without any actual folding when running under WSL.
    # This code provides backwards compatibility for versions with that bug
    for finding in findings:
        if len(finding) != 2:
            continue
        if multi_sequence:
            binding_positions = finding[0]
            binding_positions = binding_positions.split("&")
            if ")" not in binding_positions[1] and "(" not in binding_positions[1]:
                #print("This binding site is fake")
                continue
            else:
                bp_x, bp_y = convert_bracket_to_numbered_pairs(finding[0])
                #print(strands, bp_x, bp_y)
                filtered_findings.append([finding[0], finding[1], bp_x, bp_y])
        else:
            bp_x, bp_y = convert_bracket_to_numbered_pairs(finding[0])
            filtered_findings.append([finding[0], finding[1], bp_x, bp_y])
    # Linux distributions of Vienna RNA tiebreak the sorting of outputs differently than Mac. This code aims to
    # make sorting of the output consistent, which rarely affects downstream predictions by sorting by energy with
    # the string of the fold as the tiebreaker
    findings_by_energy = {}
    for finding in filtered_findings:  # Sorts by energy
        if finding[1] in findings_by_energy.keys():
            findings_by_energy[finding[1]].append(finding)
        else:
            findings_by_energy[finding[1]] = [finding]
    energy_values = [[float(i), i] for i in findings_by_energy.keys()]  # real keys may have hanging 0's
    energy_values = sorted(energy_values, key=itemgetter(1))
    sorted_findings = []
    for value in energy_values:
        dict_index = value[1]
        partitioned_findings = findings_by_energy[str(dict_index)]
        partitioned_findings = sorted(partitioned_findings, key=itemgetter(0))
        sorted_findings = sorted_findings + partitioned_findings

    return sorted_findings

def get_paramater_file(parameter):
    if parameter == 'rna1999':
        filepath = os.path.dirname(sys.modules['ostir'].__file__)+'/rna_turner1999.par'
    elif parameter == 'rna2004':
        filepath = os.path.dirname(sys.modules['ostir'].__file__)+'/rna_turner2004.par'
    else:
        filepath = ''
    return filepath

@cache
def get_paramater_object(parametersfile, temp, dangles):
    filepath = get_paramater_file(parametersfile)
    RNA.params_load(filepath)

    params = RNA.md()
    params.noLP = 1
    params.noPS = 1
    params.temperature = temp

    if dangles == "none":
        params.dangles = 0
    elif dangles == "some":
        params.dangles = 1
    elif dangles == "all":
        params.dangles = 2

    return params


def clean_temp_file(file):
    if os.path.isfile(file):
        os.remove(file)
