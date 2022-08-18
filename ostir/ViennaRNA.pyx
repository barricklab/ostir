import subprocess, time
import tempfile
from shutil import which
import warnings
import os
from operator import itemgetter
import numpy as np
from dataclasses import dataclass
import re
import RNA
from functools import cache
from cpython cimport array
import array

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
    install_location = os.path.dirname(os.path.realpath(__file__))

vienna_constants = ViennaConstants()


#Class that encapsulates all of the functions from ViennaRNA
class ViennaRNA(dict):


    def __init__(self, Sequence_List, str material = "rna37"):
        self.RT = 0.61597 #gas constant times 310 Kelvin (in units of kcal/mol)

        for seq in Sequence_List:
            if re.compile('[ATGCU]', re.IGNORECASE).match(seq) == None:
                error_string = "Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\"."
                raise ValueError(error_string)

        self["sequences"] = Sequence_List


def mfe(list sequences, str constraints, double temp , str dangles, int basepair=0):
    
    temp = float(temp)
    if temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

    cdef str seq_string = "&".join(sequences)
    constraints = None
    if constraints is None:
        input_string = seq_string + "\n" 
    else: 
        input_string = seq_string + "\n" + constraints + "\n"

    try:
        #write sequence file 
        handle, infile = tempfile.mkstemp()
        handle = os.fdopen(handle, "w")
        handle.write(input_string)
        handle.close()
        #Set arguments

        param_file = vienna_constants.material
        param_file = f"-P {param_file}" if param_file else ""

        if dangles == "none":
            dangles = " -d0"
        elif dangles == "some":
            dangles = " -d1"
        elif dangles == "all":
            dangles = " -d2"
            
        if constraints is None:
            args = f'--noPS {dangles} --noLP {param_file} "{infile}"'
        else:
            args = f'--noPS {dangles} --noLP -C {param_file} "{infile}"'

        #Call ViennaRNA C programs
        cmd = "RNAfold "
        output = subprocess.Popen(cmd + args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines = True) #newlines argument was added because subprocess.popen was appending b's to the list output (byte stuff? I dont totally understand)
        std_out = output.communicate()[0]
        #output.tochild.write(input_string)

        while output.poll() is None:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        #if debug == 1: print(output.stdout.read())

        #Skip the unnecessary output lines
        #line = output.stdout.read()
        #print(std_out)

        line = std_out
        words = line.split()
        #print("These are the mfe value for " + str(words))
        bracket_string = words[1]
        #print(bracket_string)
        #print("This is the bracket string " + str(bracket_string))
        bp_x, bp_y = convert_bracket_to_numbered_pairs(bracket_string)
        #print(strands, bp_x, bp_y)

        energy = float(words[len(words)-1].replace(")","").replace("(","").replace("\n",""))

        #print "Minimum free energy secondary structure has been calculated."
    finally:
        # Delete temp
        if os.path.exists(infile):
            os.remove(infile)

    cdef array.array mfe_basepairing_x = array.array('i', bp_x)
    cdef array.array mfe_basepairing_y = array.array('i', bp_y)
    cdef double mfe_energy = energy
            
    return mfe_basepairing_x, mfe_basepairing_y, mfe_energy


def subopt(sequences, constraints, energy_gap, temp = 37.0, dangles = "some", outputPS = False):

    #self["subopt_composition"] = strands

    if temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

    seq_string = "&".join(sequences)
    #seq_string = self["sequences"][0]

    #print(f'SEQ STRING: {seq_string}')
    constraints = None
    if constraints is None:
        input_string = seq_string + "\n"
    else:
        input_string = seq_string + "\n" + constraints + "&........." "\n"

    try:
        handle, infile = tempfile.mkstemp()
        handle = os.fdopen(handle, "w")
        handle.write(input_string)
        handle.close()

        #Set arguments
        param_file = vienna_constants.material
        param_file = f"-P {param_file}" if param_file else ""

        if dangles == "none":
            dangles = " -d0 "
        elif dangles == "some":
            dangles = " -d1 "
        elif dangles == "all":
            dangles = " -d2 "

        if outputPS:
            outputPS_str = ""
        else:
            outputPS_str = " -noPS "

        #Call ViennaRNA C programs

        cmd = "RNAsubopt"
        energy_gap = energy_gap + 2.481
        if constraints is None:
            args = f' -e {str(energy_gap)} -T {str(temp)} {dangles} --sorted --en-only --noLP {param_file} < "{infile}"'
        else:
            args = f' -e {str(energy_gap)} -T {str(temp)} {dangles} --sorted --en-only --noLP -C {param_file} < "{infile}"'

        #print(cmd + args + f" ({input_string})")
        output = subprocess.Popen(cmd + args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines = True) #newlines argument was added because subprocess.popen was appending b's to the list output (byte stuff? I dont totally understand)
        std_out = output.communicate()[0]
        #output.tochild.write(input_string)

        while output.poll() is None:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        if debug == 1:
            print(output.stdout.read())

        #Skip unnecessary line
        #line = output.std_out.readline()
        line = std_out

        subopt_energy = []
        subopt_basepairing_x = []
        subopt_basepairing_y = []


        #print(f"Line: {line}")
        findings = line.split("\n")
        findings = findings[1:]
        findings = [x.split() for x in findings]
        
        identified_findings = []
        if len(sequences) > 1:
            findings = process_fold_outputs(findings, multi_sequence=True)
        else:
            findings = process_fold_outputs(findings, multi_sequence=False)

        for finding in findings:

            subopt_energy.append(float(finding[1]))
            subopt_basepairing_x.append(finding[2])
            subopt_basepairing_y.append(finding[3])

    finally:
        if os.path.exists(infile):
            os.remove(infile)
    #print(subopt_energy[0], subopt_basepairing_x[0], subopt_basepairing_y[0])
    return subopt_energy, subopt_basepairing_x, subopt_basepairing_y


cpdef double energy(list sequences, list base_pairing_x, list base_pairing_y, float Temp, str dangles): 
    if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
    cdef str seq_string, cmd
    cdef object args


    seq_string = "&".join(sequences)
    strands = [len(seq) for seq in sequences]
    bracket_string = convert_numbered_pairs_to_bracket(strands,base_pairing_x,base_pairing_y)
    input_string = seq_string + "\n" + bracket_string + "\n"

    #print(seq_string,strands, bracket_string, input_string)
    try:
        handle, infile = tempfile.mkstemp()
        handle = os.fdopen(handle, "w")
        handle.write(input_string)
        handle.close()
    
        #Set arguments
        param_file = vienna_constants.material
        param_file = f"-P {param_file}" if param_file else ""

        if dangles == "none":
            dangles = "-d0"
        elif dangles == "some":
            dangles = "-d1"
        elif dangles == "all":
            dangles = "-d2"

        #Call ViennaRNA C programs
        cmd = "RNAeval"
        args = f' {dangles} {param_file} "{infile}"'
        #print(cmd + args)

        output = subprocess.Popen(cmd + args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines = True) #newlines argument was added because subprocess.popen was appending b's to the list output (byte stuff? I dont totally understand)


        while output.poll() is None:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        #if debug == 1: print output.fromchild.read()
        std_out = output.communicate()[0]
        #print(std_out)


        line = std_out
        words = line.split()
        #print(line, words)
        if not words:
            raise ValueError('Could not catch the output of RNAeval. Is Vienna installed and in your PATH?')
        energy = float(words[-1].replace("(", "").replace(")", "").replace("\n", ""))

    finally:
        if os.path.exists(infile):
            os.remove(infile)

    if type(energy) != float:
        print(type(energy), energy)
    return energy


cdef str convert_numbered_pairs_to_bracket(strands, bp_x, bp_y):

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

cdef convert_bracket_to_numbered_pairs(bracket_string):

    all_seq_len = len(bracket_string)
    #print(all_seq_len)
    bp_x = []
    bp_y = []
    strands = []

    for y in range(bracket_string.count(")")):
        bp_y.append([])
        #print(bp_y)

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
        bp_x = [pos+1 for pos in bp_x[:]] #Shift so that 1st position is 1
        bp_y = [pos+1 for pos in bp_y[:]] #Shift so that 1st position is 1
        #print("subopt bp_x " + str(bp_x))
        #print("subopt bp_y " + str(bp_y))
        bp_x = np.array(bp_x)
        bp_y = np.array(bp_y)

    return bp_x, bp_y

cdef process_fold_outputs(findings, multi_sequence = False):
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

cdef object get_paramater_file(object parameter):
    cdef object filepath
    cdef object install_location = vienna_constants.install_location
    if parameter == 'rna1999':
        filepath = f"{install_location}/rna_turner1999.par"
    elif parameter == 'rna2004':
        filepath = f"{install_location}/rna_turner2004.par"
    else:
        filepath = ''
    return filepath

@cache
def get_paramater_object(object parametersfile, double temp, str dangles, int noLP):
    cdef object filepath
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
