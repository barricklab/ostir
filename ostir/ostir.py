#!/usr/bin/env python

#Jeff: TODO
# Move start/end logic into run_ostir


import argparse
import math
import os
import sys
from pathlib import Path
import subprocess
import csv
import re
from shutil import which


try:
    from ostir.ostir_factory import OSTIRFactory
except ModuleNotFoundError:
    from ostir_factory import OSTIRFactory

ostir_version = '1.0.6'
oldest_vienna = '2.4.18'

# The E. coli sequence
Ecoli_anti_Shine_Dalgarno = 'ACCTCCTTA'

def run_ostir(in_seq, start=None, end=None, name=None, aSD=None, threads=1, decimal_places=4, verbose=False):
    '''Takes an RNA with optional parameters and returns binding energies.
        Keyword arguments:
        seq -- Sequence to calculate binding energies for
        outfile -- Filepath for output csv
        start -- First base to start considering start cofdons. 1-indexed. Defaults to first base
        end -- Last base to start considering start codons. 1-indexed.
                     If start_loc_1 is provided, defaults to start_loc_1
                     Otherwise, defaults to end of sequence
        name -- Returns itself, useful for tagging things for downstream processing
        aSD -- Defines anti-Shine-Dalgarno sequence: the 9 bp at the 3' end of the 16S rRNA. Defaults to the E. coli sequence.
        threads -- Defines parallel processing workers, roughly equivalent to multithreading cores
        decimal_places -- Precision of numerical output (number of places to the right of the decimal)
        verbose -- Prints debug information
    '''

    #rename to more descriptive vars
    in_start_loc_1 = start
    in_end_loc_1 = end

    # Set up a default name and sd if not provided
    if name == None:
        name = 'unnamed'

    if aSD == None:
        aSD = Ecoli_anti_Shine_Dalgarno

    # Do some checks on the inputs

    # Nucleotide character check
    nucleotides=re.compile('[^ATCGUatcgu]')
    if (nucleotides.search(aSD) != None):
        print(f"ERROR: anti-Shine-Dalgarno sequence provided ({aSD}) contains non-nucleotide characters.\n<<<Sequence ({name}) will be skipped.>>>", file=sys.stderr)
        return []

    # Length check
    if len(aSD) != 9:
        print(f"ERROR: anti-Shine-Dalgarno sequence provided ({aSD}) is not 9 bases.\n<<<Sequence ({name}) will be skipped.>>>", file=sys.stderr)
        return []

    #Upper case and convert to RNA
    aSD = aSD.upper()
    aSD = aSD.replace("T","U")

    #Clean spaces from sequence... could clean other characters too
    seq = in_seq.replace(" ", "")

    # Nucleotide character check
    if (nucleotides.search(seq) != None):
        print(f"ERROR: Input sequence contains non-nucleotide characters.\n<<<Sequence ({name}) will be skipped.>>>", file=sys.stderr)
        return []


    # Start <= end check
    if (in_start_loc_1!=None and in_end_loc_1!=None and in_end_loc_1<in_start_loc_1):
        print(f"ERROR: Start location ({in_start_loc_1}) is not less than end location ({in_end_loc_1}).\n<<<Sequence ({name}) will be skipped.>>>", file=sys.stderr)
        return []

    # Set up start and end locations to search
    start_loc_1 = in_start_loc_1
    if start_loc_1==None:
        start_loc_1 = 1
    start_loc_1 = int(start_loc_1)

    end_loc_1 = in_end_loc_1
    if end_loc_1==None:
        if in_start_loc_1==None:
            end_loc_1 = len(seq)
        else:
            end_loc_1 = start_loc_1
    end_loc_1 = int(end_loc_1)

    start_range_1 = [start_loc_1, end_loc_1]


    #for debugging
    #print(in_seq)
    #print(str(start_range_1[0]))
    #print(str(start_range_1[1]))
    #print(name)
    #print(sd)

    calcObj = OSTIRFactory(seq, start_range_1, aSD, verbose=verbose)
    calcObj.threads = threads
    calcObj.calc_dG()

    dG_total_list = calcObj.dG_total_list[:]
    dG_details = calcObj.dG_details[:]
    start_pos_list = calcObj.start_pos_list[:]
    kinetic_score_list = calcObj.kinetic_score_list[:]
    standby_site_list = calcObj.dG_standby_site_list[:]
    start_codon_list = calcObj.start_codon_list[:]

    expr_list = []
    for dG in dG_total_list:
        expr_list.append(calcObj.K * math.exp(-dG / calcObj.RT_eff))
    return_var = zip(expr_list, start_pos_list, kinetic_score_list, dG_total_list, standby_site_list, start_codon_list)
    return_var = list(return_var)
    zip_output = zip(return_var, dG_details)
    output_data_list = []
    for output in list(zip_output):
        outdata = {
            'name': name,
            #'sequence': seq,   #Don't pass back and print sequence, this will be very big for some inputs
            'start_codon': output[0][5],
            'start_position': output[0][1]+1, #Convert back to 1-indexed
            'expression': round(output[0][0], decimal_places),
            'RBS_distance_bp': output[1][5],
            'dG_total': round(output[0][3], decimal_places),
            'dG_rRNA:mRNA': round(output[1][0], decimal_places),
            'dG_mRNA': round(output[1][2], decimal_places),
            'dG_spacing': round(output[1][4], decimal_places),
            'dG_standby': round(output[1][3], decimal_places),
            'dG_start_codon': round(output[1][1], decimal_places),
        }
        output_data_list.append(outdata)
    output_data_list = sorted(output_data_list, key=lambda x: x['start_position'])

    return output_data_list

def parse_fasta(filepath):
    '''Takes a filepath to a fasta formatted file and returns a list of [header, sequence].'''
    sequences = []
    current_seq_name = None
    current_seq = ""
    with open(filepath, 'r') as infile:
        for line in infile:
            linestr = str(line)
            if linestr[0] == '>':
                if current_seq_name:
                    sequences.append([current_seq_name, current_seq])
                current_seq_name = linestr[1:].rstrip()
                current_seq = str()
                continue
            else:
                current_seq += linestr.rstrip()
    if current_seq_name:
        sequences.append([current_seq_name, current_seq])
    return sequences

def _print_output(outdict):
    '''Processes command line / csv input for pretty command line output'''
    sorted_predictions = {}
    keys = []
    for prediction in outdict:
        if prediction['name'] in sorted_predictions.keys():
            sorted_predictions[prediction['name']].append(prediction)
        else:
            sorted_predictions[prediction['name']] = [prediction]
            out_names = []
            out_names.append(prediction['name'])
            if prediction.get('i'):
                out_names.append(prediction.get('i'))
            keys.append(out_names)

    if not keys:
        print('No binding sites were identified.')
        exit(0)

    output_items = ['start_codon', 'start_position', 'expression', 'RBS_distance_bp', 'dG_total', 'dG_rRNA:mRNA', 'dG_mRNA', 'dG_spacing', 'dG_standby', 'dG_start_codon']
    row_format = "{:>16}" * (len(output_items))
    print('_________________________________________________')
    for rna in keys:
        rna = rna[0]
        if len(rna) == 2:
            print(f'Tested Sequence: {rna[1]}')
            print(f'Sequence RNA: {rna[0]}')
        elif len(rna) == 1:
            print(f'Sequence RNA: {rna[0]}')
        print(row_format.format(*output_items))
        for start in sorted_predictions[rna]:
            output_data = [start[key] for key in output_items]
            for i, data_point in enumerate(output_data):
                if type(data_point) == float:
                    data_point = format(data_point, '.4f')
                    output_data[i] = data_point
            print(row_format.format(*output_data))
        print('_________________________________________________')

def save_to_csv(column_names, outdict, outfile):
    with open(outfile, 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, column_names, lineterminator="\n")
        dict_writer.writeheader()
        dict_writer.writerows(outdict)

# Used to get rid of blank and comment lines on file input into csv.reader
class BlankCommentCSVFile:
    def __init__(self, fp):
        self.fp = fp

    def __iter__(self):
        return self

    def __next__(self):
        line = self.fp.__next__()
        if not line.strip() or line[0] == "#":
            return self.__next__()
        return line

def main():
    '''Main OSTIR function ran when called as an executable. See ostir -h.'''
    parser = argparse.ArgumentParser(description=f'OSTIR (Open Source Translation Initiation Rates) version {ostir_version}')

    parser.add_argument(
        '-i', '--input',
        action='store',
        metavar='str/filepath',
        dest='i',
        required=False,
        type=str,
        help="Input filename (FASTA/CSV) or DNA/RNA sequence. For CSV input files, there must be a 'seq' or 'sequence' column. Other columns will override any options provided at the command line if they are present: 'name/id', 'start', 'end', 'anti-Shine-Dalgarno'",
    )

    parser.add_argument(
        '-o', '--output',
        action='store',
        metavar= 'filepath',
        dest='o',
        required=False,
        type=str,
        help="Output file path. If not provided, results will output to the console",
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        dest='v',
        required=False,
        help="Prints verbose output.",
    )

    parser.add_argument(
        '-s', '--start',
        action='store',
        metavar='int',
        dest='s',
        required=False,
        type=int,
        help="Most 5' nucleotide position to consider a start codon beginning (1-indexed)",
    )

    parser.add_argument(
        '-e', '--end',
        action='store',
        metavar='int',
        dest='e',
        required=False,
        type=int,
        help="Most 3' nucleotide position to consider a start codon beginning (1-indexed)",
    )

    parser.add_argument(
        '-a', '--anti-Shine-Dalgarno',
        action='store',
        metavar='str',
        dest='a',
        required=False,
        type=str,
        help=f"anti-Shine-Dalgarno sequence: the 9 bases located at the 3' end of 16S rRNA. May be provided as DNA or RNA. Defaults to that of E. coli ({Ecoli_anti_Shine_Dalgarno}).",
    )

    parser.add_argument(
        '-p', '--print-sequence',
        action='store_true',
        dest='p',
        required=False,
        help="Include the input mRNA sequence in output CSV files",
    )

    parser.add_argument(
        '-q', '--print-anti-Shine-Dalgarno',
        action='store_true',
        dest='q',
        required=False,
        help="Include the anti-Shine-Dalgarno sequence in output CSV files",
    )

    parser.add_argument(
        '-j', '--threads',
        action='store',
        metavar='int',
        dest='j',
        required=False,
        type=int,
        help="Number of threads for multiprocessing",
    )

    parser.add_argument(
        '-t', '--type',
        action='store',
        metavar='[string|csv|fasta]',
        dest='t',
        required=False,
        type=str,
        help="Input type (overrides autodetection)",
    )

    parser.add_argument(
        '--version',
        action='store_true',
        dest='version',
        required=False,
        help="Print version and quit.",
    )

    options = parser.parse_args()

    if options.version:
        print(f'OSTIR version {ostir_version}', file=sys.stderr)
        exit(1)

    if not options.i:
        print("Input (-i) required.")
        parser.print_help()
        exit(1)

    cmd_kwargs = dict()
    cmd_kwargs['seq'] = options.i
    cmd_kwargs['detailed_out'] = True
    if options.o:
        outfile = options.o
    else:
        outfile = None
        cmd_kwargs['print_out'] = False
    if options.j:
        threads = options.j
    else:
        threads = 1
    if options.s:
        cmd_kwargs['start'] = options.s
    if options.e:
        cmd_kwargs['end'] = options.e
    if options.a:
        cmd_kwargs['aSD'] = options.a
    if options.p:
        cmd_kwargs['print_mRNA_sequence'] = options.p
    if options.q:
        cmd_kwargs['print_aSD_sequence'] = options.q

    # Check if viennaRNA is installed
    dependencies = [which('RNAfold') is not None,
                    which('RNAsubopt') is not None,
                    which('RNAeval') is not None]

    if False in dependencies:
        raise EnvironmentError('ViennaRNA is not properly installed or in PATH')

    vienna_version = subprocess.check_output(['RNAfold', '--version'])
    vienna_version = str(vienna_version.strip()).replace("'", "").split(' ')[1]
    print(f'Running OSTIR version {ostir_version} (with Vienna version: {vienna_version})', file=sys.stderr)

    # Check if the viennaRNA version is recent enough
    vienna_version_split = vienna_version.split('.')
    global oldest_vienna
    oldest_vienna_split = oldest_vienna.split('.')
    warning_string = f'The installed version of ViennaRNA ({vienna_version}) is older than what is supported ({oldest_vienna}).'
    for i in range(0, len(vienna_version_split)):
        if vienna_version_split[i] < oldest_vienna_split[i]:
            raise EnvironmentError(warning_string)
    # Output data: RNA, Codon, position, dg_total, dg rRNA:mRNA, dg mRNA, dG Spacing, dg Standby, Kinetic Score

    # Determine input type
    input_type = None
    specified_file_exists = False
    valid_string_check = re.compile('[ATGCU.-]', re.IGNORECASE)
    if options.t:
        if options.t == 'fasta':
            input_type = 'fasta'
        elif options.t == 'csv':
            input_type = 'csv'
        elif options.t == 'string':
            input_type = 'string'
        else:
            print(f'Unsupported file type {options.t}.')
            exit(1)
    elif os.path.isfile(cmd_kwargs['seq']) and not input_type:
        filepath_test = Path(cmd_kwargs['seq']).suffix
        if filepath_test == '.fasta':
            input_type = 'fasta'
        elif filepath_test == '.fa':
            input_type = 'fasta'
        elif filepath_test == '.fna':
            input_type = 'fasta'
        elif filepath_test == '.csv':
            input_type = 'csv'
        else:
            with open(cmd_kwargs['seq'], 'r') as in_file:
                specified_file_exists = True
                first_line = in_file.readline()
                if first_line[0] == '>':
                    input_type = 'fasta'
                elif 'seq' in first_line[0]:
                    input_type = 'csv'
    elif valid_string_check.match(cmd_kwargs['seq']):
        input_type = 'string'



    if input_type == None:
        if specified_file_exists:
            print(f'Unable to identify the type of file specified as inout (-i). Please define it using "-t".', file=sys.stderr)
        else:
            print(f'Fix input (-i). Provided value does not specify an existing file and is not a valid nucleotide sequence.', file=sys.stderr)
        exit(1)


    # Run OSTIR
    results = []

    ## String input ##############################################################
    if input_type == 'string':

        print(f'Reading input sequence from command line', file=sys.stderr)

        sequence = cmd_kwargs.get('seq')
        name = None
        start_loc_1 = cmd_kwargs.get('start')
        end_loc_1 = cmd_kwargs.get('end')
        aSD = cmd_kwargs.get('aSD')
        verbose = False

        output_dict_list = run_ostir(sequence, start=start_loc_1, end=end_loc_1, name=name, aSD=aSD, threads=threads, verbose=verbose)

        for output_dict in output_dict_list:
            if cmd_kwargs.get('print_mRNA_sequence'):
                output_dict['sequence'] = sequence
            if cmd_kwargs.get('print_aSD_sequence'):
                output_dict['anti-Shine-Dalgarno'] = aSD

        results.extend(output_dict_list)

    ## FASTA input ##############################################################
    elif input_type == 'fasta':
        input_file = cmd_kwargs['seq']
        print(f'Reading FASTA file {input_file}', file=sys.stderr)
        sequence_entries = parse_fasta(input_file)
        for sequence_entry in sequence_entries:
            sequence = sequence_entry[1]
            name = sequence_entry[0]
            start_loc_1 = cmd_kwargs.get('start')
            end_loc_1 = cmd_kwargs.get('end')
            aSD = cmd_kwargs.get('aSD')
            verbose = False

            output_dict_list = run_ostir(sequence, start=start_loc_1, end=end_loc_1, name=name, aSD=aSD, threads=threads, verbose=verbose)

            for output_dict in output_dict_list:
                if cmd_kwargs.get('print_mRNA_sequence'):
                    output_dict['sequence'] = sequence
                if cmd_kwargs.get('print_aSD_sequence'):
                    output_dict['anti-Shine-Dalgarno'] = aSD

            results.extend(output_dict_list)


    ## CSV input ##############################################################
    elif input_type == 'csv':
        input_file = cmd_kwargs['seq']
        print(f'Reading CSV file {input_file}', file=sys.stderr)
        reader = csv.DictReader(BlankCommentCSVFile(open(input_file, 'r', encoding='UTF-8-sig')))
        on_seq_index = 0 #used for giving names if none provided in input file
        for row in reader:
            on_seq_index = on_seq_index+1

            #lowercase the keys
            row = dict((k.lower(), v) for k, v in row.items())

            # for debugging
            #print(row)

            ## Allow 'sequence' or 'seq'
            if 'seq' in row.keys():
                sequence = row['seq']
            elif 'sequence' in row.keys():
                sequence = row['sequence']
            else:
                print(f"Required column 'sequence' or 'seq' not found for CSV file row: {row}")
                exit(1)

            # Assign a name if one is not given from name/id columns
            # If empty assign one based on the index
            name = row.get('name')
            if name == None or not name:
                name = row.get('id')
            if name == None or not name:
                name="sequence_" + str(on_seq_index)

            aSD = row.get('anti-shine-dalgarno') #remember, keys lowercased here
            start_loc_1 = row.get('start')
            end_loc_1 = row.get('end')

            # If any of these are not defined or are empty, use command line values as the defaults
            if not aSD:
                aSD = cmd_kwargs.get('aSD')
            if not start_loc_1:
                start_loc_1 = cmd_kwargs.get('start')
            if not end_loc_1:
                end_loc_1 = cmd_kwargs.get('end')

            verbose = False

            output_dict_list = run_ostir(sequence, start=start_loc_1, end=end_loc_1, name=name, aSD=aSD, threads=threads, verbose=verbose)

            for output_dict in output_dict_list:
                if cmd_kwargs.get('print_mRNA_sequence'):
                    output_dict['sequence'] = sequence
                if cmd_kwargs.get('print_aSD_sequence'):
                    output_dict['anti-Shine-Dalgarno'] = aSD

            results.extend(output_dict_list)

    ## Output - for all ways of running ##############################################################
    if outfile:
        #if the output is empty, we need to print a valid header...
        column_names = ['name', 'start_codon', 'start_position', 'expression', 'RBS_distance_bp', 'dG_total', 'dG_rRNA:mRNA', 'dG_mRNA', 'dG_spacing', 'dG_standby', 'dG_start_codon']
        if cmd_kwargs.get('print_mRNA_sequence'):
            column_names.insert(1, 'sequence')
        if cmd_kwargs.get('print_aSD_sequence'):
            column_names.insert(1, 'anti-Shine-Dalgarno')
        save_to_csv(column_names, results, outfile)
        print(f'Results written to {outfile}', file=sys.stderr)
    else:
        _print_output(results)


if __name__ == "__main__":
    main()
