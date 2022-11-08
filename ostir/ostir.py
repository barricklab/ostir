#!/usr/bin/env python
'''This module is the main module for the OSTIR package. It contains the run_ostir function for
high level execution and the CLI interface. To use, run "ostir" in the command line or import
ostir and run run_ostir.

OSTIR is distributed under GPL3. See <http://www.gnu.org/licenses/>.'''


import argparse
import os
import sys
import subprocess
import csv
import re
from shutil import which
from warnings import warn
from pathlib import Path


try:
    from ostir.ostir_factory import OSTIRFactory
except ModuleNotFoundError:
    from .ostir_factory import OSTIRFactory

OSTIR_VERSION = '1.1.0'
OLDEST_VIENNA = '2.4.18'

# The E. coli sequence
Ecoli_anti_Shine_Dalgarno = 'ACCTCCTTA'

def run_ostir(in_seq, start=None, end=None, name=None, aSD=None, threads=1, decimal_places=4, circular=False, constraints=None, verbose=False):
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
        constraints -- Folding constraints passed to vienna starting from the cuttoff window (35bp upstream of start codon)
        verbose -- Prints debug information
    '''
    #rename to more descriptive vars
    in_start_loc_1 = start
    in_end_loc_1 = end

    # Set up a default name and sd if not provided
    if name is None:
        name = 'unnamed'

    if aSD is None:
        aSD = Ecoli_anti_Shine_Dalgarno

    # Do some checks on the inputs

    # Nucleotide character check
    nucleotides=re.compile('[^ATCGUatcgu]')
    if nucleotides.search(aSD) is not None:
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
    if nucleotides.search(seq) is not None:
        print(f"ERROR: Input sequence contains non-nucleotide characters.\n<<<Sequence ({name}) will be skipped.>>>", file=sys.stderr)
        return []


    # Start <= end check
    if in_start_loc_1 is not None and in_end_loc_1 is not None and in_end_loc_1<in_start_loc_1:
        print(f"ERROR: Start location ({in_start_loc_1}) is not less than end location ({in_end_loc_1}).\n<<<Sequence ({name}) will be skipped.>>>", file=sys.stderr)
        return []

    # Set up start and end locations to search
    start_loc_1 = in_start_loc_1
    if start_loc_1 is None:
        start_loc_1 = 1
    start_loc_1 = int(start_loc_1)

    end_loc_1 = in_end_loc_1
    if end_loc_1 is None:
        if in_start_loc_1 is None:
            end_loc_1 = len(seq)
        else:
            end_loc_1 = start_loc_1
    end_loc_1 = int(end_loc_1)

    start_range_1 = [start_loc_1, end_loc_1]

    calcObj = OSTIRFactory(seq, start_range_1, aSD, constraints, circular=circular, verbose=verbose)
    calcObj.threads = threads
    calcObj.decimal_places = decimal_places
    calcObj.name = name
    calcObj.calc_dG()

    output_data_list = [result.results() for result in calcObj.results]

    output_data_list = sorted(output_data_list, key=lambda x: x['start_position'])

    del calcObj  # Send OSTIRFactory to garbage collection, clears temp files

    return output_data_list

def parse_fasta(filepath):
    '''Takes a filepath to a fasta formatted file and returns a list of [header, sequence].'''
    sequences = []
    current_seq_name = None
    current_seq = ""
    with open(filepath, 'r', encoding='utf8') as infile:
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
        if prediction['name'] in sorted_predictions:
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
        sys.exit(0)

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
                if isinstance(data_point, float):
                    data_point = format(data_point, '.4f')
                    output_data[i] = data_point
            print(row_format.format(*output_data))
        print('_________________________________________________')

def save_to_csv(column_names, outdict, outfile):
    with open(outfile, 'w', encoding='utf8') as output_file:
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
    parser = argparse.ArgumentParser(description=f'OSTIR (Open Source Translation Initiation Rates) version {OSTIR_VERSION}')

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
        '-c', '--circular',
        action='store_true',
        dest='c',
        required=False,
        help="Flag the input as circular",
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
        print(f'OSTIR version {OSTIR_VERSION}', file=sys.stderr)
        sys.exit(0)

    if not options.i:
        print("Input (-i) required.")
        parser.print_help()
        sys.exit(1)

    cmd_kwargs = dict()
    cmd_kwargs['seq'] = options.i
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
    if options.c:
        cmd_kwargs['circular'] = options.c

    # Check if viennaRNA is installed
    dependencies = [which('RNAfold') is not None,
                    which('RNAsubopt') is not None,
                    which('RNAeval') is not None]

    if False in dependencies:
        raise EnvironmentError('ViennaRNA is not properly installed or in PATH')

    vienna_version = subprocess.check_output(['RNAfold', '--version'])
    vienna_version = str(vienna_version.strip()).replace("'", "").split(' ')[1]
    print(f'Running OSTIR version {OSTIR_VERSION} (with Vienna version: {vienna_version})', file=sys.stderr)

    # Check if the viennaRNA version is recent enough
    vienna_version_split = vienna_version.split('.')
    global OLDEST_VIENNA
    OLDEST_VIENNA_split = OLDEST_VIENNA.split('.')
    warning_string = f'The installed version of ViennaRNA ({vienna_version}) is older than what is supported ({OLDEST_VIENNA}).'
    for i, split_digit in enumerate(vienna_version_split):
        if split_digit > OLDEST_VIENNA_split[i]:
            warn(f'The installed version of ViennaRNA ({vienna_version}) is new than what is last validated ({OLDEST_VIENNA}).')
            break
        if split_digit < OLDEST_VIENNA_split[i]:
            raise EnvironmentError(warning_string)
    # Output data: RNA, Codon, position, dg_total, dg rRNA:mRNA, dg mRNA, dG Spacing, dg Standby, Kinetic Score

    # Determine input type
    input_type = None
    specified_file_exists = False
    valid_string_check = re.compile('[ATGCU.-]', re.IGNORECASE)
    if options.t:  # Manual override
        if options.t == 'fasta':
            input_type = 'fasta'
        elif options.t == 'csv':
            input_type = 'csv'
        elif options.t == 'string':
            input_type = 'string'
        else:
            print(f'Unsupported file type {options.t}.')
            sys.exit(1)

    elif os.path.isfile(cmd_kwargs['seq']) and not input_type:  # Get input type from file name
        filepath_test = Path(cmd_kwargs['seq']).suffix
        if filepath_test in ['.fasta', '.fa', '.fna']:
            input_type = 'fasta'
        elif filepath_test == '.csv':
            input_type = 'csv'
        else:
            with open(cmd_kwargs['seq'], 'r', encoding='utf8') as in_file:   # Try to determine input type from file contents
                specified_file_exists = True
                first_line = in_file.readline()
                if first_line[0] == '>':
                    input_type = 'fasta'
                elif 'seq' in first_line[0]:
                    input_type = 'csv'

    elif valid_string_check.match(cmd_kwargs['seq']):  # Check to see if input is a valid sequence
        input_type = 'string'

    if input_type is None:   # Error out
        if specified_file_exists:
            print('Unable to identify the type of file specified as inout (-i). Please define it using "-t".', file=sys.stderr)
        else:
            print('Fix input (-i). Provided value does not specify an existing file and is not a valid nucleotide sequence.', file=sys.stderr)
        sys.exit(1)


    # Run OSTIR
    results = []

    ## String input ##############################################################
    if input_type == 'string':

        print('Reading input sequence from command line', file=sys.stderr)

        sequence = cmd_kwargs.get('seq')
        name = None
        start_loc_1 = cmd_kwargs.get('start')
        end_loc_1 = cmd_kwargs.get('end')
        aSD = cmd_kwargs.get('aSD')
        verbose = False
        circular = cmd_kwargs.get('circular')

        output_dict_list = run_ostir(sequence,
                                     start=start_loc_1,
                                     end=end_loc_1,
                                     name=name,
                                     aSD=aSD,
                                     threads=threads,
                                     circular=circular,
                                     verbose=verbose)

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
            circular = cmd_kwargs.get('circular')
            verbose = False

            output_dict_list = run_ostir(sequence,
                                         start=start_loc_1,
                                         end=end_loc_1,
                                         name=name,
                                         aSD=aSD,
                                         threads=threads,
                                         circular=circular,
                                         verbose=verbose)

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
                sys.exit(1)

            # Assign a name if one is not given from name/id columns
            # If empty assign one based on the index
            name = row.get('name')
            if not name:
                name = row.get('id')
            if not name:
                name="sequence_" + str(on_seq_index)

            aSD = row.get('anti-shine-dalgarno') #remember, keys lowercased here
            start_loc_1 = row.get('start')
            end_loc_1 = row.get('end')
            circular = row.get('circular')

            # If any of these are not defined or are empty, use command line values as the defaults
            if not aSD:
                aSD = cmd_kwargs.get('aSD')
            if not start_loc_1:
                start_loc_1 = cmd_kwargs.get('start')
            if not end_loc_1:
                end_loc_1 = cmd_kwargs.get('end')
            if not circular:
                circular = cmd_kwargs.get('circular')

            verbose = False

            output_dict_list = run_ostir(sequence,
                                         start=start_loc_1,
                                         end=end_loc_1,
                                         name=name,
                                         aSD=aSD,
                                         threads=threads,
                                         circular=circular,
                                         verbose=verbose)

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
