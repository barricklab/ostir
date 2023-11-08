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
import importlib.util
from shutil import which
from warnings import warn
from pathlib import Path
from .shortcuts import from_csv


try:
    from ostir.ostir_factory import OSTIRFactory
except ModuleNotFoundError:
    from .ostir_factory import OSTIRFactory

# Optional dependency
if importlib.util.find_spec("rich") is not None:
    from rich import print as rprint
else:
    rprint = None

OSTIR_VERSION = '1.1.1'
OLDEST_VIENNA = '2.4.18'

# The E. coli sequence
Ecoli_anti_Shine_Dalgarno = 'ACCTCCTTA'

def run_ostir(in_seq,
              start=None,
              end=None,
              name=None,
              aSD=None,
              threads=1,
              decimal_places=4,
              circular=False,
              constraints=None,
              verbosity=0):
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
        verbosity -- Controls verbosity of output. 0 is silent, 1 is normal, 2 is verbose
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
        error_message = f"ERROR: anti-Shine-Dalgarno sequence provided ({aSD}) contains non-nucleotide characters.\n<<<Sequence ({name}) will be skipped.>>>"
        if not rprint:
            print(error_message, file=sys.stderr)
        else:
            rprint(error_message, style="bold red")
        return []

    # Length check
    if len(aSD) != 9:
        error_message = f"ERROR: anti-Shine-Dalgarno sequence provided ({aSD}) is not 9 bases.\n<<<Sequence ({name}) will be skipped.>>>"
        if not rprint:
            print(error_message, file=sys.stderr)
        else:
            rprint(error_message, style="bold red")
        return []

    #Upper case and convert to RNA
    aSD = aSD.upper()
    aSD = aSD.replace("T","U")

    #Clean spaces from sequence... could clean other characters too
    seq = in_seq.replace(" ", "")

    # Nucleotide character check
    if nucleotides.search(seq) is not None:
        error_message = f"ERROR: Input sequence contains non-nucleotide characters.\n<<<Sequence ({name}) will be skipped.>>>"
        if not rprint:
            print(error_message, file=sys.stderr)
        else:
            rprint(error_message, style="bold red")
        return []


    # Start <= end check
    if in_start_loc_1 is not None and in_end_loc_1 is not None and in_end_loc_1<in_start_loc_1:
        error_message = f"ERROR: Start location ({in_start_loc_1}) is not less than end location ({in_end_loc_1}).\n<<<Sequence ({name}) will be skipped.>>>"
        if not rprint:
            print(error_message, file=sys.stderr)
        else:
            rprint(error_message, style="bold red")
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

    calcObj = OSTIRFactory(seq, start_range_1, aSD, constraints, circular=circular, verbosity=verbosity)
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
            if 'sequence' in prediction:
                out_names.append(prediction['sequence'])
            else:
                out_names.append(None)
            if 'anti-Shine-Dalgarno' in prediction:
                out_names.append(prediction['anti-Shine-Dalgarno'])
            else:
                out_names.append(None)
            keys.append(out_names)

    if not keys:
        if not rprint:
            print('No binding sites were identified.')
        else:
            rprint('[bold red]No binding sites were identified.[/bold red]')
        sys.exit(0)
    output_items = ['start_codon', 'start_position', 'expression', 'RBS_distance_bp', 'dG_total', 'dG_rRNA:mRNA', 'dG_mRNA', 'dG_spacing', 'dG_standby', 'dG_start_codon']
    if not rprint:
        row_format = "{:>16}" * (len(output_items))
        print('_________________________________________________')
        for rna in keys:
            print(f'\nSample: {rna[0]}')
            if rna[1]:
                print(f'Tested Sequence: {rna[1]}')
            if rna[2]:
                print(f'Sequence RNA: {rna[2]}')
            print(row_format.format(*output_items))
            for start in sorted_predictions[rna[0]]:
                output_data = [start[key] for key in output_items]
                for i, data_point in enumerate(output_data):
                    if isinstance(data_point, float):
                        data_point = format(data_point, '.4f')
                        output_data[i] = data_point
                print(row_format.format(*output_data))
            print('_________________________________________________')
    else:
        from rich.table import Table
        from rich.console import Console
        console = Console()
        for rna in keys:
            print(f'\nSample: {rna[0]}')
            if rna[1]:
                print(f'Tested Sequence: {rna[1]}')
            if rna[2]:
                print(f'Sequence RNA: {rna[2]}')
            table = Table()
            for item in output_items:
                table.add_column(item)
            for start in sorted_predictions[rna[0]]:
                output_data = [start[key] for key in output_items]
                for i, data_point in enumerate(output_data):
                    if isinstance(data_point, float):
                        data_point = format(data_point, '.4f')
                        output_data[i] = data_point
                table.add_row(str(start['start_codon']),
                              str(start['start_position']),
                                str(start['expression']),
                                str(start['RBS_distance_bp']),
                                str(start['dG_total']),
                                str(start['dG_rRNA:mRNA']),
                                str(start['dG_mRNA']),
                                str(start['dG_spacing']),
                                str(start['dG_standby']),
                                str(start['dG_start_codon']))

            console.print(table)
            
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
        '-v', '--versity',
        action='store',
        metavar= 'int',
        dest='v',
        required=False,
        default="1",
        help="Sets the verbosity level. Default 0 is quiet, 1 is normal, 2 is verbose",
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
        print_string = f'OSTIR version {OSTIR_VERSION}'
        if not rprint:
            print(print_string, file=sys.stdout)
        else:
            rprint(print_string)
        sys.exit(0)

    if not options.i:
        print_string = "Input (-i) required."
        if not rprint:
            print(print_string, file=sys.stderr)
        else:
            rprint(print_string, file=sys.stderr)
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
        cmd_kwargs['threads'] = threads
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
    if options.v:
        if options.v not in ["0", "1", "2"]:
            print_string = "Verbosity (-v) must be 0, 1, or 2."
            print(options.v)
            if not rprint:
                print(print_string, file=sys.stderr)
            else:
                rprint(print_string, file=sys.stderr)
            parser.print_help()
            sys.exit(1)
        cmd_kwargs['verbosity'] = int(options.v)
        verbosity = cmd_kwargs['verbosity']

    # Check if viennaRNA is installed
    dependencies = [which('RNAfold') is not None,
                    which('RNAsubopt') is not None,
                    which('RNAeval') is not None]

    if False in dependencies:
        raise EnvironmentError('ViennaRNA is not properly installed or in PATH')

    vienna_version = subprocess.check_output(['RNAfold', '--version'])
    vienna_version = str(vienna_version.strip()).replace("'", "").split(' ')[1]
    if verbosity > 0:
        print_string = f'Running OSTIR version {OSTIR_VERSION} (with Vienna version: {vienna_version})'
        if not rprint:
            print(print_string, file=sys.stdout)
        else:
            rprint(print_string)

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
            error_string = f'Unsupported file type {options.t}.'
            if not rprint:
                print(error_string, file=sys.stderr)
            else:
                rprint(error_string, file=sys.stderr)
            sys.exit(1)
    elif os.path.isfile(options.i) and not input_type:  # Get input type from file name
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
            error_string = f'Unable to identify the type of file specified as input (-i). Please define it using "-t".'
            if not rprint:
                print(error_string, file=sys.stderr)
            else:
                rprint(error_string, file=sys.stderr)
        else:
            error_string = f'Unable to identify the type of input (-i). Please define it using "-t".'
            if not rprint:
                print(error_string, file=sys.stderr)
            else:
                rprint(error_string, file=sys.stderr)
        sys.exit(1)


    # Run OSTIR
    results = []
    ## String input ##############################################################
    if input_type == 'string':
        if verbosity > 0:
            print_string = 'Reading input sequence from command line'
            if not rprint:
                print(print_string, file=sys.stdout)
            else:
                rprint(print_string)

        sequence = cmd_kwargs.get('seq')
        name = None
        start_loc_1 = cmd_kwargs.get('start')
        end_loc_1 = cmd_kwargs.get('end')
        aSD = cmd_kwargs.get('aSD')
        verbosity = cmd_kwargs.get('verbosity')
        circular = cmd_kwargs.get('circular')

        output_dict_list = run_ostir(sequence,
                                     start=start_loc_1,
                                     end=end_loc_1,
                                     name=name,
                                     aSD=aSD,
                                     threads=threads,
                                     circular=circular,
                                     verbosity=verbosity)

        for output_dict in output_dict_list:
            if cmd_kwargs.get('print_mRNA_sequence'):
                output_dict['sequence'] = sequence
            if cmd_kwargs.get('print_aSD_sequence'):
                output_dict['anti-Shine-Dalgarno'] = aSD

        results.extend(output_dict_list)

    ## FASTA input ##############################################################
    elif input_type == 'fasta':
        input_file = cmd_kwargs['seq']
        if verbosity > 0:
            print_string = f'Reading FASTA file {input_file}'
            if not rprint:
                print(print_string, file=sys.stdout)
            else:
                rprint(print_string)
        sequence_entries = parse_fasta(input_file)
        for sequence_entry in sequence_entries:
            sequence = sequence_entry[1]
            name = sequence_entry[0]
            start_loc_1 = cmd_kwargs.get('start')
            end_loc_1 = cmd_kwargs.get('end')
            aSD = cmd_kwargs.get('aSD')
            circular = cmd_kwargs.get('circular')
            verbosity = cmd_kwargs.get('verbosity')

            output_dict_list = run_ostir(sequence,
                                         start=start_loc_1,
                                         end=end_loc_1,
                                         name=name,
                                         aSD=aSD,
                                         threads=threads,
                                         circular=circular,
                                         verbosity=verbosity)

            for output_dict in output_dict_list:
                if cmd_kwargs.get('print_mRNA_sequence'):
                    output_dict['sequence'] = sequence
                if cmd_kwargs.get('print_aSD_sequence'):
                    output_dict['anti-Shine-Dalgarno'] = aSD

            results.extend(output_dict_list)


    ## CSV input ##############################################################
    elif input_type == 'csv':
        input_file = cmd_kwargs['seq']
        if verbosity > 0:
            print_string = f'Reading CSV file {input_file}'
            if not rprint:
                print(print_string, file=sys.stdout)
            else:
                rprint(print_string)

        result = from_csv(csv_file=input_file,
                        threads=cmd_kwargs.get('threads'),
                        asd=cmd_kwargs.get('aSD'),
                        start=cmd_kwargs.get('start'),
                        end=cmd_kwargs.get('end'),
                        verbosity=cmd_kwargs.get('verbosity'),
                        circular=cmd_kwargs.get('circular'),)

        # We also need the metadata from the CSV file
        run_metadata = []
        with open(input_file, 'r', encoding='UTF-8-sig') as csv_file:
            reader = csv.DictReader(BlankCommentCSVFile(csv_file))
            for row in reader:
                run_metadata.append(row)

        results = []
        for key, finding in result.items():
            result_set = []
            for value in finding:
                value_dict = value.results()
                if cmd_kwargs.get('print_mRNA_sequence'):
                    value_dict['sequence'] = run_metadata[key].get('sequence') or run_metadata[key].get('seq')
                if cmd_kwargs.get('print_aSD_sequence'):
                    value_dict['anti-Shine-Dalgarno'] = run_metadata[key].get('anti-Shine-Dalgarno') or run_metadata[key].get('asd') or cmd_kwargs.get('aSD')
                result_set.append(value_dict)
            results.extend(result_set)
        

    ## Output - for all ways of running ##############################################################
    if outfile:
        #if the output is empty, we need to print a valid header...
        column_names = ['name', 'start_codon', 'start_position', 'expression', 'RBS_distance_bp', 'dG_total', 'dG_rRNA:mRNA', 'dG_mRNA', 'dG_spacing', 'dG_standby', 'dG_start_codon']
        if cmd_kwargs.get('print_mRNA_sequence'):
            column_names.insert(1, 'sequence')
        if cmd_kwargs.get('print_aSD_sequence'):
            column_names.insert(1, 'anti-Shine-Dalgarno')
        save_to_csv(column_names, results, outfile)
        if cmd_kwargs['verbosity'] > 0:
            print_string = f'Results written to {outfile}'
        if cmd_kwargs['verbosity'] > 0:
            if not rprint:
                print(print_string, file=sys.stdout)
            else:
                rprint(print_string)
    else:
        _print_output(results)


if __name__ == "__main__":
    main()
