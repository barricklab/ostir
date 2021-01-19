#!/usr/bin/env python

import argparse
import math
import os
from pathlib import Path
import subprocess
import itertools
import csv
import copy
import re
try:
    from ostir.ostir_factory import OSTIRFactory
except ModuleNotFoundError:
    from ostir_factory import OSTIRFactory

ostir_version = '0.0.2 (In-Development)'

def run_ostir(seq, outfile=None, start_loc=0, end_loc=None, name=None, sd='ACCTCCTTA', threads=1, verbose=False, pos_index=0):
    '''Takes an RNA with optional parameters and returns binding energies.
        Keyword arguments:
        seq -- Sequence to calculate binding energies for
        outfile -- Filepath for output csv
        start_loc -- First base to start considering start codons. Defaults to first base
        end_loc -- Last base to start considering start codons. Defaults to end of sequence
        name -- Returns itself, useful for tagging things for downstream processing
        sd -- Defines anti-Shine-Dalgarno sequence. Defaults to that of E. coli's
        threads -- Defines parallel processing workers, roughly equivalent to multithreading cores
        verbose -- Prints debug information
    '''
    mRNA = seq
    if end_loc == None:
        end_loc = len(mRNA)
    start_range = [start_loc-pos_index, end_loc-pos_index]


    calcObj = OSTIRFactory(mRNA, start_range, sd, verbose=verbose)
    calcObj.calc_dG()
    calcObj.threads = threads

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
        outdata = {'RNA': seq,
                   'codon': output[0][5],
                   'start_pos': output[0][1]+pos_index,
                   'dG_total': output[0][3],
                   'dG_rRNA:mRNA': output[1][0],
                   'dG_mRNA': output[1][2],
                   'dG_Spacing': output[1][4],
                   'Spacing': str(output[1][5]) + ' bp',
                   'dG_Standby': output[1][3],
                   'dG_Start_Codon': output[1][1],
                   'Expression': output[0][0],
                   'name': name
                   }
        output_data_list.append(outdata)
    output_data_list = sorted(output_data_list, key=lambda x: x['start_pos'])

    if outfile:
        save_to_csv(output_data_list, outfile)

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
        if prediction['RNA'] in sorted_predictions.keys():
            sorted_predictions[prediction['RNA']].append(prediction)
        else:
            sorted_predictions[prediction['RNA']] = [prediction]
            out_names = []
            out_names.append(prediction['RNA'])
            if prediction.get('i'):
                out_names.append(prediction.get('i'))
            keys.append(out_names)

    if not keys:
        print('No binding sites were identified.')
        exit(0)

    output_items = ['start_pos', 'codon', 'Expression', 'dG_total', 'dG_rRNA:mRNA', 'dG_mRNA', 'dG_Spacing', 'Spacing', 'dG_Standby', 'dG_Start_Codon']
    row_format = "{:>15}" * (len(output_items))
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

def save_to_csv(outdict, outfile):
    csv_keys = outdict[0].keys()
    with open(outfile, 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, csv_keys)
        dict_writer.writeheader()
        dict_writer.writerows(outdict)


def main():
    '''Main OSTIR function ran when called as an executable. See ostir -h.'''
    parser = argparse.ArgumentParser(description='Open Source Transcription Initiation Rates')

    parser.add_argument(
        '-i', '--input',
        action='store',
        metavar='str/filepath',
        dest='i',
        required=False,
        type=str,
        help="input DNA/RNA. Required if not using --validate.",
    )

    parser.add_argument(
        '-o', '--output',
        action='store',
        metavar= 'filepath',
        dest='o',
        required=False,
        type=str,
        help="Output filepath. If not provided, results will output to the console.",
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
        help="Defines most 5' position to consider start codons.",
    )

    parser.add_argument(
        '-e', '--end',
        action='store',
        metavar='int',
        dest='e',
        required=False,
        type=int,
        help="Defines most 3' position to consider start codons.",
    )

    parser.add_argument(
        '-r', '--rRNA',
        action='store',
        metavar='str',
        dest='r',
        required=False,
        type=int,
        help="Defines rRNA anti-Shine-Dalgarno sequence. Defaults to that of E. coli's.",
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
        metavar='[seq|csv|fasta]',
        dest='t',
        required=False,
        type=str,
        help="Input filetype",
    )

    parser.add_argument(
        '--validate',
        action='store_true',
        dest='validate',
        required=False,
        help="Runs a consistency test to ensure proper installation.",
    )

    options = parser.parse_args()

    if options.validate:
        from ostir.ostir_validation import verify_ostir_install
        verify_ostir_install()
        exit(0)
    elif not options.i:
        subprocess.run(['ostir', '-h'])
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
    if options.v:
        cmd_kwargs['verbose'] = options.v
    if options.s:
        cmd_kwargs['start'] = options.s
    if options.e:
        cmd_kwargs['end'] = options.e
    if options.r:
        cmd_kwargs['sd'] = options.r

    vienna_version = subprocess.check_output(['RNAfold', '--version'])
    vienna_version = str(vienna_version.strip()).replace("'", "").split(' ')[1]
    print(f'Running osTIR Version {ostir_version} (with Vienna version: {vienna_version})')

    # Output data: RNA, Codon, position, dg_total, dg rRNA:mRNA, dg mRNA, dG Spacing, dg Standby, Kinetic Score

    # Determine file input type
    input_type = None
    valid_string_check = re.compile('[ATGCU.-]', re.IGNORECASE)
    if options.t:
        if options.t == 'fasta':
            input_type = 'fasta'
        elif options.t == 'csv':
            input_type = 'csv'
        elif options.t == 'seq':
            input_type = 'string'
        else:
            print(f'Unsupported file type {options.t}.')
            exit(1)
    elif os.path.isfile(cmd_kwargs['seq']) and not input_type:
        filepath_test = Path(cmd_kwargs['seq']).suffix
        if filepath_test == '.fasta':
            input_type = 'fasta'
        elif filepath_test == '.csv':
            input_type = 'csv'
        else:
            with open(cmd_kwargs['seq'], 'r') as in_file:
                first_line = in_file.readline()
                if first_line[0] == '>':
                    input_type = 'fasta'
                elif 'seq' in first_line[0]:
                    input_type = 'csv'
    elif valid_string_check.match(cmd_kwargs['seq']):
        input_type = 'string'

    if input_type == None:
        print(f'Unable to identify the input type. Please define this using the flag "-t"')
        exit(1)


    # Run OSTIR
    if input_type == 'fasta':
        sequences = parse_fasta(cmd_kwargs['seq'])
        result = []
        for sequence in sequences:
            cmd_kwargs['seq'] = sequence[1]
            if 'start' in cmd_kwargs.keys():
                if cmd_kwargs['start']:
                    start_loc = int(cmd_kwargs['start'])
                else:
                    start_loc = 1
            else:
                start_loc = 1
            if 'end' in cmd_kwargs.keys():
                if cmd_kwargs['end']:
                    end_loc = int(cmd_kwargs['end'])
                else:
                    end_loc = start_loc+1
            elif 'start' in cmd_kwargs.keys():
                end_loc = start_loc+1
            else:
                end_loc = len(cmd_kwargs['seq'])
            name = sequence[0]
            verbose = False
            if 'sd' in cmd_kwargs.keys():
                sd = cmd_kwargs['sd']
            else:
                sd = None
            output_dict = run_ostir(cmd_kwargs['seq'], outfile, start_loc,
                                    end_loc, name, sd, threads, verbose, pos_index=1)
            result.append(output_dict)
        result = list(itertools.chain.from_iterable(result))

        if outfile:
            save_to_csv(result, outfile)
        else:
            _print_output(result)

    elif input_type == 'csv':
        csv_keys = []
        csv_values = []
        with open(cmd_kwargs['seq'], 'r') as csv_input:
            reader = csv.reader(csv_input, delimiter=',')
            for i1, row in enumerate(reader):
                if i1 == 0:
                    csv_keys = row
                else:
                    for i2, item in enumerate(row):
                        if item == '':
                            item = None
                        elif '\\ufeff' in item:
                            item = item.replace('\\ufeff', '')
                        if csv_keys[i2] == 'seq':
                            item = item.replace(' ', '')
                        if item:
                            cmd_kwargs[csv_keys[i2]] = item
                    csv_values.append(copy.deepcopy(cmd_kwargs))

        results = []
        for csv_input in csv_values:
            sequence = csv_input['seq']
            if 'start' in csv_input.keys():
                if csv_input['start']:
                    start_loc = int(csv_input['start'])
                else:
                    start_loc = 1
            else:
                start_loc = 1
            if 'end' in csv_input.keys():
                if csv_input['end']:
                    end_loc = int(csv_input['end'])
                else:
                    end_loc = start_loc+1
            elif 'start' in csv_input.keys():
                end_loc = start_loc+1
            else:
                end_loc = len(sequence)
            if csv_input.get('name'):
                name = csv_input.get('name')
            else:
                name = None
            verbose = False
            if 'sd' in csv_input.keys():
                sd = csv_input['sd']
            else:
                sd = None
            output_dict = run_ostir(sequence, outfile, start_loc,
                                    end_loc, name,sd, threads, verbose, pos_index=1)

            results.append(output_dict)

        results = list(itertools.chain.from_iterable(results))
        if outfile:
            save_to_csv(results, outfile)
        else:
            _print_output(results)

    elif input_type == 'string':
        if 'start' in cmd_kwargs.keys():
            if cmd_kwargs['start']:
                start_loc = int(cmd_kwargs['start'])
            else:
                start_loc = 1
        else:
            start_loc = 1
        if 'end' in cmd_kwargs.keys():
            if cmd_kwargs['end']:
                end_loc = int(cmd_kwargs['end'])
            else:
                end_loc = start_loc+1
        elif 'start' in cmd_kwargs.keys():
            end_loc = start_loc+1
        else:
            end_loc = len(cmd_kwargs['seq'])
        name = None
        verbose = False
        if 'sd' in cmd_kwargs.keys():
            sd = cmd_kwargs['sd']
        else:
            sd = None

        output_dict = run_ostir(cmd_kwargs['seq'], outfile, start_loc,
                                end_loc, name, sd, threads, verbose, pos_index=1)
        if outfile:
            save_to_csv(output_dict, outfile)
        else:
            _print_output(output_dict)


if __name__ == "__main__":
    main()
