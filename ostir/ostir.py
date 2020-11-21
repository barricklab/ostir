#!/usr/bin/env python

import argparse
import math
from ostir.ostir_factory import OSTIRFactory
import os
from pathlib import Path
import subprocess
import concurrent.futures
import itertools
import csv
import copy

ostir_version = 1.0

def calc_dG_from_file(handle, output, verbose=True, parameters={}):
    from Bio import SeqIO

    records = SeqIO.parse(handle, "fasta")
    First = True
    export_PDF = True

    # for i in range(30):
    #    records.next()

    for record in records:

        mRNA = record.seq.tostring().upper()

        # Set any defaults
        start_range = [0, len(mRNA)]
        name = record.description.split(" ")[0]

        # Create instance of RBS Calculator
        test = OSTIRFactory(mRNA, start_range, name)

        # Examine kvars dictionary and pull out any options. Assign them to instanced class.
        for (key, value) in list(parameters.items()):

            if key == "cutoff":
                test.cutoff = value
            elif key == "start_range":
                test.start_range = value
            elif key == "rRNA":
                test.rRNA = value
            elif key == "energy_cutoff":
                test.energy_cutoff = value
            elif key == "standby_site_length":
                test.standby_site_length = value
            elif key == "dangles":
                test.dangles = value
            elif key == "export_PDF":
                export_PDF = value

        test.calc_dG()
        test.print_dG(test.infinity, print_expression=verbose)
        test.save_data(output, First)

        if First:
            First = False

        if export_PDF:
            num_structs = len(test.mRNA_rRNA_uncorrected_structure_list)
            for (structure, counter) in zip(test.mRNA_rRNA_uncorrected_structure_list, list(range(num_structs))):
                index = structure["MinStructureID"]
                structure.export_PDF(index, name, filename=name + "_rRNA" + "_" + str(counter) + ".pdf",
                                     program="subopt")

            num_structs = len(test.mRNA_structure_list)
            for (structure, counter) in zip(test.mRNA_structure_list, list(range(num_structs))):
                structure.export_PDF(0, name, filename=name + "_mRNA" + "_" + str(counter) + ".pdf")

    output.close()


def calc_dG_pre_post_RBS(pre_list, post_list, RBS_list, name_list, output, verbose=True, parameters={}):

    First = True

    for (pre, post, RBS, name) in zip(pre_list, post_list, RBS_list, name_list):

        mRNA = pre + RBS + post

        start_range = [0, len(mRNA)]

        # Create instance of RBS Calculator
        test = OSTIRFactory(mRNA, start_range, name)

        # Examine kvars dictionary and pull out any options. Assign them to instanced class.
        for (key, value) in list(parameters.items()):

            if key == "cutoff":
                test.cutoff = value
            elif key == "start_range":
                test.start_range = value
            elif key == "rRNA":
                test.rRNA = value
            elif key == "energy_cutoff":
                test.energy_cutoff = value
            elif key == "standby_site_length":
                test.standby_sitRBSe_length = value
            elif key == "dangles":
                test.dangles = value
            elif key == "export_PDF":
                export_PDF = value

        test.calc_dG()
        if verbose: test.print_dG(test.infinity)

        test.save_data(output, First)
        if First:
            First = False

    output.close()


def run_ostir(seq, constraint_str=None, outfile=None, start_loc=0, end_loc=None, i=None, verbose=False,
                          detailed_out=False, print_out=False, sd=None):
    mRNA = seq
    if end_loc == None:
        end_loc = len(mRNA)
    start_range = [start_loc, end_loc]

    if not sd:
        sd = 'ACCTCCTTA'

    calcObj = OSTIRFactory(mRNA, start_range, sd, constraint_str, verbose=verbose)
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

    if outfile:
        with open(outfile, 'w') as out_file:
            for (expr, start_pos, ks, dG, dstandby, codon) in zip(expr_list, start_pos_list, kinetic_score_list, dG_total_list,
                                                           standby_site_list, start_codon_list):
                out_file.writelines(['1\n', f"{start_pos} {expr} {ks}, {i}\n"])
        if detailed_out:
            second_out = outfile + '.detailed'
            with open(second_out, 'w') as outfile:
                pass
    if print_out:
        calcObj.print_dG()
    if detailed_out:
        return return_var, dG_details
    else:
        return return_var

def parse_fasta(filepath):
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

def _parallelizer(input):

    def _valid_entry(key_to_test):
        if key_to_test not in input.keys():
            return False
        elif input[key_to_test] == '':
            return False
        elif not input[key_to_test]:
            return False
        else:
            return True

    #unpack args:
    seq = input['seq']
    if 'constraint_str' in input.keys():
        constraint_str = input['constraint_str']
    else:
        constraint_str = None
    if 'outfile' in input.keys():
        outfile = input['outfile']
    else:
        outfile = None
    if _valid_entry('start'):
        start_loc = int(input['start'])-1
    else:
        start_loc = 0
    if _valid_entry('end'):
        end_loc = int(input['end'])-1
    elif _valid_entry('start'):
        end_loc = int(input['start'])+1
    else:
        end_loc = len(seq)
    if 'i' in input.keys():
        i = input['i']
    else:
        i = None
    if 'verbose' in input.keys():
        verbose = input['verbose']
    else:
        verbose = None
    if 'detailed_out' in input.keys():
        detailed_out = input['detailed_out']
    else:
        detailed_out = None
    if 'print_out' in input.keys():
        print_out = input['print_out']
    else:
        print_out = False
    if 'sd' in input.keys():
        sd = input['sd']
    else:
        sd = None


    normal, detailed= ostir(seq, constraint_str, outfile, start_loc, end_loc, i, verbose, detailed_out, print_out, sd)
    normal = list(normal)
    zip_output = zip(normal, detailed)
    output_data_list = []
    for output in list(zip_output):
        outdata = {'RNA': seq,
                   'codon': output[0][5],
                   'start_pos': output[0][1]+1,
                   'dG_total': output[0][3],
                   'dG_rRNA:mRNA': output[1][0],
                   'dG_mRNA': output[1][2],
                   'dG_Spacing': output[1][4],
                   'dG_Standby': output[1][3],
                   'dG_Start_Codon': output[1][1],
                   'Expression': output[0][0]
                   }
        output_data_list.append(outdata)
    return output_data_list

def _print_output(outlist):
    sorted_predictions = {}
    keys = []
    for prediction in outlist:
        if prediction['RNA'] in sorted_predictions.keys():
            sorted_predictions[prediction['RNA']].append(prediction)
        else:
            sorted_predictions[prediction['RNA']] = [prediction]
            keys.append(prediction['RNA'])

    output_items = ['start_pos', 'codon', 'Expression', 'dG_total', 'dG_rRNA:mRNA', 'dG_mRNA', 'dG_Spacing', 'dG_Standby', 'dG_Start_Codon']
    row_format = "{:>15}" * (len(output_items))
    print('_________________________________________________')
    for rna in keys:
        print(f'Tested Sequence: {rna}')
        print(row_format.format(*output_items))
        for start in sorted_predictions[rna]:
            output_data = [start[key] for key in output_items]
            for i, data_point in enumerate(output_data):
                if type(data_point) == float:
                    data_point = format(data_point, '.4f')
                    output_data[i] = data_point
            print(row_format.format(*output_data))
        print('_________________________________________________')
        pass


def main():
    parser = argparse.ArgumentParser(description='Open Source Transcription Initiation Rates')

    parser.add_argument(
        '-i', '--input',
        action='store',
        dest='i',
        required=False,
        type=str,
        help="input DNA/RNA. Usage: -i [sequence]. Required if not using --validate.",
    )

    parser.add_argument(
        '-o', '--output',
        action='store',
        dest='o',
        required=False,
        type=str,
        help="Output filepath. If not provided, results will output to the console. Usage: -o [filepath]",
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        dest='v',
        required=False,
        help="Prints verbose output. Usage: -v",
    )

    parser.add_argument(
        '-s', '--start',
        action='store',
        dest='s',
        required=False,
        type=int,
        help="Defines most 5' position to consider start codons Usage: -s [int]",
    )

    parser.add_argument(
        '-e', '--end',
        action='store',
        dest='e',
        required=False,
        type=int,
        help="Defines most 3' position to consider start codons Usage: -e [int]",
    )

    parser.add_argument(
        '-r', '--rRNA',
        action='store',
        dest='r',
        required=False,
        type=int,
        help="Defines rRNA anti-Shine-Dalgarno sequence. Defaults to that of E. coli's. Usage: -r [sequence]",
    )

    parser.add_argument(
        '-j', '--threads',
        action='store',
        dest='j',
        required=False,
        type=int,
        help="Number of threads for multiprocessing",
    )

    parser.add_argument(
        '--validate',
        action='store_true',
        dest='validate',
        required=False,
        help="Runs a consistency test to ensure proper installation. Usage: --validate",
    )

    options = parser.parse_args()

    if options.validate:
        from ostir.ostir_validation import verify_ostir_install
        verify_ostir_install()
        exit()
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
        cores = options.j
    else:
        cores = 1
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

    if os.path.isfile(cmd_kwargs['seq']):
        if Path(cmd_kwargs['seq']).suffix == '.fasta':
            sequences = parse_fasta(cmd_kwargs['seq'])
            concurrent_futures_args = []
            for sequence in sequences:
                cmd_kwargs['seq'] = sequence[1]
                concurrent_futures_args.append(copy.deepcopy(cmd_kwargs))
            with concurrent.futures.ThreadPoolExecutor(max_workers=cores) as multiprocessor:
                result = multiprocessor.map(_parallelizer, concurrent_futures_args)
            result = [x for x in result]
            result = list(itertools.chain.from_iterable(result))
            if outfile:
                csv_keys = result[0].keys()
                with open(outfile, 'w')  as output_file:
                    dict_writer = csv.DictWriter(output_file, csv_keys)
                    dict_writer.writeheader()
                    dict_writer.writerows(result)
            else:
                _print_output(result)
        elif Path(cmd_kwargs['seq']).suffix == '.csv':
            csv_keys = []
            csv_values = []
            with open(cmd_kwargs['seq'], 'r') as csv_input:
                reader = csv.reader(csv_input, delimiter=',')
                for i1, row in enumerate(reader):
                    if i1 == 0:
                        csv_keys = row
                    else:
                        single_input = {}
                        for i2, item in enumerate(row):
                            if item == '':
                                item = None
                            if csv_keys[i2] == 'seq':
                                item = item.replace(' ', '')
                            cmd_kwargs[csv_keys[i2]] = item
                        csv_values.append(copy.deepcopy(cmd_kwargs))
            with concurrent.futures.ThreadPoolExecutor(max_workers=cores) as multiprocessor:
                result = multiprocessor.map(_parallelizer, csv_values)
            result = [x for x in result]
            result = list(itertools.chain.from_iterable(result))
            if outfile:
                csv_keys = result[0].keys()
                with open(outfile, 'w')  as output_file:
                    dict_writer = csv.DictWriter(output_file, csv_keys)
                    dict_writer.writeheader()
                    dict_writer.writerows(result)
            else:
                _print_output(result)

        else:
            raise ValueError('Input file is not a supported filetype')
    else:
        result = _parallelizer(**cmd_kwargs)
        if outfile:
            csv_keys = result[0].keys()
            with open(outfile, 'w')  as output_file:
                dict_writer = csv.DictWriter(output_file, csv_keys)
                dict_writer.writeheader()
                dict_writer.writerows(result)
        else:
            _print_output(result)


if __name__ == "__main__":
    main()