import subprocess
from ostir import ostir
from shutil import which
import concurrent.futures
import json
import csv
import os



expected_vienna_version = '2.4.15'

class ValidationError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def generate_validation_table():
    '''Runs osTIR with a hardcoded input file to generate the consistency test standard.'''
    RNAs = []
    starts = []
    csv_keys = []
    file_location = os.path.dirname(os.path.realpath(__file__))
    with open(f'{file_location}/salis_2009_input.csv', 'r') as csv_input:
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
                        RNAs.append(item)
                    elif csv_keys[i2] == 'start':
                        if item == '':
                            item = None
                        else:
                            item = int(item)
                        starts.append(item)

    input_sequences = list(zip(RNAs, starts))



    def _get_rbs_details(data_in):
        '''Internal. Runs vienna on each '''
        seq, start = data_in
        findings, details = ostir.run_ostir(seq, start_loc=start, detailed_out=True)
        return seq, list(findings), details

    with concurrent.futures.ThreadPoolExecutor(1) as multiprocessor:
        results = multiprocessor.map(_get_rbs_details, input_sequences)

    results = [x for x in results]

    result_table = {}
    for result in results:
        result_table[result[0]] = {'findings': result[1],
                                   'details': result[2]}
    with open(f'{file_location}/ostir_validation_table.json', 'w') as outfile:
        json.dump(result_table, fp=outfile, indent=4)



def verify_ostir_install():
    file_location = os.path.dirname(os.path.realpath(__file__))

    # Check for missing critical ViennaRNA executables
    dependencies = [which('RNAfold') is not None,
                    which('RNAsubopt') is not None,
                    which('RNAeval') is not None]
    if False in dependencies:
        raise ValidationError(f'Vienna installation was missing critical programs. Is ViennaRNA installed?')

    #  Check ViennaRNA version
    vienna_version = subprocess.check_output(['RNAfold', '--version'])
    vienna_version = str(vienna_version.strip()).replace("'", "").split(' ')[1]
    if vienna_version != expected_vienna_version:
        raise ValidationError(f'Vienna version did not match {expected_vienna_version} (was {vienna_version}).')

    # Check to make sure predefined RNAs match expected
    with open(f'{file_location}/ostir_validation_table.json', 'r') as infile:
        validation_table = json.load(infile)
    multiprocessing_in = [[index, value['details'], value['findings']] for index, value in validation_table.items()]

    def _verify_predefined_RNAs(table_info):
        findings, details = ostir.run_ostir(table_info[0], detailed_out=True)
        if table_info[2] == details and table_info[1] == findings:
            return False
        else:
            return True

    with concurrent.futures.ThreadPoolExecutor() as multiprocessor:
        result = multiprocessor.map(_verify_predefined_RNAs, multiprocessing_in)
    result = [x for x in result]

    if False in result:  # If true, validation has failed
        raise ValidationError(f'At least 1 RNA test did not match expected.')

    print('OSTIR validation successful')

if __name__ == '__main__':
    verify_ostir_install()
