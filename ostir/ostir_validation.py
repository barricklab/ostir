import subprocess
from ostir import ostir
from shutil import which
import concurrent.futures
import json
import csv
import os



expected_vienna_version = '2.4.17'

class ValidationError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


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
        raise UserWarning(f'Vienna version did not match {expected_vienna_version} (was {vienna_version}, '
                          f'this may invalidate test results).')

    # Check to make sure predefined RNAs match expected
    with open(f'{file_location}/ostir_validation_table.json', 'r') as infile:
        validation_table = json.load(infile)
    multiprocessing_in = [[index, value['details'], value['findings']] for index, value in validation_table.items()]

    def _verify_predefined_RNAs(table_info):
        findings = ostir.run_ostir(table_info[0], detailed_out=True)
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
