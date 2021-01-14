import subprocess
from ostir import ostir
from shutil import which
import concurrent.futures
import json
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
    RNA_checks = [[index, value['findings']] for index, value in validation_table.items()]

    results = []
    for RNA, finding in RNA_checks:
        return_results = ostir.run_ostir(RNA)
        if return_results == finding:
            results.append(1)
        else:
            results.append(0)
        pass

    result = [x for x in results]

    if 0 in results:  # If true, validation has failed
        raise ValidationError(f'At least 1 RNA test did not match expected value.')

    print('OSTIR validation successful')

if __name__ == '__main__':
    verify_ostir_install()
