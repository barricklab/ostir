import subprocess
from RBS_Calculator_Vienna import RBS_Calculator_Vienna
from shutil import which
import concurrent.futures
import json



expected_vienna_version = '2.4.15'

class ValidationError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def generate_validation_table(): #TODO: Set this up with Salis input
    RNAs = [
    ]

    def _get_rbs_details(seq):
        findings, details = RBS_Calculator_Vienna(seq, detailed_out=True)
        return seq, list(findings), details

    with concurrent.futures.ThreadPoolExecutor() as multiprocessor:
        results = multiprocessor.map(_get_rbs_details, RNAs)

    results = [x for x in results]

    result_table = {}
    for result in results:
        result_table[result[0]] = {'findings': result[1],
                                   'details': result[2]}
    with open('vienna_validation_table.json', 'w') as outfile:
        json.dump(result_table, fp=outfile, indent=4)



def verify_ostir_install():
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
    with open('vienna_validation_table.json', 'r') as infile:
        validation_table = json.load(infile)
    multiprocessing_in = [[index, value['details'], value['findings']] for index, value in validation_table.items()]

    def _verify_predefined_RNAs(table_info):
        findings, details = RBS_Calculator_Vienna(table_info[0], detailed_out=True)
        if table_info[2] == details and table_info[1] == findings:
            return False
        else:
            return True

    with concurrent.futures.ThreadPoolExecutor() as multiprocessor:
        result = multiprocessor.map(_verify_predefined_RNAs, multiprocessing_in)
    result = [x for x in result]

    if False in result:  # If true, validation has failed
        raise ValidationError(f'At least 1 RNA test did not match expected.')

    print('RNA Calculator validation successful')



if __name__ == '__main__':
    verify_ostir_install()
