import subprocess
from shutil import which
import ostir
import concurrent.futures
import json
import csv
import os

expected_vienna_version = '2.4.17'

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


    results = []
    for RNA in RNAs:
        results.append([RNA, ostir.run_ostir(RNA, posindex=1)])

    result_table = {}
    for result in results:
        result_table[result[0]] = {'findings': result[1]}
    with open(f'{file_location}/ostir_validation_table.json', 'w') as outfile:
        json.dump(result_table, fp=outfile, indent=4)
    print(f'Updated validation table output to {file_location}/ostir_validation_table.json')


if __name__ == '__main__':
    generate_validation_table()
