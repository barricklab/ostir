#!/usr/bin/env python

import unittest
import os
import subprocess
import csv
import shutil

import ostir.ostir as ostir

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
output_directory_path = os.path.join(THIS_DIR, 'output')


def csv_are_identical(csv_file_path_1, csv_file_path_2):
    'Utility function for comparing output, Returns whether files are identical as bool'
    failed = False
    try:
        with open(csv_file_path_1) as csv1:
            reader = csv.reader(csv1)
            csv_1_contents = [row for row in reader]
    except:
        print(f'Failed to open {str(csv_file_path_1)}')
        failed = True

    try:
        with open(csv_file_path_2) as csv2:
            reader = csv.reader(csv2)
            csv_2_contents = [row for row in reader]
    except:
        print(f'Failed to open {str(csv_file_path_2)}')
        failed = True

    if len(csv_1_contents) != len(csv_2_contents):
        print(f'CSVs are different lengths')
        failed = True

    differences = []
    for i in range(0, len(csv_1_contents)):
        if csv_1_contents[i] != csv_2_contents[i]:
            # print(csv_1_rows[i], "\n!=\n",csv_1_rows[i])
            differences.append(i)
    if differences:
        print(f'CSVs have {len(differences)} different values at rows {differences}')
        failed = True

    if failed == True:
        return False
    else:
        return True


##############################################################################################
# Unit tests (function calls)
##############################################################################################
class test_unit_run_ostir_one_sequence(unittest.TestCase):
    def test_unit_run_ostir_one_sequence(self):
        """
        Test run_ostir on one sequence
        """

        print("\n" + "Unit test: run_ostir_one_sequence")

        test_result = ostir.run_ostir("ACUUCUAAUUUAUUCUAUUUAUUCGCGGAUAUGCAUAGGAGUGCUUCGAUGUCAU")

        # print(test_result)

        expected_result = [
            {
                'name': 'unnamed',
                'start_codon': 'AUG',
                'start_position': 31,
                'expression': 2.1121,
                'RBS_distance_bp': 16,
                'dG_total': 16.3277,
                'dG_rRNA:mRNA': -0.481,
                'dG_mRNA': -7.2,
                'dG_spacing': 10.8027,
                'dG_standby': 0.0,
                'dG_start_codon': -1.194,
            },
            {
                'name': 'unnamed',
                'start_codon': 'GUG',
                'start_position': 41,
                'expression': 845.1532,
                'RBS_distance_bp': 8,
                'dG_total': 1.3491,
                'dG_rRNA:mRNA': -3.581,
                'dG_mRNA': -3.6,
                'dG_spacing': 1.4049,
                'dG_standby': 0.0,
                'dG_start_codon': -0.0748,
            },
            {
                'name': 'unnamed',
                'start_codon': 'AUG',
                'start_position': 49,
                'expression': 12448.1756,
                'RBS_distance_bp': 5,
                'dG_total': -5.375,
                'dG_rRNA:mRNA': -7.981,
                'dG_mRNA': -3.6,
                'dG_spacing': 0.0,
                'dG_standby': -0.2,
                'dG_start_codon': -1.194,
            }
        ]

        self.assertEqual(test_result, expected_result)


class test_unit_run_ostir_one_sequence_new_aSD(unittest.TestCase):
    def test_unit_run_ostir_one_sequence_new_aSD(self):
        """
        Test run_ostir on one sequence
        """

        print("\n" + "Unit test: run_ostir_one_sequence_new_aSD")

        test_result = ostir.run_ostir("ACUUCUAAUUUAUUCUAUUUAUUCGCGGAUAUGCAUAGGAGUGCUUCGAUGUCAU", start=31,
                                      aSD="ACGTCCCTA")

        # print(test_result)

        expected_result = [
            {
                'name': 'unnamed',
                'start_codon': 'AUG',
                'start_position': 31,
                'expression': 367.8105,
                'RBS_distance_bp': 4,
                'dG_total': 3.4289,
                'dG_rRNA:mRNA': -2.581,
                'dG_mRNA': -7.2,
                'dG_spacing': 0.0039,
                'dG_standby': 0.0,
                'dG_start_codon': -1.194
            },
        ]

        # print(list(test_result))
        # print(list(expected_result))
        self.assertEqual(test_result, expected_result)

class test_unit_constraints(unittest.TestCase):
    def test_unit_constraints(self):
        """
        Tests running with arbitrary constraints
        """
        print("\n" + "Unit test: constraints")
        test_result = ostir.run_ostir("ACUUCUAAUUUAUUCUAUUUAUUCGCGGAUAUGCAUAGGAGUGCUUCGAUGUCAU", start=31,
                                      aSD="ACGTCCCTA", constraints=".....xxxx.")
        expected_result = [{'name': 'unnamed',
                            'start_codon': 'AUG',
                            'start_position': 31,
                            'expression': 886.8033,
                            'RBS_distance_bp': 4,
                            'dG_total': 1.2289,
                            'dG_rRNA:mRNA': -2.581,
                            'dG_mRNA': -5.0,
                            'dG_spacing': 0.0039,
                            'dG_standby': 0.0,
                            'dG_start_codon': -1.194}]

        self.assertEqual(test_result, expected_result)


##############################################################################################
# Integration tests (command line calls)
##############################################################################################

def setUpModule():
    'Handles creation of the output directory and removal of any stale files from a previous test run.'
    # Delete existing directory
    try:
        shutil.rmtree(output_directory_path)
        print(f"Previous 'output' directory removed: {output_directory_path}")
    except OSError as e:
        pass

    # Create fresh output directory
    try:
        os.mkdir(output_directory_path)
        print(f"Created 'output' directory: {output_directory_path}")
    except OSError as e:
        print(f"Failed to create 'output' directory: {output_directory_path}")
        exit(1)


class test_integration_command_line_FASTA_input(unittest.TestCase):
    def test_integration_command_line_FASTA_input(self):
        input_path = os.path.join(THIS_DIR, 'input', 'command_line_FASTA_input.fa')
        output_path = os.path.join(THIS_DIR, 'output', 'command_line_FASTA_input.csv')
        expected_path = os.path.join(THIS_DIR, 'expected', 'command_line_FASTA_input.csv')
        the_command = f'ostir -j 4 -i "{input_path}" -o "{output_path}"'
        print("\n" + the_command)
        subprocess.call(the_command, shell=True)
        self.assertEqual(True, csv_are_identical(output_path, expected_path))


class test_integration_command_line_string_input(unittest.TestCase):
    def test_integration_command_line_string_input(self):
        output_path = os.path.join(THIS_DIR, 'output', 'command_line_string_input.csv')
        expected_path = os.path.join(THIS_DIR, 'expected', 'command_line_string_input.csv')
        input_sequence = "TTCTAGAAAAAAAATAAGGAGGTAAAATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATG"
        the_command = f'ostir -j 4 -p -i "{input_sequence}" -o "{output_path}"'
        print("\n" + the_command)
        subprocess.call(the_command, shell=True)
        self.assertEqual(True, csv_are_identical(output_path, expected_path))


class test_integration_command_line_CSV_input(unittest.TestCase):
    def test_integration_command_line_CSV_input(self):
        input_path = os.path.join(THIS_DIR, 'input', 'command_line_CSV_input.csv')
        output_path = os.path.join(THIS_DIR, 'output', 'command_line_CSV_input.csv')
        expected_path = os.path.join(THIS_DIR, 'expected', 'command_line_CSV_input.csv')
        the_command = f'ostir -j 4 -i "{input_path}" -o "{output_path}"'
        print("\n" + the_command)
        subprocess.call(the_command, shell=True)
        self.assertEqual(True, csv_are_identical(output_path, expected_path))


class test__integration_command_line_CSV_input_alternate_columns_and_defaults(unittest.TestCase):
    def test_integration_command_line_CSV_input(self):
        input_path = os.path.join(THIS_DIR, 'input', 'command_line_CSV_input_alternate_columns_and_defaults.csv')
        output_path = os.path.join(THIS_DIR, 'output', 'command_line_CSV_input_alternate_columns_and_defaults.csv')
        expected_path = os.path.join(THIS_DIR, 'expected', 'command_line_CSV_input_alternate_columns_and_defaults.csv')
        the_command = f'ostir -j 4 -a TCTGAAGAC -p -q -i "{input_path}" -o "{output_path}"'
        print("\n" + the_command)
        subprocess.call(the_command, shell=True)
        self.assertEqual(True, csv_are_identical(output_path, expected_path))


class test_integration_Salis2009(unittest.TestCase):
    def test_integration_Salis2009(self):
        """
        Test of running Salis2009 sequences
        """

        input_path = os.path.join(THIS_DIR, 'input', 'Salis2009.csv')
        output_path = os.path.join(THIS_DIR, 'output', 'Salis2009.csv')
        expected_path = os.path.join(THIS_DIR, 'expected', 'Salis2009.csv')
        the_command = f'ostir -j 8 -p -i "{input_path}" -o "{output_path}"'
        print("\n" + the_command)
        subprocess.call(the_command, shell=True)
        self.assertEqual(True, csv_are_identical(output_path, expected_path))


class test_integration_T7_genome(unittest.TestCase):
    def test_integration_T7_genome(self):
        """
        Test of running T7 genome sequence
        """

        input_path = os.path.join(THIS_DIR, 'input', 'T7_genome.fasta')
        output_path = os.path.join(THIS_DIR, 'output', 'T7_genome.csv')
        expected_path = os.path.join(THIS_DIR, 'expected', 'T7_genome.csv')
        the_command = f'ostir -j 8 -i "{input_path}" -o "{output_path}"'
        print("\n" + the_command)
        subprocess.call(the_command, shell=True)
        self.assertEqual(True, csv_are_identical(output_path, expected_path))

if __name__ == "__main__":
    unittest.main()


