"""
Shortcuts for common tasks
"""

import sys
import csv
from itertools import groupby
from .data_classes import ostir_task
from .ostir_farm import ostir_farm

# ----- Process CSV's ----

def from_csv(csv_file,
             threads=1,
             verbosity=0,
             **kwargs):
    """Parses a given CSV and runs OSTIR on it. Returns a dictionary of results. Additional kwargs are used as defaults."""
    asd = kwargs.get('asd', 'ACCTCCTTA') or 'ACCTCCTTA'
    start = kwargs.get('start', 0) or 0
    end = kwargs.get('end', None)  # None means the end of the sequence
    circular = kwargs.get('circular', False)
    constraints = kwargs.get('constraints', None)
    name = kwargs.get('name', None) # Prefex for unnamed sequences

    tasks = []

    with open(csv_file, 'r', encoding='UTF-8-sig') as f:
        reader = csv.DictReader(BlankCommentCSVFile(f))
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
                error_string = f"Required column 'sequence' or 'seq' not found for CSV file row: {row}"
                print(error_string, file=sys.stderr)
                sys.exit(1)

            sequence = sequence.replace(' ', '') #remove spaces

            # Assign a name if one is not given from name/id columns
            # If empty assign one based on the index
            name = row.get('name')
            if not name:
                name = row.get('id') or "sequence_" + str(on_seq_index)

            aSD = row.get('anti-shine-dalgarno', asd) or asd #remember, keys lowercased here

            end_loc_1 = row.get('end', end) or end or row.get('start', start) or len(sequence)
            start_loc_1 = row.get('start', start) or start or 1
            circular = row.get('circular', circular) or circular
            constraints = row.get('constraints', constraints) or constraints

            if not circular:
                pass
            elif circular.lower() in ['true', 't', 'yes', 'y', '1']:
                circular = True
            elif circular.lower() in ['false', 'f', 'no', 'n', '0']:
                circular = False
            else:
                print(f"Invalid value for circular: {circular}", file=sys.stderr)
                sys.exit(1)

            # Make an ostir_task from the collected data
            task = ostir_task(sequence,
                              name=name,
                              rrna=aSD,
                              start=int(start_loc_1),
                              end=int(end_loc_1),
                              circular=circular,
                              constraints=constraints)
            tasks.append(task)

        # If any of these are not defined or are empty, use command line values as the defaults

    results = ostir_farm(tasks, threads=threads, verbosity=verbosity)
    # group results by the first value in the list and sort by start position
    results = {k: list(v) for k, v in groupby(results, lambda x: x[0])}
    for key, value in results.items():
        results[key] = sorted([data for _, data in value], key=lambda x: x.start_position)
        

    # Insure the dictionary is sorted by key
    results = {k: results[k] for k in sorted(results)}

    return results

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

# ----- Process FASTA's ----

def from_fasta():
    pass

# ----- Process strings ----

def from_string():
    pass

# ----- Process GenBank's ----

def from_genbank():
    pass



# ---- Simplify outputs ----
def simplify_results(results):
    """Takes a dictionary of results and returns a simplified dictionary of results"""
    simplified_results = {}
    for key, value in results.items():
        simplified_results[key] = value.simplify()
    return simplified_results

