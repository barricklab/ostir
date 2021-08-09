from Bio import SeqIO
from ostir.ostir import run_ostir, save_to_csv
import sys
from os.path import exists
import os
import csv
from datetime import datetime

def ostir_from_gb(gb):
    genbank = SeqIO.read(gb, "genbank")
    sequence = str(genbank.seq)
    length = len(sequence)
    strand = "-"
    bp_to_iter = 50000
    start_pos = 0
    end_pos = start_pos + bp_to_iter
    print(f'{length} bases to go through')
    while True:
        outpath = f"ecoli/{start_pos}_{end_pos}_{strand}.csv"
        if not exists(outpath):
            print(f'Running ostir on E. coli {start_pos} to {end_pos} on {strand} strand')
            results = run_ostir(in_seq=sequence, start=start_pos, end=end_pos, name=None, aSD=None, threads=8, decimal_places=4,
                                verbose=False)
            column_names = ['name', 'start_codon', 'start_position', 'expression', 'RBS_distance_bp', 'dG_total',
                            'dG_rRNA:mRNA', 'dG_mRNA', 'dG_spacing', 'dG_standby', 'dG_start_codon']
            save_to_csv(column_names, results, outpath)
            print(f'[{datetime.now()}] Results written to {outpath}', file=sys.stderr)

        if end_pos == length:
            if strand == "-":
                break
            else:
                strand = "-"
                start_pos = 0
                end_pos = start_pos + bp_to_iter
                sequence = str(genbank.seq.reverse_complement())
                continue

        start_pos = end_pos+1
        end_pos = end_pos + bp_to_iter if (end_pos + bp_to_iter) <= length else length
    print('Finished OSTIR on the genome')

def start_codon_csv(gb):
    genbank = SeqIO.read(gb, "genbank")

    forward_CDS = {}
    reverse_CDS = {}
    forward_starts_ostir = []
    reverse_starts_ostir = []

    print('Building dataframe for start codons')
    print(f'{len(os.listdir("ecoli"))} files to go through')


    for filenumber, file in enumerate(os.listdir('ecoli')):  # Get start positions from CSVs
        print(f'Fetching start for file {filenumber+1} of {len(os.listdir("ecoli"))}')
        with open(f'ecoli/{file}') as infile:
            csv_file = csv.reader(infile)
            if '+' in file:
                target_list = forward_starts_ostir
                strand = '+'
            elif '-' in file:
                target_list = reverse_starts_ostir
                strand = '-'
            else:
                continue
            for i, row in enumerate(csv_file):
                if i == 0:
                    continue
                data = [int(row[2]),  # Start position
                        float(row[3]),  # Expression estimate
                        strand]
                target_list.append(data)

    gb_length = len(str(genbank.seq))

    reverse_starts_ostir = [[gb_length-start_position[0], start_position[1], start_position[2]] for
                            start_position in reverse_starts_ostir]

    print("Parsing genbank file")
    for record in SeqIO.parse(gb, 'genbank'):  # Get CDS positions from gb file
        for feature in record.features:
            if feature.type == "CDS":
                data = [int(feature.location.start), int(feature.location.end), feature.qualifiers["gene"]]
                if feature.strand == 1:
                    target_dict = forward_CDS
                    target_dict[str(int(feature.location.start))] = data
                elif feature.strand == -1:
                    target_dict = reverse_CDS
                    target_dict[str(int(feature.location.end))] = data
                else:
                    raise ValueError('Found tripple stranded DNA?')


    length_forward = len(forward_starts_ostir)
    total_length = length_forward + len(reverse_starts_ostir)
    print(f'{total_length} positions to check for CDS')

    last_time = datetime.now()
    matches = 0

    csv_out = open('compiled_data.csv', 'w')
    csv_out.write('start_position,expression,is_annotated_start,strand,offset_base,offset_codon\n')

    def write_csv(content):
        csv_out.write(f'{content["start_position"]},{content["expression"]},{content["feature"]},{content["strand"]},{content["bases_off"]},{content["codons_off"]}\n')
        return True

    for i, start in enumerate(forward_starts_ostir):
        if (i+1) % 5000 == 0:
            time_difference = datetime.now() - last_time
            print(f"[{datetime.now()}] Checked {i+1} of {total_length} positions ({matches} matches, time delta {time_difference}))")
            last_time = datetime.now()
        found_match = False
        potential_matches = []
        for n in range(-30, 31):
            if n != 0 and n % 3 != 0:  # Reject out of frame annotations
                continue
            else:
                if str(start[0] + n - 1) in forward_CDS.keys():
                    potential_matches.append(n)
        if potential_matches:
            n = min(potential_matches, key=abs)  # Finds the closest match
            write_csv({
                        "start_position": start[0] - 1,  # Start position of the OSTIR predicted CDS, NOT Genbank
                        "expression": start[1],  # Relative expression
                        "feature": forward_CDS[str(start[0] + n - 1)][2],  # Feature ID
                        "strand": "+",
                        "bases_off": n,
                        "codons_off": int(n/3),
                    })
            matches += 1
        else:
            write_csv({
                "start_position": start[0] - 1,  # Start position
                "expression": start[1],  # Relative expression
                "feature": "NA",  # Feature ID
                 "strand": "+",
                "bases_off": "NA",
                "codons_off": "NA",
            })
    for i, start in enumerate(reverse_starts_ostir):
        if (i+1+length_forward) % 5000 == 0:
            time_difference = datetime.now() - last_time
            print(f"[{datetime.now()}] Checked {i+1+length_forward} of {total_length} positions ({matches} matches, time delta {time_difference})")
            last_time = datetime.now()
        potential_matches = []
        for n in range(-30, 31):
            if n != 0 and n % 3 != 0:  # Reject out of frame annotations
                continue
            else:
                if str(start[0] + n + 1) in reverse_CDS.keys():
                    potential_matches.append(n)
        if potential_matches:
            n = min(potential_matches, key=abs)  # Finds the closest match
            write_csv({
                    "start_position": start[0] + 1,  # Start position of the OSTIR predicted CDS, NOT Genbank
                    "expression": start[1],  # Relative expression
                    "feature": reverse_CDS[str(start[0]+n+1)][2],  # Feature ID
                    "strand": "-",
                    "bases_off": -n,   # Considers strandedness -> More negative = more towards 5' end from RBS
                    "codons_off": int(n/3),
                })
            matches += 1
        else:
            write_csv({
                "start_position": start[0] + 1,  # Start position
                "expression": start[1],  # Relative expression
                "feature": "NA",  # Feature ID
                "strand": "-",
                "bases_off": "NA",
                "codons_off": "NA",
            })

    csv_out.close()



if __name__ == "__main__":
    filepath = "NC_000913.3.gb"
    ostir_from_gb(filepath)
    start_codon_csv(filepath)
