# OSTIR (Open Source Translation Initiation Rates)

`OSTIR` is a
Python package for predicting the rates at which ribosomes will bind to and initiate
translation from different start codons in bacterial mRNAs. It uses the ViennaRNA software
suite to perform the necessary free energy calculations. The code builds on the last open
source version of the
[RBS calculator](https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0)
that implements the calculations described in [Salis 2009](https://doi.org/10.1038/nbt.1568).

`OSTIR` includes several improvements in usability. It supports multi-FASTA
input with command line parameters or CSV input that can define
parameters on a per-sequence basis. Additionally, `OSTIR` supports multi-threaded
execution, accelerating use cases that require the analysis of very large sequences.

## Installation
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ostir/README.html)

From Bioconda:
- Run `conda install -c bioconda ostir`

From Pip:
- Download and install ViennaRNA, following the instructions [here](https://www.tbi.univie.ac.at/RNA/).
- Run `pip install ostir`

From Source:
- Clone this repository
- Navigate to the `./ostir/` directory
- Install/update the necessary packaging tools with `python3 -m pip install --user --upgrade setuptools wheel`
- Package the code with `python3 setup.py sdist bdist_wheel`
- Install the package with `python3 -m pip install ./`
- Download and install ViennaRNA, following the instructions [here](https://www.tbi.univie.ac.at/RNA/).
- To test your install run `python -m unittest`

## Command Line Usage

OSTIR can be executed via the included command-line script, `ostir`.

```
usage: ostir [-h] -i str/filepath [-o filepath] [-v] [-s int] [-e int] [-a str] [-p] [-q] [-j int] [-t [string|csv|fasta]]

Open Source Transcription Initiation Rates

optional arguments:
  -h, --help            show this help message and exit
  -i str/filepath, --input str/filepath
                        Input filename (FASTA/CSV) or DNA/RNA sequence. For CSV input files, there must be a 'seq' or 'sequence'
                        column. Other columns will override any options provided at the command line if they are present:
                        'name/id', 'start', 'end', 'anti-Shine-Dalgarno'.
  -o filepath, --output filepath
                        Output file path. If not provided, results will output to the console.
  -v, --verbose         Prints verbose output.
  -s int, --start int   Most 5' position to consider a start codon beginning.
  -e int, --end int     Most 3' position to consider a start codon beginning
  -a str, --anti-Shine-Dalgarno str
                        anti-Shine-Dalgarno sequence: the 9 bases located at the 3' end of 16S rRNA. May be provided as DNA or
                        RNA. Defaults to that of E. coli (ACCTCCTTA).
  -p, --print-sequence  Include the input mRNA sequence in output CSV files
  -q, --print-anti-Shine-Dalgarno
                        Include the anti-Shine-Dalgarno sequence in output CSV files
  -j int, --threads int
                        Number of threads for multiprocessing
  -t [string|csv|fasta], --type [string|csv|fasta]
                        Input type (overrides autodetection)
```

### Example of Command Line Input/Output

Example command for specifying the mRNA sequence to search as a parameter:
```bash
ostir -i TTCTAGATGAGAATAAGGTTATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGT 
```

Output of this command:
```bash
_________________________________________________
     start_codon  start_position      expression RBS_distance_bp        dG_total    dG_rRNA:mRNA         dG_mRNA      dG_spacing      dG_standby  dG_start_codon
             ATG               7          3.4040              -2         15.1346         -1.9810         -1.1000         17.2096          0.0000         -1.1940
             ATG              21        643.9384               4          2.0289         -5.2810         -8.5000          0.0039          0.0000         -1.1940
             ATG              54          0.0558               1         25.4101         -4.3810        -13.8000         17.1851          0.0000         -1.1940
             ATG              72        227.5882               4          4.6289         -0.6810         -6.5000          0.0039          0.0000         -1.1940
_________________________________________________
```

### Example of CSV Input/Output

Example CSV input file (file: `input.csv`):
| id              | seq                                                                                           | anti-Shine-Dalgarno |
|-----------------|-----------------------------------------------------------------------------------------------|---------------------|
| first_sequence  | TTCTAGActttaatttaacgtaaataaggaagtcattATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGA | ACCTCCTTA           |
| second_sequence | TTCTAGActttaatttaacgtaaataaggaagtcattATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGA |                     |
|                 | TTCTAGActttaatttaacgtaaataaggaagtcattATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGA | CCCCCCCCC           |

This example demonstrates that if a column corresponding to a parameter is missing or blank for one input sequence, then the default value for that parameter will be substituted for the missing value. Additional examples of CSV input files can be found in [tests/input](tests/input). 

Example command that specifies CSV input and output files, assigns `TCTGAAGAC` as the default anti-Shine-Dalgarno sequence, includes the anti-Shine-Dalgarno sequence as a column in the output, and uses four threads for parallelization:
```bash
ostir -j 4 -a TCTGAAGAC -q -i input.csv -o output.csv
```

Example CSV output (file: `output.csv`):
| name            | anti-Shine-Dalgarno | start_codon | start_position | expression | RBS_distance_bp | dG_total | dG_rRNA:mRNA | dG_mRNA | dG_spacing | dG_standby | dG_start_codon |
|-----------------|---------------------|-------------|----------------|------------|-----------------|----------|--------------|---------|------------|------------|----------------|
| first_sequence  | ACCTCCTTA           | ATG         | 38             | 643.9384   | 4               | 2.0289   | -8.381       | -11.6   | 0.0039     | 0.0        | -1.194         |
| first_sequence  | ACCTCCTTA           | ATG         | 71             | 0.0558     | 1               | 25.4101  | -4.381       | -13.8   | 17.1851    | 0.0        | -1.194         |
| first_sequence  | ACCTCCTTA           | ATG         | 89             | 227.5882   | 4               | 4.6289   | -0.681       | -6.5    | 0.0039     | 0.0        | -1.194         |
| second_sequence | TCTGAAGAC           | ATG         | 38             | 34.6638    | 7               | 9.3332   | -1.881       | -11.6   | 0.8082     | 0.0        | -1.194         |
| second_sequence | TCTGAAGAC           | ATG         | 71             | 40.1662    | 6               | 8.9649   | -3.981       | -13.8   | 0.3399     | 0.0        | -1.194         |
| second_sequence | TCTGAAGAC           | ATG         | 89             | 460.8985   | 6               | 2.8649   | -3.381       | -6.5    | 0.3399     | -0.6       | -1.194         |
| sequence_3      | CCCCCCCCC           | ATG         | 38             | 44.2111    | 5               | 8.725    | -1.681       | -11.6   | 0.0        | 0.0        | -1.194         |
| sequence_3      | CCCCCCCCC           | ATG         | 71             | 0.0017     | 22              | 34.1706  | -1.681       | -13.8   | 23.2456    | 0.0        | -1.194         |


### Output Column Descriptions
- `name`: Name of input sequence.
- `anti-Shine-Dalgarno`: Anti-Shine-Dalgarno sequence.
- `start_codon`: Start codon for predicted translation initiation rate.
- `start_position`: Nucleotide position of first start codon base in input sequence. (1-indexed).
- `expression`: Predicted translation initation rate at this start codon.
- `RBS_distance_bp`: Number of nucleotides between the predicted ribosome-binding site (RBS) and the start codon.
- `dG_total`: Total change in free energy (ΔG) for translation initiation at this RBS.
- `dG_rRNA:mRNA`: Free energy term for ribosome binding to the mRNA.
- `dG_mRNA`: Free energy term for unfolding mRNA secondary structures that overlap the RBS and start codon region.
- `dG_spacing`: Free energy term that accounts for the effect of spacing (the number of nucleotides) between an RBS and start codon.
- `dG_standby`: Free energy term for unfolding mRNA secondary structures that occlude the standby site (the four bases upstream of the RBS).
- `dG_start_codon`: Free energy term for initiator tRNA binding to the start codons.

## Python Module Usage

OSTIR can also be called from within a user's Python script via the `run_ostir` function. This function returns a list of
the translation initiation rates (expression levels) predicted for each start codon in the sequence.

Example usage:
```python3
from ostir.ostir import run_ostir

seq = "ACUUCUAAUUUAUUCUAUUUAUUCGCGGAUAUGCAUAGGAGUGCUUCGAUGUCAU"
results = run_ostir(seq, name="my_sequence", aSD="ACCTCCTTA", threads=8)
print(results)

results = run_ostir(seq, start=31, end=31, aSD="ACCCCCTTA", verbose=True)
print(results)
```

Parameters:
- `seq`: mRNA sequence to search for translation initiation sites. (REQUIRED)
- `start`: Most 5' position to consider a start codon beginning. (1-indexed)
- `end`: "Most 3' position to consider a start codon beginning. (1-indexed)"
- `name`: Name or id of the sequence to include in output.
- `aSD`: anti-Shine-Dalgarno sequence. Defaults to that of *E. coli*.
- `threads`: Number of parallel processes to launch during prediction.
- `verbose`: Prints debug information
