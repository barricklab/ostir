# OSTIR

## Open Source Translation Initiation Rates

`OSTIR` (Open Source Translation Initiation Rates) is a
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

### Installation

From Conda:
- Run `conda install --bioconda ostir`

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

### Usage

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

Arguments:
- `seq`: mRNA sequence to search for translation initiation sites (REQUIRED)
- `start`: Most 5' position to consider a start codon beginning (1-indexed)
- `end`: "Most 3' position to consider a start codon beginning (1-indexed)"
- `name`: Name or id of the sequence to include in output.
- `aSD`: anti-Shine-Dalgarno sequence. Defaults to that of E. coli.
- `threads`: Number of parallel processes to launch during prediction.
- `verbose`: Prints debug information
