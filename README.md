# OSTIR (Open Source Translation Initiation Rates)
[![Status](https://github.com/barricklab/ostir/actions/workflows/package_and_test.yml/badge.svg)](https://github.com/barricklab/ostir/actions/workflows/package_and_test.yml) [![status](https://joss.theoj.org/papers/43d02b32408a161e608b886f63e753c1/status.svg)](https://joss.theoj.org/papers/43d02b32408a161e608b886f63e753c1)


`OSTIR` is a
Python package for predicting the rates at which ribosomes will bind to and initiate
translation from different start codons in bacterial mRNAs. It uses the ViennaRNA Package to perform 
the necessary free energy calculations. The code builds on the last open source version of the
[RBS calculator](https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0).

`OSTIR` includes several improvements in usability. It supports multi-FASTA
input with command line parameters or CSV input that can define
parameters on a per-sequence basis. Additionally, `OSTIR` supports multi-threaded
execution, accelerating the analysis of very large sequences.

### [Please see the OSTIR Wiki for full documentation](https://github.com/barricklab/ostir/wiki)

# Quickstart

## Installation
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ostir/README.html)

`OSTIR` is a Python module and associated command line script. We recommend installing `OSTIR` using [Bioconda](https://bioconda.github.io/) on Linux or macOS. This will automatically install `OSTIR` and all of its dependencies, including [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) and the required Python modules.

From Bioconda (recommended; Linux, macOS):
- Run `conda install -c bioconda ostir`

From Pip (for experts; Linux, macOS, Windows):
- Download and install ViennaRNA, following the instructions [here](https://www.tbi.univie.ac.at/RNA/).
- Run `pip install ostir`

For information on installing for development see the [Wiki Documentation](https://github.com/barricklab/ostir/wiki/Installation).

## Command Line Usage

Print OSTIR help:
```
ostir -h
```

Run OSTIR on a sequence provided at the command line and print output to the console:
```
ostir -i TTCTAGATGAGAATAAGGTTATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGT
```

Run OSTIR on all sequences provided in a FASTA file and print output to a CSV file:
```
ostir -i input.fasta -o output.csv
```

More options and examples are described in the [Wiki Documentation](https://github.com/barricklab/ostir/wiki/Command-Line-Usage).

## Python Module Usage

Run OSTIR on a sequence inside of a Python script:

```python3
from ostir import run_ostir

seq = "ACUUCUAAUUUAUUCUAUUUAUUCGCGGAUAUGCAUAGGAGUGCUUCGAUGUCAU"
results = run_ostir(seq, name="my_sequence", threads=8)
print(results)
```

More options and examples are described in the [Wiki Documentation](https://github.com/barricklab/ostir/wiki/Python-Module-Usage).

