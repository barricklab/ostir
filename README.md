# OSTIR

## Open Source Translation Initiation Rates

`OSTIR` (Open Source Translation Initiation of Ribosomes) is an open source
Python package that integrates RNA spacing considerations and the
ViennaRNA software suite to perform free energy calculations
for binding events involved during protein translation initiation. This work is
derived from the
[last open source version of the calculator](https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0)
described in [Salis 2009](https://doi.org/10.1038/nbt.1568).

`OSTIR` was written with a user-friendliness and high throughput in mind.
Unlike related tools, this work allows for user installation through Bioconda or PyPi (requires installing dependencies),
requiring no commercial licensing or export to a webserver. `OSTIR` supports multi-FASTA
input with command line parameters or CSV input with support for defining
parameters on a per-sample basis. Additionally, `OSTIR` supports asynchronous
processing, accelerating use cases that require screening very large mRNAs.

### Installation:

This package was tested on Ubuntu 20.04.1 LTS (under WSL). Installation steps and support may vary depending on platform. It is recommended
that you run `ostir --validate` after installation to ensure dependencies are properly installed and that a consistency test is passed.

From Conda (coming soon):
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

OSTIR is executable as a command-line python script, `ostir`.

Flags:
- `-h, --help`: Shows a help message
- `-i str/filepath`, --input str/filepath: Defines the input sequence, FASTA filepath, or CSV filepath.
- `[-o filepath, --output filepath]`: Defines the output filepath (csv formatted). If not provided, results will output
  to console
- `[-v, --verbose]`: Prints additional run information to the console
- `[-s int, --start int]`: Defines the most 5' position to consider start codons. Defaults to the first base.
- `[-e int, --end int]`: Defines the most 3' position to consider start codons. If `-s` is set, this defaults to 3 bases
  upstream, otherwise defaults to the end of the sequence
- `[-r str, --rRNA str]`: Defines the rRNA anti-Shine Dalgarno sequence. Defaults to that of E. coli
- `[j int, --threads int]'`: Defines how many mRNAs will be analyzed in parallel at a time
- `[t str, --type str]'`: Defines input type (seq|fasta|csv). If not provided, OSTIR will attempt to autodetect.

CSV-based inputs support overriding command-line flags `start`, `end`, and `sd` (Anti Shine Dalgarno) through csv
columns. The column `seq` is required. `name` is also supported, which simply returns itself in the output.

OSTIR can also be imported into a python script. Most users will want to call the `run_ostir` function.
Running OSTIR this way is pythonic and both expects and returns zero indexed start/end positions unless otherwise defined.

```python3
from ostir.ostir import run_ostir
run_ostir(seq, outfile, start_loc, end_loc, name, sd, threads, verbose, posindex)
```

Arguments:
- `seq`: Sequence to calculate binding energies for
- `outfile`: Filepath for output csv
- `start_loc`: First base to start considering start codons. Defaults to first base
- `end_loc`: Last base to start considering start codons. Defaults to end of sequence
- `name`: Returns itself, useful for tagging things for downstream processing.
- `sd`: Defines anti-Shine-Dalgarno sequence. Defaults to that of E. coli's
- `threads`: Defines parallel processing workers, roughly equivalent to multithreading cores
- `verbose`: Prints debug information
- `posindex`: Determines indexing for input/return start/end positions. Generally 0 or 1.
