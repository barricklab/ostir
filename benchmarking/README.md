This is a benchmarking script for OSTIR and related RBS calculators.

Successfully running the benchmark with the RBS Calculator 1.0 requires installing NUPACK. To do so, first somehow aquire the old 2.1 version of NUPACK and the rna1999.dG parameter file. Unpack NUPACK and name it nupack, then place it in the root of this directory. Put the rna1999.dG file in the nupack/parameters directory. You may have to compile NUPACK as well.

This benchmark handles insalling each calculator using Micromamba. The following calculators are tested:

- OSTIR (Latest)
- OSTIR (1.0.6, the version at time of publication)
- RBS Calculator 1.0

The benchmarking script is run with `python benchmark_ostir.py`. It requires `progress`, `rich`, `psutil`, and `ostir` to be installed in the executing python environment as well as `git`.
The MG1655 genome can be found [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_904425475.1/)
