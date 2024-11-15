import os
import subprocess
import datetime
import asyncio
import math
import os
import shutil

from ostir.ostir import parse_fasta

import psutil

from tempfile import mkdtemp
from rich import print
from progress.spinner import Spinner
from progress.bar import Bar
from time import sleep
from contextlib import contextmanager

conda_namespace = "benchmark_ostir"
conda_path = "/usr/bin/micromamba"

starting_dir = os.path.dirname(os.path.abspath(__file__))

# importing librarie

def get_memory_usage(pid):
    try:
        process = psutil.Process(pid)
        memory = process.memory_info().rss
        children = process.children(recursive=True)
        for child in children:
            try:
                memory += child.memory_info().rss
            except psutil.NoSuchProcess:
                pass
        return memory
    except psutil.NoSuchProcess:
        return 0

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join([complement[base] for base in seq[::-1]])

@contextmanager
def rhobust_tempdir():
    temp_dir = mkdtemp()
    try:
        yield temp_dir
    finally:
        current_dir = os.getcwd()
        assert not current_dir.startswith(temp_dir ), f"Current directory '{current_dir}' is within '{temp_dir}' or its subdirectories."
        shutil.rmtree(temp_dir)
        assert not os.path.exists(temp_dir), f"Temporary directory '{temp_dir}' failed to be removed."


class Benchmarker():
    def __init__(self):
        pass

    def benchmark_t7(self, n=1):
        print("Benchmarking T7 genome")
        t7_genome = os.path.abspath("T7_genome.fasta")
        self.benchmark(t7_genome, n)

    def benchmark_mg1655(self, n=1):
        print("Benchmarking MG1655 genome")
        mg1655_genome = os.path.abspath("MG1655.fna")
        self.benchmark(mg1655_genome, n)

    def benchmark_ecoli(self):
        pass

    def execute(self):
        starttime = datetime.datetime.now()
        max_memory_usage = 0
        for command in self.commands:
            p = subprocess.Popen(command, shell=False, env=self.env,
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            poll = p.poll()
            while poll is None:
                poll = p.poll()
                max_memory_usage = max(max_memory_usage, get_memory_usage(p.pid))
                sleep(0.001)
        endtime = datetime.datetime.now()
        time = (endtime - starttime).total_seconds()
        return time, max_memory_usage

class BenchmarkerRBSCalc(Benchmarker):

    def benchmark(self, input_path, no=1):
        print(f"Setting up RBS Calculator 1.0")

        with rhobust_tempdir() as temp_dir:

            # Setup
            environment = os.environ.copy()
            nupack_bin_path = os.path.abspath("nupack/bin")
            nupack_base_path = os.path.abspath("nupack")
            environment["PATH"] = f"{nupack_bin_path}:{nupack_base_path}:{environment['PATH']}"
            environment["NUPACKHOME"] = nupack_base_path
            self.env = environment

            patch_path = os.path.abspath("rbscalc_bugfixes.patch")

            os.chdir(temp_dir)
            url = "https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0.git"
            subprocess.run(f"git clone {url} --quiet", shell=True)
            os.chdir("Ribosome-Binding-Site-Calculator-v1.0")
            subprocess.run(f"git apply {patch_path} --quiet", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            input_string = parse_fasta(input_path)[0][1]
            reverse_complement_string = reverse_complement(input_string)

            # write rc to fasta
            with open(f"{temp_dir}/rc.fasta", "w") as f:
                f.write(f">{input_path}\n{reverse_complement_string}")


            self.commands = [f"micromamba run -n {conda_namespace} python {temp_dir}/Ribosome-Binding-Site-Calculator-v1.0/Run_RBS_Calculator.py {input_path}".split(),
                f"micromamba run -n {conda_namespace} python {temp_dir}/Ribosome-Binding-Site-Calculator-v1.0/Run_RBS_Calculator.py {temp_dir}/rc.fasta".split()]

            # Run benchmark
            times = []
            mem = []
            with Bar('Running benchmark', max=no) as bar:
                for _ in range(no):
                    with rhobust_tempdir() as temp_dir2:
                        os.chdir(temp_dir2)
                        result = self.execute()
                        times.append(result[0])
                        mem.append(result[1])
                        bar.next()
                        os.chdir(starting_dir)
            bar.finish()

            # Calculate averages and standard deviations
            avg_time = sum(times) / len(times)
            std_time = math.sqrt(sum([(x - avg_time) ** 2 for x in times]) / len(times))
            avg_mem = sum(mem) / len(mem)
            std_mem = math.sqrt(sum([(x - avg_mem) ** 2 for x in mem]) / len(mem))

            print(f"Average time: {avg_time} +- {std_time} seconds")
            print(f"Average memory: {avg_mem / 1024 / 1024} +- {std_mem / 1024 / 1024} MB")
            os.chdir(starting_dir)

class BenchmarkerOstir(Benchmarker):
    def __init__(self, commit="latest"):

        self.commit = commit
        self.env = None
        super().__init__()

    def benchmark(self, input_path, no=1):
        print(f"Setting up OSTIR ({self.commit})")
        self.commands = [f'{conda_path} run -n {conda_namespace} ostir -i {input_path} -o /dev/null -t fasta -j 8'.split()]
        with rhobust_tempdir() as temp_dir:

            # Setup
            os.chdir(temp_dir)
            url = "https://github.com/barricklab/ostir.git"
            subprocess.run(f"git clone {url} --quiet", shell=True)
            os.chdir("ostir")
            if self.commit != "latest":
                subprocess.run(f"git checkout {self.commit} --quiet", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            subprocess.run(f"micromamba env update -n {conda_namespace} -f environment.yml -q -y", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run(f"micromamba run -n {conda_namespace} pip install . --quiet", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Run unittests

            # Run benchmark
            times = []
            mem = []
            with Bar('Running benchmark', max=no) as bar:
                for _ in range(no):
                    with rhobust_tempdir() as temp_dir2:
                        os.chdir(temp_dir2)
                        result = self.execute()
                        times.append(result[0])
                        mem.append(result[1])
                        bar.next()
                        os.chdir(starting_dir)
            bar.finish()

            # Calculate averages and standard deviations
            avg_time = sum(times) / len(times)
            std_time = math.sqrt(sum([(x - avg_time) ** 2 for x in times]) / len(times))
            avg_mem = sum(mem) / len(mem)
            std_mem = math.sqrt(sum([(x - avg_mem) ** 2 for x in mem]) / len(mem))

            print(f"Average time: {avg_time} +- {std_time} seconds")
            print(f"Average memory: {avg_mem / 1024 / 1024} +- {std_mem / 1024 / 1024} MB")
            os.chdir(starting_dir)



def setup_micromamba():
    print("Setting up Micromamba")
    # delete existing environment
    subprocess.run(f"micromamba env remove -n {conda_namespace} -q -y",
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # create new environment
    subprocess.run(f"micromamba create -n {conda_namespace} -q -y",
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

if __name__ == "__main__":
    no_t7 = 100
    no_mg1655 = 3
    print("----------")
    os.chdir(starting_dir)
    setup_micromamba()
    ostir_benchmark = BenchmarkerOstir()
    ostir_benchmark.benchmark_t7(no_t7)
    #os.chdir(starting_dir)
    ostir_benchmark.benchmark_mg1655(no_mg1655)

    print("----------\n")
    os.chdir(starting_dir)
    setup_micromamba()
    subprocess.run(f"micromamba install -n {conda_namespace} viennarna=2.4.18 -y -c bioconda",
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ostir_benchmark = BenchmarkerOstir(commit = "79eb4b9")
    ostir_benchmark.benchmark_t7(100)
    os.chdir(starting_dir)
    ostir_benchmark.benchmark_mg1655(no_t7)
    print("----------\n")
    os.chdir(starting_dir)
    setup_micromamba()
    subprocess.run(f"micromamba install -n {conda_namespace} python=2 -y",
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rbscalc_benchmark = BenchmarkerRBSCalc()
    rbscalc_benchmark.benchmark_t7(no_t7)
    os.chdir(starting_dir)
    rbscalc_benchmark.benchmark_mg1655(no_mg1655)

    exit()
