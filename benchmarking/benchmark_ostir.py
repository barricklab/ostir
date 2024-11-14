import os
import subprocess
import datetime
import asyncio
import math
import os

import psutil

from tempfile import TemporaryDirectory
from rich import print
from progress.spinner import Spinner
from progress.bar import Bar

conda_namespace = "benchmark_ostir"
conda_path = "/usr/bin/micromamba"

# importing librarie

def get_memory_usage(pid):
    try:
        process = psutil.Process(pid)
        return process.memory_info().rss
    except:
        return 0


class Benchmarker():
    def __init__(self):
        pass

    def benchmark_t7(self):
        t7_genome = os.path.abspath("T7_genome.fasta")
        self.benchmark(t7_genome, 5)

    def benchmark_ecoli(self):
        pass

class BenchmarkerRBSCalc(Benchmarker):

    def ensure_installed(self):
        pass


class BenchmarkerOstir(Benchmarker):
    def __init__(self, commit="latest"):

        self.commit = commit
        super().__init__()

    def benchmark(self, input_path, no=1):
        print(f"Setting up OSTIR ({self.commit})")
        with TemporaryDirectory() as temp_dir:

            # Setup
            os.chdir(temp_dir)
            url = "https://github.com/barricklab/ostir.git"
            subprocess.run(f"git clone {url} --quiet", shell=True)
            os.chdir("ostir")
            if self.commit != "latest":
                subprocess.run(f"git checkout {self.commit} --quiet", shell=True)

            subprocess.run(f"micromamba env update -n {conda_namespace} -f environment.yml -q -y", shell=True)
            subprocess.run(f"micromamba run -n {conda_namespace} pip install . --quiet", shell=True)

            # Run unittests

            # Run benchmark
            times = []
            mem = []
            with Bar('Running benchmark', max=no) as bar:
                for _ in range(no):
                    try:
                        os.mkdir("output")
                    except FileExistsError:
                        pass
                    result = self.execute(input_path)
                    times.append(result[0])
                    mem.append(result[1])
                    bar.next()
            bar.finish()

            # Calculate averages and standard deviations
            avg_time = sum(times) / len(times)
            std_time = math.sqrt(sum([(x - avg_time) ** 2 for x in times]) / len(times))
            avg_mem = sum(mem) / len(mem)
            std_mem = math.sqrt(sum([(x - avg_mem) ** 2 for x in mem]) / len(mem))

            print(f"Average time: {avg_time} +- {std_time} seconds")
            print(f"Average memory: {avg_mem / 1024 / 1024} +- {std_mem / 1024 / 1024} MB")


    @staticmethod
    def execute(input_path):
        starttime = datetime.datetime.now()
        p = subprocess.Popen(f'{conda_path} run -n {conda_namespace} ostir -i {input_path} -o output/output -t fasta -j 8 -v 0'.split(), shell=False)
        poll = p.poll()
        max_memory_usage = 0
        while poll is None:
            poll = p.poll()
            max_memory_usage = max(max_memory_usage, get_memory_usage(p.pid))
        endtime = datetime.datetime.now()
        time = (endtime - starttime).total_seconds()
        return time, max_memory_usage


def setup_micromamba():
    print("Setting up Micromamba")
    # delete existing environment
    subprocess.run(f"micromamba env remove -n {conda_namespace} -q -y", shell=True)
    # create new environment
    subprocess.run(f"micromamba create -n {conda_namespace} -q -y", shell=True)

if __name__ == "__main__":
    starting_dir = os.getcwd()
    setup_micromamba()
    ostir_benchmark = BenchmarkerOstir()
    ostir_benchmark.benchmark_t7()

    exit()
