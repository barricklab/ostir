from timeit import timeit

if __name__ == "__main__":
    setup = "from __main__ import my_function"
    t7_genome_string = ""
    with open('input/T7_genome.fasta', 'r', encoding='utf8') as f:
        for line in f:
            if line.startswith('>'):
                continue
            t7_genome_string += line.strip()
    setup = "from ostir import run_ostir"
    print(timeit(f"run_ostir('{t7_genome_string}')", setup=setup, number=100))
