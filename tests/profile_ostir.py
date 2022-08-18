import pstats, cProfile

from ostir.ostir import run_ostir

def test_ostir(in_seq="TTCTAGATGAGAATAAGGTTATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGT", time=1000):
    for i in range(time):
        run_ostir(in_seq)


cProfile.runctx('test_ostir()', globals(), locals(), 'output/ostir.profile')
s = pstats.Stats('output/ostir.profile')
s.strip_dirs().sort_stats('time').print_stats()