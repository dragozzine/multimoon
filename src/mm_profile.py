import cProfile
import pstats
import io


pr = cProfile.Profile()
pr.enable()

import mm_run

pr.disable()
s1 = io.StringIO()
ps1 = pstats.Stats(pr, stream=s1).sort_stats('cumtime')
ps1.print_stats()

with open('test_cum.txt', 'w+') as f:
    f.write(s1.getvalue())

s2 = io.StringIO()    
ps2 = pstats.Stats(pr, stream=s2).sort_stats('tottime')
ps2.print_stats()

with open('test_tot.txt', 'w+') as f:
    f.write(s2.getvalue())
