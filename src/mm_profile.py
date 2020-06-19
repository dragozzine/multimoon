import cProfile
import pstats
import io


pr = cProfile.Profile()
pr.enable()

import mm_run

pr.disable()
s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
ps.print_stats()

with open('test.txt', 'w+') as f:
    f.write(s.getvalue())
