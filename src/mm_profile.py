# OR
# python -m cProfile -o output.pstats mm_run.py
# python gprof2dot.py -f pstats output.pstats | dot -Tpng -o profile.png
# then look at profile.png

import cProfile
import pstats
import io


pr = cProfile.Profile()
pr.enable()

if __name__ == '__main__':
    import mm_run
    mm_run.run()

pr.disable()
import os, glob
folder = max(glob.glob(os.path.join('../runs/', '*/')), key=os.path.getmtime)
filename = folder+'profile.prof'
pr.dump_stats(filename)
print('profile.prof dumped to ', filename)

#pr.dump_stats('profile1.prof')

s1 = io.StringIO()
ps1 = pstats.Stats(pr, stream=s1).sort_stats('cumtime')
ps1.print_stats()

#with open('test_cum.txt', 'w+') as f:
 #   f.write(s1.getvalue())

s2 = io.StringIO()    
ps2 = pstats.Stats(pr, stream=s2).sort_stats('tottime')
ps2.print_stats()

#with open('test_tot.txt', 'w+') as f:
 #   f.write(s2.getvalue())
