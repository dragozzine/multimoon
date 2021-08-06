import numpy as np
import rebound
import matplotlib.pyplot as plt
from schwimmbad import MPIPool

ErisData = [1.6, 0.1163]
global sims
def loop1(j):
    run_num=100
    sims = np.zeros(run_num)
    a_array = np.linspace(1.0,3.0,num=run_num, endpoint=True)
    m_array = np.geomspace(1.6e-4,1.6e-2,num=run_num, endpoint=True)
    for k in range(run_num):
        sim = rebound.Simulation()
        sim.add(m=ErisData[0])
        a_num = a_array[j]
        m_num = m_array[k]
        sim.add(m = m_num, r = 0.0100, a = a_num, e = 0.01489, inc = 48, Omega = 50, omega = 131.89, M = 50)
        sim.add(m = 1.6e-2, r = 0.0450, a = 3.7211, e = 0.01, inc = 61.25, Omega = 260, omega = 139, M = 80)

        sim.move_to_com()
#sim.status()
        rebound.OrbitPlot(sim)
        Noutputs = 1000

        year = 2.*np.pi # One year in units where G=1
        times = np.linspace(0.,100000.*year, Noutputs)
        x = np.zeros((2,Noutputs))
        y = np.zeros((2,Noutputs))

        sim.integrator = "ias15" # IAS15 is the default integrator, so we actually don't need this line
        sim.move_to_com()        # We always move to the center of momentum frame before an integration
        ps = sim.particles       # ps is now an array of pointers and will change as the simulation runs

        for i,time in enumerate(times):
            sim.integrate(time)
            x[0][i] = ps[1].x   # This stores the data which allows us to plot it later
            y[0][i] = ps[1].y
            x[1][i] = ps[2].x
            y[1][i] = ps[2].y
    
        if np.sqrt(x[0][-1]**2+y[0][-1]**2) > 20.0:
            sims[k] = False
        elif np.sqrt(x[1][-1]**2+y[1][-1]**2) > 20.0:
            sims[k] = False
        else:
            sims[k] = True
        print(sims[k], k, j)
    return sims
        
with MPIPool() as pool:   
    sims = pool.map(loop1, range(100))
    sims = np.transpose(sims)
    print(sims)
    

    plt.imshow(sims, cmap='binary',extent=[1.0,3.0,18,20])
    plt.xlabel('Kilometers (in 10,000)')
    plt.ylabel('Mass (1.6e_) kg')

    plt.savefig('sims_binary_chart.png')

    np.savetxt('sims_long.txt',sims)