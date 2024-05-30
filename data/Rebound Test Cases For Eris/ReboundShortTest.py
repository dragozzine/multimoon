import numpy as np
import rebound
import matplotlib.pyplot as plt
from schwimmbad import MPIPool
from matplotlib.colors import LogNorm


'''
import numpy as np
import rebound
import matplotlib.pyplot as plt

ErisData = [1.6, 0.1163]

run_num = 100
sims = np.zeros((run_num,run_num))
color_map = np.zeros((run_num,run_num))
a_array = np.linspace(1.0,3.0,num=run_num, endpoint=True)
m_array = np.geomspace(1e-4,1e-2,num=run_num, endpoint=True)
for j in range(run_num):
    for k in range(run_num):
        sim = rebound.Simulation()
        sim.add(m=ErisData[0])
        a_num = a_array[j]
        m_num = m_array[k]
        sim.add(m = m_num, r = 0.0260, a = a_num, e = 0.243, inc = 33.961, Omega = 134.664, omega = 114, M = 210)
        sim.add(m = 2.65e-2, r = 0.0350, a = 3.7250, e = 0.006, inc = 62.097, Omega = 139.16, omega = 200, M = 137.277)

        sim.move_to_com()

np.savetxt('color.txt', color_map)
'''

ErisData = [1.616, 0.1163]
global sims
global color_map
def loop1(j):
    run_num=100
    if j%5==0:
        print(j)
    sims = np.zeros(run_num)
    color_map = np.zeros(run_num)
    a_array = np.linspace(1.0,3.0,num=run_num, endpoint=True)
    m_array = np.geomspace(1e-4,1.6e-2,num=run_num, endpoint=True)
    for k in range(run_num):
        sim = rebound.Simulation()
        sim.add(m=ErisData[0])
        a_num = a_array[j]
        m_num = m_array[k]
        sim.add(m = m_num, r = 0.0260, a = a_num, e = 0.243, inc = 33.961/180*np.pi, Omega = 134.664/180*np.pi, omega = 114/180*np.pi, M = 210/180*np.pi)
        sim.add(m = 2.65e-2, r = 0.0350, a = 3.7250, e = 0.006, inc = 62.097/180*np.pi, Omega = 139.16/180*np.pi, omega = 200/180*np.pi, M = 137.277/180*np.pi)

        sim.move_to_com()
#sim.status()
        rebound.OrbitPlot(sim)
        Noutputs = 10

        year = 2.*np.pi # One year in units where G=1
        times = np.linspace(0.,10**4*year, Noutputs)
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
            
        dist_arr = np.sqrt(x[0]**2+y[0]**2)
        max_dist = np.max(dist_arr)
        min_dist = np.min(dist_arr)
        
        delta_a = (max_dist-a_num)/a_num #- (a_num-min_dist)/a_num
        color_map[k] = delta_a
        
        if np.sqrt(x[0][-1]**2+y[0][-1]**2) > 20.0:
            sims[k] = False
        elif np.sqrt(x[1][-1]**2+y[1][-1]**2) > 20.0:
            sims[k] = False    
        else:
            sims[k] = True
        #print('k j ', k, j)

    return sims, color_map
        
with MPIPool() as pool:   
    run_num = 100
    data = pool.map(loop1, range(run_num))
    sims = np.zeros((run_num,run_num))
    color_map = np.zeros((run_num,run_num))
    for i in range(len(data)):
        sims[i] = data[i][0]
        color_map[i] = data[i][1]
    
    sims = np.transpose(sims)
    color_map = np.transpose(color_map)
    #print(sims, color_map)
    

    plt.imshow(np.flip(sims,axis=0), cmap='binary',extent=[1.0,3.0,18,20])
    plt.xlabel('Kilometers (in 10,000)')
    plt.ylabel('Mass (1.6e_) kg')
    plt.savefig('sims_binary_chart_short.pdf')
    
    plt.imshow(np.flip(color_map,axis=0), extent=[1.0,3.0,18,20], norm=LogNorm())
    plt.xlabel('Kilometers (in 10,000)')
    plt.ylabel('Mass (1.6e_) kg')
    plt.title('Delta_a')
    plt.colorbar()
    plt.show()

    plt.savefig('delta_a_color_short.pdf')

    np.savetxt('sims_long_short.txt',sims)
    np.savetxt('delta_a_arr_short.txt', color_map)