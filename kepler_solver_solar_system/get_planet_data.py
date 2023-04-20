import rebound
import numpy as np


data_r = np.loadtxt("radii")
sim = rebound.Simulation()
sim.units = ('yr','AU','Msun')
sim.add(["Sun","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

m = []
a = []
e = []

for i in range(1,8):
    orbit = sim.particles[i].calculate_orbit(primary=sim.particles[0])
    m.append(sim.particles[i].m)
    a.append(orbit.a)
    e.append(orbit.e)


np.savetxt("planet_data",np.c_[m,a,e,data_r[1:8]])

#for i in range(0,len(m)):

