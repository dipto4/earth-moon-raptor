import numpy as np
from ctypes import *
import os
import rebound

script_dir = os.path.abspath(os.path.dirname(__file__))
lib_path = os.path.join(script_dir, "twobodyPy.so")

libkepler = cdll.LoadLibrary(lib_path)

libkepler.integrate2body.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                c_double, POINTER(c_double), POINTER(c_double)]

m = [1.0,3.0036060914862233e-6,3.695250762067745e-8]
mtot = m[0] + m[1] + m[2]
xs = [0.0, 0.0, 0.0]
xe = [-7.4922108555337696e-01 ,6.8064940902596427e-01 ,0.0000000000000000e+00]
xm = [-7.5182886856733511e-01, 6.8137668169823307e-01,0.0]
xc = [6.2773119385896337e+08, -3.1165291138762122e+08, 7.1332041263377357e+08]

vxs = [0.0,0.0,0.0]
vxe = [-4.2254777889689503e+00, -4.5462293821324415e+00, 0.0000000000000000e+00]
vxm = [-4.2835996742299418e+00, -4.7427589358312829e+00, 0.0000000000000000e+00]
vxc = [-1.3243256423883512e-01, 6.5729268752177161e-02, -1.5047473509394596e-01]

xb = []
vxb = []
xci = [0.0,0.0,0.0]
vxci = [0.0,0.0,0.0]
for i in range(0,3):
    xbi = (m[0]*xs[i] + m[1] * xe[i] + m[2] * xm[i])/mtot
    vxbi = (m[0]*vxs[i] + m[1] * vxe[i] + m[2] * vxm[i])/mtot
    xb.append(xbi)
    vxb.append(vxbi)
#run the rebound simulation first
sim = rebound.Simulation()
sim.units = ('yr', 'AU', 'Msun')

sim.add(m=mtot,x=xb[0],y=xb[1],z=xb[2],vx=vxb[0],vy=vxb[1],vz=vxb[2])
sim.add(m=0.0,x=xc[0],y=xc[1],z=xc[2],vx=vxc[0],vy=vxc[1],vz=vxc[2])

sim.move_to_com()
orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
#print(orbit.a,orbit.e)

#print(orbit.T)
t_end = orbit.T-20
#print(t_end)
sim.integrate(t_end)

#now run the kepler solution
xb_c = (c_double * 3)(*xb)
xc_c = (c_double * 3)(*xc)
xci_c = (c_double *3)(*xci)

vxb_c = (c_double * 3)(*vxb)
vxc_c = (c_double * 3)(*vxc)
vxci_c = (c_double *3)(*vxci)

coll_flag = c_int(0)

libkepler.integrate2body(xb_c, xc_c, vxb_c, vxc_c, c_double(t_end), xci_c, vxci_c)
bary = sim.particles[0]
comet = sim.particles[1]

dist=np.sqrt((xci_c[0])**2 +(xci_c[1])**2 + (xci_c[2])**2)
dist2=np.sqrt((comet.x)**2 + (comet.y)**2 + (comet.z)**2)
print(dist,dist2)

print("------------posvel-info-----------")
print("{:.16e} {:.16e} {:.16e}".format(comet.x,comet.y,comet.z))
print("{:.16e} {:.16e} {:.16e}".format(xci_c[0],xci_c[1],xci_c[2]))

print("{:.16e} {:.16e} {:.16e}".format(comet.vx,comet.vy,comet.vz))
print("{:.16e} {:.16e} {:.16e}".format(vxci_c[0],vxci_c[1],vxci_c[2]))

print(bary.x,bary.y,bary.z)
print(bary.vx,bary.vy,bary.vz)
print(xb_c[0],xb_c[1],xb_c[2])
print(vxb_c[0],vxb_c[1],vxb_c[2])
print("-----debug info--------")

orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
print(orbit.a,orbit.e)

sim2 = rebound.Simulation()
sim2.units = ('yr', 'AU', 'Msun')
sim2.add(m=mtot,x=xb_c[0],y=xb_c[1],z=xb_c[2],vx=vxb_c[0],vy=vxb_c[1],vz=vxb_c[2])
sim2.add(m=0.0,x=xci_c[0],y=xci_c[1],z=xci_c[2],vx=vxci_c[0],vy=vxci_c[1],vz=vxci_c[2])
orbit = sim2.particles[1].calculate_orbit(primary=sim2.particles[0])
print(orbit.a,orbit.e)
#print(comet.x,comet.y,comet.z)
#print(xci_c[0],xci_c[1],xci_c[2])

#print(comet.vx,comet.vy,comet.vz)
#print(vxci_c[0],vxci_c[1],vxci_c[2])

#for i in xc_c:
#    print(i)

#for i in vxc_c:
#    print(i)
