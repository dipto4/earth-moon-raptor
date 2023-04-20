import numpy as np
from ctypes import *
import os
import time

libkepler = cdll.LoadLibrary('./solarsystemPy.so')
libkepler.Py_cart2jacobi.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                     POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int,
                                     POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                     POINTER(c_double), POINTER(c_double), POINTER(c_double)]

libkepler.Py_jacobi2cart.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                     POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int,
                                     POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                     POINTER(c_double), POINTER(c_double), POINTER(c_double)]

libkepler.integrate.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                POINTER(c_double), POINTER(c_double),c_int, c_double,POINTER(c_int)]

def integrate(x,y,z,vx,vy,vz,m,r,t_end):
    N = len(x)
    x_c = (c_double * N) (*x)
    y_c = (c_double * N) (*y)
    z_c = (c_double * N) (*z)
    vx_c = (c_double * N) (*vx)
    vy_c = (c_double * N) (*vy)
    vz_c = (c_double * N) (*vz)
    m_c = (c_double * N) (*m)
    r_c = (c_double * N) (*r)
    collb_flag = c_int(0)
    libkepler.integrate(x_c,y_c,z_c,vx_c,vy_c,vz_c,m_c,r_c,c_int(N),c_double(t_end),byref(collb_flag))
    for i in range(0,N):
        x[i] = x_c[i]
        y[i] = y_c[i]
        z[i] = z_c[i]
        vx[i] = vx_c[i]
        vy[i] = vy_c[i]
        vz[i] = vz_c[i]
    return collb_flag
    #print(x_c[-1],y_c[-1],z_c[-1])
if __name__ == "__main__":
    #for testing purposes
    #mass = np.array([5.,4.,3.,2.])
    #x = np.array([0.,1.,2.,4.])
    #y = np.array([0.,1.,-1.,3.])
    #z = np.array([0.,0.,0.,0.])

    #vx = np.array([0.,1.,-1.,-1.])
    #vy = np.array([0.,2.,0.,2.])
    #vz = np.array([0.,0.,0.,0.])
    data = np.loadtxt("solar_system_starman.dat")
    r = np.loadtxt("radii")
    mass = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    vx = data[:,4]
    vy = data[:,5]
    vz = data[:,6]
    x -= x[0]
    y -= y[0]
    z -= z[0]
    vx -= vx[0]
    vy -= vy[0]
    vz -= vz[0]
    start_time = time.time()
    for i in range(0,100):
        integrate(x,y,z,vx,vy,vz,mass,r,100)
    print(time.time()-start_time)

