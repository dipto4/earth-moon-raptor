import numpy as np
from ctypes import *
import os


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
                                POINTER(c_double), c_int, c_double]
def cart2jacobi(mass,x,y,z,vx,vy,vz):
    N = len(mass)
    mass_c = (c_double * N) (*mass)
    x_c = (c_double * N) (*x)
    y_c = (c_double * N) (*y)
    z_c = (c_double * N) (*z)
    vx_c = (c_double * N) (*vx)
    vy_c = (c_double * N) (*vy)
    vz_c = (c_double * N) (*vz)


    eta = (c_double*N) (*np.zeros(N))
    x_star = (c_double *N) (*np.zeros(N))
    y_star = (c_double *N) (*np.zeros(N))
    z_star = (c_double *N) (*np.zeros(N))
    vx_star = (c_double *N) (*np.zeros(N))
    vy_star = (c_double *N) (*np.zeros(N))
    vz_star = (c_double *N) (*np.zeros(N))

    libkepler.Py_cart2jacobi(mass_c, x_c, y_c,z_c,
                          vx_c,vy_c,vz_c,c_int(N),
                          eta,x_star,y_star,z_star,
                          vx_star,vy_star,vz_star)

    xx_star = np.array((x_star,y_star,z_star))
    xx_star = xx_star.T
    vv_star = np.array((vx_star,vy_star,vz_star))
    vv_star = vv_star.T
    return np.array(eta), xx_star, vv_star


def jacobi2cart(eta,x_star,v_star,mass):
    N = len(mass)
    mass_c = (c_double * N) (*mass)
    x_c = (c_double *N) (*np.zeros(N))
    y_c = (c_double *N) (*np.zeros(N))
    z_c = (c_double *N) (*np.zeros(N))
    vx_c = (c_double *N) (*np.zeros(N))
    vy_c = (c_double *N) (*np.zeros(N))
    vz_c = (c_double *N) (*np.zeros(N))

    eta_c = (c_double * N) (*eta)
    x_star_c = (c_double * N) (*x_star[:,0])
    y_star_c = (c_double * N) (*x_star[:,1])
    z_star_c = (c_double * N) (*x_star[:,2])

    vx_star_c = (c_double * N) (*v_star[:,0])
    vy_star_c = (c_double * N) (*v_star[:,1])
    vz_star_c = (c_double * N) (*v_star[:,2])

    libkepler.Py_jacobi2cart(mass_c, x_c, y_c,z_c,
                          vx_c,vy_c,vz_c,c_int(N),
                          eta_c,x_star_c,y_star_c,z_star_c,
                          vx_star_c,vy_star_c,vz_star_c)

    xyz = np.array((x_c,y_c,z_c))
    xyz = xyz.T
    vel = np.array((vx_c,vy_c,vz_c))
    vel = vel.T

    return xyz,vel
    pass

def integrate(x,y,z,vx,vy,vz,m,t_end):
    N = len(x)
    x_c = (c_double * N) (*x)
    y_c = (c_double * N) (*y)
    z_c = (c_double * N) (*z)
    vx_c = (c_double * N) (*vx)
    vy_c = (c_double * N) (*vy)
    vz_c = (c_double * N) (*vz)
    m_c = (c_double * N) (*m)
    libkepler.integrate(x_c,y_c,z_c,vx_c,vy_c,vz_c,m_c,c_int(N),c_double(t_end))

if __name__ == "__main__":
    #for testing purposes
    #mass = np.array([5.,4.,3.,2.])
    #x = np.array([0.,1.,2.,4.])
    #y = np.array([0.,1.,-1.,3.])
    #z = np.array([0.,0.,0.,0.])

    #vx = np.array([0.,1.,-1.,-1.])
    #vy = np.array([0.,2.,0.,2.])
    #vz = np.array([0.,0.,0.,0.])
    data = np.loadtxt("solar_system.dat")
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

    integrate(x,y,z,vx,vy,vz,mass,200)

    xyz = np.array((x,y,z))

    xyz = xyz.T
    vel = np.array((vx,vy,vz))
    vel = vel.T
    #print(xyz)
    #print(vel)
    #eta, x_star, v_star= cart2jacobi(mass,x,y,z,vx,vy,vz)
    #print(eta)
    #print(x_star)
    #print(v_star)

    #xyz_r, vel_r = jacobi2cart(eta,x_star,v_star,mass)
    #print(xyz_r-xyz)
    #print(vel_r-vel)
